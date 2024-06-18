#### Libraries ####
library(tidyverse)
library(magrittr)
library(Biostrings)
library(bioseq)
library(patchwork)

#### Settings ####
alpha <- 0.05

#### Functions ####
preprocess_codons <- function(x, pos){
  mutate(x, locus = str_c(Seqid, position, sep = '-')) %>%
    filter(locus %in% lfmm_results$locus)
}

find_snp_gene_overlap <- function(genes, snp){
  
  gene_range <- with(genes, split(IRanges(gene_start, gene_end), Seqid))
  snp_range <- with(snp, split(IRanges(position, width = 1), chromosome))
  
  the_overlap <- findOverlaps(gene_range, snp_range, type = 'any', select = 'all')
  
  as.list(the_overlap) %>%
    map_dfr(as_tibble, .id = 'Seqid') %>%
    dplyr::rename(gene_row = queryHits,
                  snp_row = subjectHits) %>%
    left_join(group_by(genes, Seqid) %>% mutate(gene_row = row_number()),
              by = c('Seqid', 'gene_row')) %>%
    full_join(group_by(snp, chromosome) %>% mutate(snp_row = row_number()),
              by = c('Seqid' = 'chromosome', 'snp_row'))%>%
    select(-gene_row, -snp_row) %>%
    dplyr::rename(chromosome = Seqid,
                  gene_strand = Strand) %>%
    arrange(chromosome, position) %>%
    relocate(starts_with('gene'), .after = last_col())
}

#### Read in Data ####
genome <- readBStringSet('../genome/acerv_genome.fasta')
names(genome) <- str_remove(names(genome), ' Acropora cervicornis isolate K2 Acerv_scaffold_[0-9]+, whole genome shotgun sequence')

lfmm_results <- read_csv('../Results/lfmm_results.csv', show_col_types = FALSE) %>%
  filter(pop_group == 'all') %>%
  select(-pop_group) %>%
  mutate(position = as.integer(position),
         locus = str_c(chromosome, position, sep = '-'), .after = 'position') 


snp_data <- read_delim('../variant_calling/genotyping/genolike.mafs.gz', delim = '\t', col_types = 'cicccnnnni') %>%
  select(chromo, position, major, minor, ref) %>%
  dplyr::rename(chromosome = chromo) %>%
  mutate(locus = str_c(chromosome, position, sep = '-'),
         alt = if_else(ref == major, minor, major))

gene_locations <- read_csv('../genome/genomic_gene_id_locations.csv', show_col_types = FALSE)
structural_annotations <- read_csv('../genome/genomic_structural_annotation.csv', show_col_types = FALSE)
functional_annotations <- filter(structural_annotations, structure_type == 'CDS') %>%
  select(Seqid, gene_id, structure_attributes) %>%
  distinct
  

codon_sequences <- read_csv_chunked('../genome/genomic_codon_positions.csv.gz', 
                 callback = DataFrameCallback$new(preprocess_codons),
                 chunk_size = 100000, col_types = 'cccciiciic')

#### Set up chromosome sequence ####
chromosome_groups <- tibble(chromosome = names(genome),
       length = width(genome)) %>%
  mutate(chrom_group = row_number() %% 2 == 1) %>%
  mutate(concat_start_position = cumsum(lag(length, default = 0L)))

#### Join SNP info ####
full_snp_data <- left_join(lfmm_results, snp_data,
                 by = c('chromosome', 'position', 'locus')) %>%
  find_snp_gene_overlap(genes = gene_locations) %>%
  left_join(codon_sequences,
            by = c('chromosome' = 'Seqid',
                   'position',
                   'locus',
                   'gene_id',
                   'gene_strand' = 'Strand')) %>%
  
  mutate(across(c(ref, alt, dna_base, codon_seq), dna)) %>%
  dplyr::rename(ref_codon_seq = codon_seq) %>%
  mutate(alt_codon_seq = seq_replace_position(ref_codon_seq, position_in = codon_position, 
                                              position_out = codon_position, 
                                              replacement = alt),
         ref_aa = seq_translate(ref_codon_seq),
         alt_aa = seq_translate(alt_codon_seq),
         is_synonomous = ref_aa == alt_aa,
         is_stop = alt_aa == aa("*")) %>%
  mutate(mutation_type = case_when(is_stop ~ 'Nonsense',
                                   is_synonomous ~ 'Synonomous',
                                   !is_synonomous ~ 'Non-Synonomous',
                                   !is.na(gene_id) ~ 'Non-Coding in Gene',
                                   is.na(gene_id) ~ 'Non-Coding outside Gene',
                                   TRUE ~ 'unidentified_mutation'))
write_csv(full_snp_data, '../intermediate_files/snps_with_mutationType.csv')

#### SNP Table ####
overall_lfmm_summary <- full_snp_data %>% 
  select(chromosome, es, fdr, gene_id, mutation_type) %>%
  group_by(mutation_type, chromosome, gene_id) %>%
  
  summarise(n_snps = n(),
            sig_snps = sum(fdr < alpha),
            
            n_positive = sum(es > 0),
            sig_positive = sum(es > 0 & fdr < alpha),
            
            n_negative = sum(es < 0),
            sig_negative = sum(es < 0 & fdr < alpha),
            
            .groups = 'drop_last') %>%
  # group_by(codon_position, mutation_type, chromosome) %>%
  summarise(across(where(is.integer), sum), .groups = 'drop_last') %>%
  summarise(across(where(is.integer), sum), .groups = 'drop')
write_csv(overall_lfmm_summary, '../Results/overall_lfmm_summary.csv')  

gene_lfmm_summary <- full_snp_data %>% 
  select(chromosome, es, fdr, gene_id, codon_position, mutation_type) %>%
  # filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>%
  filter(any(fdr < alpha)) %>%
  ungroup %>%
  group_by(chromosome, gene_id, mutation_type, codon_position) %>%
  
  summarise(n_snps = n(),
            sig_snps = sum(fdr < alpha),
            
            n_positive = sum(es > 0),
            sig_positive = sum(es > 0 & fdr < alpha),
            
            n_negative = sum(es < 0),
            sig_negative = sum(es < 0 & fdr < alpha),
            
            .groups = 'drop_last') %>%
  summarise(across(where(is.integer), sum), .groups = 'drop') %>%
  left_join(select(functional_annotations, gene_id, structure_attributes),
            by = 'gene_id') 
write_csv(gene_lfmm_summary, '../Results/gene_lfmm_summary.csv')  

gene_info <- full_snp_data %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>%
  filter(any(fdr < alpha)) %>%
  select(gene_id) %>%
  distinct %>%
  left_join(functional_annotations,
            by = 'gene_id')
write_csv(gene_info, '../Results/gene_info_summary.csv')

#### P-Value/ES over Genome ####
#arrange so significant points are plotted on top and so the mutation types plot on top in order or interest
preprocess_genome_snps <- full_snp_data %>%
  left_join(chromosome_groups, 
            by = 'chromosome') %>%
  mutate(coloration = case_when(fdr < alpha ~ 'significant',
                                chrom_group ~ 'A',
                                !chrom_group ~ 'B')) %>%
  mutate(mutation_type = if_else(str_detect(mutation_type, 'Non-Coding'), 'Non-Coding', mutation_type) %>%
           factor(levels = c('Nonsense', 'Non-Synonomous', 'Synonomous', 'Non-Coding'))) %>%
  mutate(position_concat = position + concat_start_position) %>%
  arrange(coloration, desc(mutation_type)) 

L <- nrow(lfmm_results)
aggregate_p <- function(p){
  w <- rep(1/L, length = length(p))
  p_hmp <- p.hmp(p, w, L)
  adj_p <- p_hmp / sum(w)
  tibble(p_hmp = p_hmp, adj_p = adj_p)
}

make_p_es_plots <- function(data, universal_p = FALSE, add_annotations = FALSE){
  p_plot <- data %>%
    # sample_frac(0.01) %>%
    ggplot(aes(x = position_concat / 1e6, y = -log10(p), size = fdr < alpha, 
               colour = coloration, shape = mutation_type)) +
    geom_point() +
    scale_size_manual(values = c('TRUE' = 1.5, 'FALSE' = 0.2)) +
    scale_colour_manual(values = c('significant' = 'firebrick2', 'A' = 'grey60', 
                                   'B' = 'grey40')) +
    scale_shape_manual(values = c('Nonsense' = 'square',
                                  'Non-Synonomous' = 'diamond',
                                  'Synonomous' = 'triangle',
                                  'Non-Coding' = 'circle')) +
    scale_x_continuous(labels = scales::comma_format()) +
    guides(size = 'none',
           colour = 'none',
           shape = guide_legend(override.aes = list(size = 2))) +
    labs(x = 'Position (MB)',
         y = '-log10(p value)',
         shape = 'Mutation Type') +
    theme_classic() +
    theme(axis.text = element_text(face = 'plain', colour = 'black', size = 24),
          axis.title = element_text(face = 'plain', colour = 'black', size = 24),
          panel.border = element_rect(colour = 'black', fill = NA),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  es_plot <- data %>%
    # sample_frac(0.01) %>%
    ggplot(aes(x = position_concat / 1e6, y = es, size = fdr < alpha, 
               colour = coloration, shape = mutation_type)) +
    geom_point() +
    scale_size_manual(values = c('TRUE' = 1.5, 'FALSE' = 0.2)) +
    scale_colour_manual(values = c('significant' = 'firebrick2', 'A' = 'grey60', 
                                   'B' = 'grey40')) +
    scale_shape_manual(values = c('Nonsense' = 'square',
                                  'Non-Synonomous' = 'diamond',
                                  'Synonomous' = 'triangle',
                                  'Non-Coding' = 'circle')) +
    scale_x_continuous(labels = scales::comma_format()) +
    guides(size = 'none',
           colour = 'none',
           shape = guide_legend(override.aes = list(size = 2))) +
    labs(x = 'Position (MB)',
         y = 'Effect Size',
         shape = 'Mutation Type') +
    theme_classic() +
    theme(axis.text = element_text(face = 'plain', colour = 'black', size = 24),
          axis.title = element_text(face = 'plain', colour = 'black', size = 24),
          panel.border = element_rect(colour = 'black', fill = NA))
  
  
  if(universal_p){
    L <- nrow(lfmm_results); library(harmonicmeanp)
    significant_chuncks <- preprocess_genome_snps %>%
      group_by(chromosome) %>%
      summarise(aggregate_p(p)) %>%
      filter(adj_p < alpha) %>%
      left_join(preprocess_genome_snps, by = 'chromosome') %>%
      group_by(chromosome) %>%
      filter(position == min(position) |
               position == max(position)) %>%
      select(chromosome, position_concat) %>%
      mutate(which_spot = c('start', 'stop')) %>%
      pivot_wider(names_from = 'which_spot',
                  values_from = 'position_concat') %>%
      ungroup
    
    
    p_plot <- p_plot +
      geom_segment(data = significant_chuncks, 
                   aes(x = start/1e6, xend = stop/1e6,
                       y = 1.1 * max(-1 * log(data$p, base = 10)), 
                       yend = 1.1 * max(-1 * log(data$p, base = 10))), 
                   inherit.aes = FALSE)
    
    es_plot <- es_plot +
      geom_segment(data = significant_chuncks, 
                   aes(x = start/1e6, xend = stop/1e6,
                       y = 1.1 * max(data$es), 
                       yend = 1.1 * max(data$es)), 
                   inherit.aes = FALSE)
  }
  
  if(add_annotations){
    
    annotation_data <- data %>%
      filter(fdr < alpha) %>%
      group_by(chromosome) %>% #count(mutation_type)
      summarise(position_concat = median(position_concat),
                p = max(-1 * log(p, base = 10)),
                n_snp = n(),
                n_in_gene = sum(!is.na(gene_id)),
                n_functional = sum(mutation_type == 'Non-Synonomous'),
                .groups = 'drop') %>%
      mutate(annotation = str_c(n_snp, n_in_gene, n_functional, sep = ', '),
             annotation = str_c(' (', annotation, ')', sep = ''))
      
    library(ggrepel) #1.1 * max(-1 * log(data$p, base = 10))
    p_plot <- p_plot +
      geom_text(data = annotation_data, 
                aes(x = position_concat / 1e6, 
                    y = p,
                    label = annotation),
                inherit.aes = FALSE, hjust = 0,
                vjust = 0)
  }
  
  list(p_plot, es_plot)
}

chromosome_plot <- make_p_es_plots(preprocess_genome_snps, 
                                   universal_p = TRUE,
                                   add_annotations = TRUE) %>%
  wrap_plots() +
  plot_layout(guides = 'collect', ncol = 1) +
  plot_annotation(tag_levels = 'A') &
  theme(legend.title = element_text(colour = 'black', size = 24),
        legend.text = element_text(colour = 'black', size = 20),
        plot.tag = element_text(size = 28))
ggsave('../Results/all_contig_p_es.png', plot = chromosome_plot, height = 10, width = 15)
