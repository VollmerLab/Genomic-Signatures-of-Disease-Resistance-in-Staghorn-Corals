#### Linkage Block plotting ####

## Goal - function inputs = genomic range. output = dataframe with all mutation types, what genes are included, mutation type broken out by gene
## Goal - function inputs = genomic range & buffer = plot with gene diagram, p-value, es of all snps and panels breaking out each linkage region

#### Libraries ####
library(IRanges)
library(bioseq)
library(tidyverse)
library(magrittr)
library(patchwork)

rename <- dplyr::rename

alpha <- 0.05
linkage_buffer <- 1000
window_buffer  <- 1000
pgs_p_cutoff <- 0.0001
position_divisor <- 1e6
only_fdr_sig_clumps <- TRUE

#### Functions ####
preprocess_codons <- function(x, pos){
  mutate(x, locus = str_c(Seqid, position, sep = '-')) %>%
    filter(locus %in% lfmm_results$locus)
}

collapse_clumps <- function(lower, upper, buffer){
  linkage_ranges <- IRanges(start = lower - (buffer / 2), end = upper + (buffer / 2))
  
  collapsed_groups <- findOverlaps(linkage_ranges, IRanges::reduce(linkage_ranges)) %>% subjectHits()
  tibble(start = lower, stop = upper, plot_group = str_c('P', collapsed_groups)) %>%
    group_by(plot_group) %>%
    summarise(min = min(start), max = max(stop),
              .groups = 'drop')
}

create_gene_diagram <- function(gene_info, region_start, region_stop, start, stop){
  range_of_interest <- c(start, stop)
  
  gene_info %>%
    mutate(min_value = if_else(Strand == '+', 0, -0.25),
           max_value = if_else(Strand == '+', 0.25, 0),
           gene_mid = (gene_end + gene_start) / 2,
           gene_id = str_remove(gene_id, 'Acer_0+'),
           ymid = c(min_value + max_value) / 2) %>%
    
    ggplot(aes(xmin = gene_start / position_divisor, 
               xmax = gene_end / position_divisor, 
               ymin = min_value, ymax = max_value * (Strand == '+'))) +
    
    geom_rect(fill = 'grey50') +
    geom_text(aes(x = gene_mid / position_divisor, y = ymid, label = gene_id),
              hjust = 0.5, vjust = 0.5, colour = 'white') +
    
    geom_segment(data = tibble(start = start, stop = stop, 
                               y = if_else(any(gene_info$Strand == '-'), -0.26, -0.01)),
                 inherit.aes = FALSE, 
                 aes(x = start / position_divisor, xend = stop / position_divisor,
                     y = 0, yend = 0), linewidth = 2, colour = 'black', linetype = 'dashed') +
    
    geom_segment(data = tibble(start = region_start, stop = region_stop, 
                               y = if_else(any(gene_info$Strand == '-'), -0.26, -0.01)),
                 inherit.aes = FALSE, colour = 'black',
                 aes(x = start / position_divisor, xend = stop / position_divisor,
                     y = 0, yend = 0), linewidth = 2) +
    
    scale_x_continuous(limits = (range_of_interest + c(-1, 1)) / position_divisor) +
    scale_y_continuous(limits = c(-0.25, 0.25)) +
    theme_void()   +
    theme(axis.text = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text())
  
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

aggregate_p <- function(p){
  library(harmonicmeanp)
  L <- nrow(lfmm_results)
  w <- rep(1/L, length = length(p))
  p_hmp <- p.hmp(p, w, L)
  adj_p <- p_hmp / sum(w)
  tibble(p_hmp = p_hmp, adj_p = adj_p)
}

tabulate_functional_mutations <- function(snp_info, clumps){
  
  unified_p <- snp_info %>%
    group_by(gene_id) %>%
    summarise(aggregate_p(p))
  
  all <- snp_info %>%
    group_by(gene_id, clump_id) %>%
    summarise(S = sum(mutation_type == 'Synonomous'),
              N = sum(mutation_type == 'Non-Synonomous'),
              .groups = 'drop_last') %>%
    summarise(n_clump = n(),
              clump_names = str_c(clump_id, collapse = '; '),
              across(where(is.integer), sum),
              .groups = 'drop') 
  
  fdr_sig <- snp_info %>%
    group_by(gene_id, clump_id) %>%
    filter(any(fdr < alpha)) %>%
    summarise(S = sum(mutation_type == 'Synonomous'),
              N = sum(mutation_type == 'Non-Synonomous'),
              .groups = 'drop_last') %>%
    summarise(n_clump = n(),
              clump_names = str_c(clump_id, collapse = '; '),
              across(where(is.integer), sum),
              .groups = 'drop') 
  
  in_pgs_clump <- snp_info %>%
    filter(clump_id %in% clumps) %>%
    group_by(gene_id, clump_id) %>%
    summarise(S = sum(mutation_type == 'Synonomous'),
              N = sum(mutation_type == 'Non-Synonomous'),
              .groups = 'drop_last') %>%
    summarise(n_clump = n(),
              clump_names = str_c(clump_id, collapse = '; '),
              across(where(is.integer), sum),
              .groups = 'drop')
  
  bind_rows(all = all, fdr = fdr_sig, pgs_inclusion = in_pgs_clump, .id = 'cutoff') %>%
    filter(!is.na(gene_id)) %>%
    complete(gene_id, cutoff = c('all', 'fdr', 'pgs_inclusion'), fill = list(n_clump = 0L, S = 0L, N = 0L)) %>%
    mutate(dNdS = N / S) %>%
    pivot_wider(names_from = 'cutoff',
                values_from = c('n_clump', 'clump_names', 'S', 'N', 'dNdS'),
                names_vary = 'slowest') %>%
    left_join(unified_p, by = 'gene_id') %>%
    left_join(functional_annotations, by = 'gene_id') %>%
    select(-starts_with('has_'), -blast_details, -interproscan_results,
           -p_hmp, -adj_p) %>%
    select(gene_id, everything())
}

make_es_p_plot <- function(plot_option, locus_info, linkage_start, linkage_stop, 
                           range_start, range_stop, plot_region = TRUE, 
                           clump_inclusion = NULL){
  #plot option must either be p or ES - to plot subset of loci filter the locus_info to specified subset
  
  range_of_interest <- c(range_start, range_stop)
  
  region_of_interest <- tibble(region_start = linkage_start,
                               region_end = linkage_stop)
  
  if(all(!is.null(clump_inclusion))){
    locus_info <- filter(locus_info, clump_id %in% clump_inclusion)
  }
  
  if(plot_option == 'p'){
    base_plot <- locus_info %>%
      arrange(-1 * log(p, base = 10)) %>%
      
      ggplot(aes(x = position / position_divisor, y = -1 * log(p, base = 10)))
  } else if(plot_option %in% c('es', 'functional')) {
    
    base_plot <- locus_info %>%
      arrange(-1 * log(p, base = 10)) %>%
      
      ggplot(aes(x = position / position_divisor, y = es))
    
  } else {
    return(NULL)
  }
  
  
  if(plot_region){
    if(linkage_start != linkage_stop){
      interest_added <- base_plot +
        geom_rect(data = region_of_interest, inherit.aes = FALSE,
                  aes(xmin = region_start / position_divisor, xmax = region_end / position_divisor,
                      ymin = -Inf, ymax = Inf), fill = 'grey50')
    } else {
      interest_added <- base_plot +
        geom_vline(xintercept = linkage_start / position_divisor, colour = 'grey50')
    }
  } else {
    interest_added <- base_plot
  }
  
  if(plot_option %in% c('es', 'functional')){
    interest_added <- interest_added +
      geom_hline(yintercept = 0, linetype = 'dashed')
  }
  
  
  if(plot_option == 'functional'){
    points_added <- interest_added +
      geom_point(aes(colour = -log(p, base = 10),
                     size = coding,
                     shape = mutation_type)) +
      binned_scale(aesthetics = "color",
                   scale_name = "stepsn",
                   palette = function(x) c("grey50",  RColorBrewer::brewer.pal(4, 'Reds')[2:4]),
                   breaks = c(2, 4, 6),
                   limits = c(0, 8),
                   show.limits = FALSE,
                   guide = "colorsteps",
                   oob = scales::squish)
  } else {
    points_added <- interest_added +
      geom_point(aes(colour = fdr < alpha,
                     size = coding,
                     shape = mutation_type)) +
      scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'black'), drop = FALSE) +
      guides(colour = 'none')
  }
  
  points_added +
    scale_x_continuous(limits = range_of_interest / position_divisor) +
    scale_size_manual(values = c('TRUE' = 2, 'FALSE' = 0.5), drop = FALSE) +
    scale_shape_manual(values = c('Non-Coding' = 'circle',
                                  'Synonomous' = 'square',
                                  'Non-Synonomous' = 'triangle',
                                  'Nonsense' = 'cross'), 
                       drop = FALSE) +
    guides(size = 'none',
           shape = guide_legend(override.aes = list(size = 3))) +
    labs(x = 'Position (MB)',
         y = if_else(plot_option == 'p', '-log10(p)', 'Effect Size'),
         # colour = '-log10(p)',
         colour = bquote(-log[10](p)),
         shape = 'Mutation Type') +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_text(colour = 'black', size = 16),
          panel.border = element_rect(colour = 'black', fill = NA),
          axis.text = element_text(colour = 'black', size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.text = element_text(colour = 'black', size = 16),
          legend.title = element_text(colour = 'black', size = 14))
}


#### Data ####
lfmm_results <- read_csv('../Results/lfmm_results.csv', show_col_types = FALSE) %>%
  filter(pop_group == 'all',
         model_family == 'norm') %>%
  select(-pop_group, -model_family) %>%
  mutate(locus = str_c(chromosome, position, sep = '-'), .after = 'position') 

gene_locations <- read_csv('../../Bioinformatics/genome_annotation/r5_gene_id_locations.csv', show_col_types = FALSE) %>%
  filter(gene_type != 'tRNA')
structural_annotations <- read_csv('../../Bioinformatics/genome_annotation/r5_structural_annotation.csv', show_col_types = FALSE)
functional_annotations <- read_rds('../../Bioinformatics/genome_annotation/r5_functional_annotations.rds.gz') %>%
  dplyr::rename(gene_id = qseqid)
codon_sequences <- read_csv_chunked('../../Bioinformatics/genome_annotation/r5_codon_positions.csv.gz', 
                                    callback = DataFrameCallback$new(preprocess_codons),
                                    chunk_size = 100000, col_types = 'cccciiciic')

linkage_blocks <- read_csv('../intermediate_files/lfmm_snp_clumps.csv', show_col_types = FALSE)
threshold_linkage_inclusion <- read_csv('../intermediate_files/clumping_thersholding_results.csv', show_col_types = FALSE) %>%
  filter(threshold_p == pgs_p_cutoff) %>%
  select(chromosome, clump_id)

snp_data <- read_delim('../../Bioinformatics/variant_calling/18October2022/genolike.mafs.gz', delim = '\t', col_types = 'cicccnnnni') %>%
  select(chromo, position, major, minor, ref) %>%
  dplyr::rename(chromosome = chromo) %>%
  mutate(locus = str_c(chromosome, position, sep = '-'),
         alt = if_else(ref == major, minor, major))

#### Summarise SNPs in genes retained by PGS calculation ####
all_snp_info <- linkage_blocks %>%
  left_join(lfmm_results, 
            by = c('chromosome', 'locus')) %>%
  left_join(select(snp_data, -major, -minor),
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
                                              position_out = codon_position, replacement = alt),
         ref_aa = seq_translate(ref_codon_seq),
         alt_aa = seq_translate(alt_codon_seq),
         is_synonomous = ref_aa == alt_aa,
         is_stop = alt_aa == aa("*")) %>%
  mutate(mutation_type = case_when(is_stop ~ 'Nonsense',
                                   is_synonomous ~ 'Synonomous',
                                   !is_synonomous ~ 'Non-Synonomous',
                                   !is.na(gene_id) ~ 'Non-Coding in Gene',
                                   is.na(gene_id) ~ 'Non-Coding outside Gene',
                                   TRUE ~ 'unidentified_mutation')) %>%
  mutate(mutation_type = if_else(str_detect(mutation_type, 'Non-Coding'), 'Non-Coding', mutation_type) %>%
           factor(levels = c('Nonsense', 'Non-Synonomous', 'Synonomous', 'Non-Coding')) %>% fct_rev,
         coding = mutation_type %in% c('Non-Synonomous', 'Synonomous')) %>%
  
  select(chromosome, clump_id, locus, position, is_index, es, z, p, fdr, mutation_type, gene_id) %>%
  # mutate(fdr = fdr < alpha) %>%
  left_join(select(functional_annotations, gene_id, starts_with('swissprot')), 
            by = 'gene_id') %>%
  mutate(chromosome = str_remove_all(chromosome, '[A-Za-z_]+') %>% as.integer(),
         clump_id = str_remove_all(clump_id, '[A-Za-z_]+'),
         locus_id = str_c(chromosome, position, clump_id, sep = '_')) %>%
  distinct %>%
  select(locus_id, everything()) %>%
  arrange(chromosome, position) %>%
  
  group_by(chromosome, clump_id) %>%
  mutate(min_linkage_fdr = min(fdr)) %>%
  ungroup %>%
  select(locus_id:is_index, min_linkage_fdr, everything()) %>%
  rename(linkage_group = clump_id)
write_csv(all_snp_info, '../intermediate_files/all_snp_info.csv')

all_snp_info %>%
  filter(fdr < alpha) %>%
  count(chromosome, linkage_group) %>%
  group_by(chromosome) %>%
  summarise(n_lg = n(),
            n_snp = sum(n))

gene_summary_overall <- all_snp_info %>%
  left_join(mutate(threshold_linkage_inclusion, 
                   chromosome = str_extract(chromosome, '[0-9]+') %>% as.integer,
                   clump_id = str_extract(clump_id, '[0-9]+'),
                   in_plink = TRUE),
            by = c('chromosome', 'linkage_group' = 'clump_id')) %>%
  mutate(in_plink = (in_plink & min_linkage_fdr < alpha)) %>%
  group_by(chromosome, gene_id, 
           across(starts_with('swissprot')),
           mutation_type) %>%
  summarise(total_snps = n_distinct(locus),
            plink_fdr = n_distinct(locus[in_plink]),
            lea_snps = n_distinct(locus[fdr < alpha]),
            .groups = 'drop') %>%
  ungroup %>%
  pivot_wider(names_from = mutation_type,
              values_from = c('total_snps', 'plink_fdr', 'lea_snps'),
              values_fill = 0L) %>%
  
  pivot_longer(cols = c(starts_with('total_'), starts_with('plink_'), starts_with('lea_')),
               names_to = c('grouping', '.value'),
               names_pattern = '([a-z_]+)_(.*)') %>%
  rowwise %>%
  mutate(Total = sum(c_across(c(where(is.integer), -chromosome)))) %>%
  ungroup %>%
  pivot_wider(names_from = grouping, 
              values_from = c(where(is.integer), -chromosome),
              names_vary = 'slowest',
              names_glue = '{grouping}_{.value}') 

gene_summary_overall %>%
  filter(plink_fdr_Total > 0) %>%
  write_csv('../Results/summarized_clumpGene_mutations.csv', na = '')

#### Summarize SNPs in PGS calculation ####
clump_snp_summary <- threshold_linkage_inclusion %>%
  left_join(linkage_blocks, by = c('chromosome', 'clump_id')) %>%
  left_join(lfmm_results, 
            by = c('chromosome', 'locus')) %>%
  left_join(select(snp_data, -major, -minor),
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
                                              position_out = codon_position, replacement = alt),
         ref_aa = seq_translate(ref_codon_seq),
         alt_aa = seq_translate(alt_codon_seq),
         is_synonomous = ref_aa == alt_aa,
         is_stop = alt_aa == aa("*")) %>%
  mutate(mutation_type = case_when(is_stop ~ 'Nonsense',
                                   is_synonomous ~ 'Synonomous',
                                   !is_synonomous ~ 'Non-Synonomous',
                                   !is.na(gene_id) ~ 'Non-Coding in Gene',
                                   is.na(gene_id) ~ 'Non-Coding outside Gene',
                                   TRUE ~ 'unidentified_mutation')) %>%
  mutate(mutation_type = if_else(str_detect(mutation_type, 'Non-Coding'), 'Non-Coding', mutation_type) %>%
           factor(levels = c('Nonsense', 'Non-Synonomous', 'Synonomous', 'Non-Coding')) %>% fct_rev,
         coding = mutation_type %in% c('Non-Synonomous', 'Synonomous')) %>%
  
  select(chromosome, clump_id, locus, position, is_index, es, z, p, fdr, mutation_type, gene_id) %>%
  # mutate(fdr = fdr < alpha) %>%
  left_join(select(functional_annotations, gene_id, starts_with('swissprot')), 
            by = 'gene_id') %>%
  mutate(chromosome = str_remove_all(chromosome, '[A-Za-z_]+') %>% as.integer(),
         clump_id = str_remove_all(clump_id, '[A-Za-z_]+'),
         locus_id = str_c(chromosome, position, clump_id, sep = '_')) %>%
  distinct %>%
  select(locus_id, everything()) %>%
  arrange(chromosome, position) %>%
  
  group_by(chromosome, clump_id) %>%
  mutate(min_linkage_fdr = min(fdr)) %>%
  ungroup %>%
  select(locus_id:is_index, min_linkage_fdr, everything()) %>%
  rename(linkage_group = clump_id)

write_csv(clump_snp_summary, '../Results/clump_snp_summary.csv')
# filter(clump_snp_summary, min_linkage_fdr < alpha) %>% View

#### Summarize Clump SNP summary ####

#filter lfmm fdr significant
#chromosome, lg, gene_id, count mutation types

#for each ID number of each mutation type that are FDR significant
#within each PGS contig, lg, gene_id, count mutation types, (all snps)

clump_summary_squared <- clump_snp_summary %>% 
  group_by(chromosome, linkage_group, gene_id, 
           across(starts_with('swissprot')), 
           min_linkage_fdr, mutation_type) %>%
  summarise(total_snps = n_distinct(locus),
            plink_fdr = n_distinct(locus[min_linkage_fdr < alpha]),
            lea_snps = n_distinct(locus[fdr < alpha]),
            .groups = 'drop') %>%
  ungroup %>%
  pivot_wider(names_from = mutation_type,
              values_from = c('total_snps', 'plink_fdr', 'lea_snps'),
              values_fill = 0L) %>%
  
  pivot_longer(cols = c(starts_with('total_'), starts_with('plink_'), starts_with('lea_')),
               names_to = c('grouping', '.value'),
               names_pattern = '([a-z_]+)_(.*)') %>%
  rowwise %>%
  mutate(Total = sum(c_across(c(where(is.integer), -chromosome)))) %>%
  ungroup %>%
  pivot_wider(names_from = grouping, 
              values_from = c(where(is.integer), -chromosome),
              names_vary = 'slowest',
              names_glue = '{grouping}_{.value}') 

clump_summary_squared <- clump_summary_squared %>%
  group_by(chromosome, gene_id, across(starts_with('swissprot'))) %>%
  summarise(across(where(is.integer), sum), .groups = 'drop') 
write_csv(clump_summary_squared, '../Results/summarized_clumpGene_mutations.csv', na = '')


#### Go from Linkage Blocks to Regions of interest ####
# Step 1 - get range of each important linkage block
# Step 2 - combine overlapping linkage blocks

non_overlapping_clumps <- threshold_linkage_inclusion %>%
  left_join(linkage_blocks, by = c('chromosome', 'clump_id')) %>%
  left_join(select(lfmm_results, chromosome, position, locus, fdr), 
            by = c('chromosome', 'locus')) %>%
  
  {if(only_fdr_sig_clumps) group_by(., clump_id) %>%
      filter(any(fdr < alpha)) %>%
      ungroup else .} %>%
  
  select(-is_index, -locus, -fdr) %>%
  group_by(chromosome, clump_id) %>%
  filter(position == min(position) | position == max(position)) %>%
  mutate(location = if_else(position == min(position), 'min', 'max')) %>%
  pivot_wider(names_from  = 'location',
              values_from = 'position') %>%
  mutate(max = if_else(is.na(max), min, max)) %>%
  
  group_by(chromosome) %>%
  # filter(n() > 1) %>%
  summarise(collapse_clumps(min, max, linkage_buffer), 
            clumps = list(clump_id),
            .groups = 'drop') %>%
  rename(region_start = min,
         region_end = max)

non_overlapping_clumps %>%
  arrange(str_extract(chromosome, '[0-9]+') %>% as.integer) %>%
  mutate(non_overlapping_clump = str_c('nonOverlap', row_number(), sep = '_')) %>%
  select(chromosome, clumps, non_overlapping_clump) %>%
  unnest(clumps)

#### Summarize Index SNPs ####
index_snp_summary <- threshold_linkage_inclusion %>%
  left_join(linkage_blocks, by = c('chromosome', 'clump_id')) %>%
  left_join(lfmm_results, 
            by = c('chromosome', 'locus')) %>%
  filter(is_index) %>%
  left_join(select(snp_data, -major, -minor),
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
                                              position_out = codon_position, replacement = alt),
         ref_aa = seq_translate(ref_codon_seq),
         alt_aa = seq_translate(alt_codon_seq),
         is_synonomous = ref_aa == alt_aa,
         is_stop = alt_aa == aa("*")) %>%
  mutate(mutation_type = case_when(is_stop ~ 'Nonsense',
                                   is_synonomous ~ 'Synonomous',
                                   !is_synonomous ~ 'Non-Synonomous',
                                   !is.na(gene_id) ~ 'Non-Coding in Gene',
                                   is.na(gene_id) ~ 'Non-Coding outside Gene',
                                   TRUE ~ 'unidentified_mutation')) %>%
  mutate(mutation_type = if_else(str_detect(mutation_type, 'Non-Coding'), 'Non-Coding', mutation_type) %>%
           factor(levels = c('Nonsense', 'Non-Synonomous', 'Synonomous', 'Non-Coding')) %>% fct_rev,
         coding = mutation_type %in% c('Non-Synonomous', 'Synonomous')) %>%
  
  select(chromosome, clump_id, locus, position, es, z, p, fdr, mutation_type, gene_id) %>%
  # mutate(fdr = fdr < alpha) %>%
  left_join(select(functional_annotations, gene_id, starts_with('swissprot')), 
            by = 'gene_id') %>%
  mutate(chromosome = str_remove_all(chromosome, '[A-Za-z_]+') %>% as.integer(),
         clump_id = str_remove_all(clump_id, '[A-Za-z_]+'),
         locus_id = str_c(chromosome, position, clump_id, sep = '_')) %>%
  distinct %>%
  select(locus_id, everything()) %>%
  arrange(chromosome, position)

write_csv(index_snp_summary, '../Results/index_snp_summary.csv')

#### Join Non-overlapping regions with genes in those regions ####
## 

gene_clumps <- non_overlapping_clumps %>%
  left_join(nest(gene_locations, gene_data = -c(Seqid)),
            by = c('chromosome' = 'Seqid')) %>%
  unnest(gene_data) %>%
  filter(gene_end >= region_start - (window_buffer / 2),
         gene_start <= region_end + (window_buffer / 2)) %>%
  filter(gene_type != 'tRNA') %>%
  nest(gene_data = c(gene_id, gene_start, gene_end, Strand, gene_type)) %>%
  rowwise %>%
  mutate(plot_start = min(c(region_start, gene_data$gene_start)),
         plot_end = max(c(region_end, gene_data$gene_end))) %>%
  # mutate(plot_start = min(c(gene_data$gene_start)),
  #        plot_end = max(c(gene_data$gene_end))) %>%
  ungroup 

gene_clumps$gene_data[[2]] %>%
  left_join(functional_annotations, by = 'gene_id') %>%
  select(gene_id, gene_type, swissprot_name)



#### Get SNP Info for each Joined Region ####
snp_info_joined <- gene_clumps %>%
  left_join(lfmm_results,
            by = c('chromosome')) %>%
  filter(position > plot_start,
         position < plot_end) %>%
  left_join(linkage_blocks,
            by = c('chromosome', 'locus')) %>%
  left_join(select(snp_data, -major, -minor),
            by = c('chromosome', 'position', 'locus')) %>%
  
  find_snp_gene_overlap(genes = gene_locations) %>% # increase # rows
  left_join(codon_sequences,
            by = c('chromosome' = 'Seqid',
                   'position',
                   'locus',
                   'gene_id',
                   'gene_strand' = 'Strand')) %>%
  
  mutate(across(c(ref, alt, dna_base, codon_seq), dna)) %>%
  dplyr::rename(ref_codon_seq = codon_seq) %>%
  mutate(alt_codon_seq = seq_replace_position(ref_codon_seq, position_in = codon_position, 
                                              position_out = codon_position, replacement = alt),
         ref_aa = seq_translate(ref_codon_seq),
         alt_aa = seq_translate(alt_codon_seq),
         is_synonomous = ref_aa == alt_aa,
         is_stop = alt_aa == aa("*")) %>%
  mutate(mutation_type = case_when(is_stop ~ 'Nonsense',
                                   is_synonomous ~ 'Synonomous',
                                   !is_synonomous ~ 'Non-Synonomous',
                                   !is.na(gene_id) ~ 'Non-Coding in Gene',
                                   is.na(gene_id) ~ 'Non-Coding outside Gene',
                                   TRUE ~ 'unidentified_mutation')) %>%
  mutate(mutation_type = if_else(str_detect(mutation_type, 'Non-Coding'), 'Non-Coding', mutation_type) %>%
           factor(levels = c('Nonsense', 'Non-Synonomous', 'Synonomous', 'Non-Coding')) %>% fct_rev,
         coding = mutation_type %in% c('Non-Synonomous', 'Synonomous')) %>%
  nest(locus_data = c(locus, position, es, z, p, fdr, 
                      clump_id, is_index, ref, alt, 
                      gene_id, gene_start, gene_end, gene_type, gene_strand,
                      cds_id, structure_phase, dna_base, codon, codon_position,
                      ref_codon_seq, alt_codon_seq, ref_aa, alt_aa, 
                      is_synonomous, is_stop, mutation_type, coding)) 

gene_info_summary <- snp_info_joined %>%
  rowwise(chromosome, plot_group) %>%
  summarise(tabulate_functional_mutations(locus_data, clumps),
            .groups = 'drop') %>%
  unnest(entap_results, keep_empty = TRUE) %>%
  select(chromosome, plot_group, gene_id, ends_with('all'), ends_with('fdr'), ends_with('pgs_inclusion'), everything()) %>%
  arrange(gene_id)


gene_info_summary %>%
  filter(gene_id == 'Acer_00005579') %>%
  select(contains('KEGG'))


gene_info_summary %>%
  filter(n_clump_fdr > 0) %>%
  select(chromosome, gene_id)

write_csv(gene_info_summary, '../Results/linkage_group_genes.csv')


#### Make Overall with functional changes ####
chromosome_plots <- snp_info_joined %>%
  
  # filter(chromosome == 'Acerv_scaffold_115') %>%
  
  expand_grid(tibble(plot_type = c('gene_diagram', 'p', 'es', 'nonSyn_all', 'nonSyn_contributePGS'))) %>%
  distinct %>%
  rowwise %>%
  mutate(locus_data = case_when(str_detect(plot_type, 'functional') ~ list(filter(locus_data, 
                                                                                  mutation_type != 'Non-Coding')),
                                str_detect(plot_type, 'nonSyn') ~ list(filter(locus_data, 
                                                                                  mutation_type == 'Non-Synonomous')),
                                TRUE ~ list(locus_data)),
         plot = case_when(plot_type == 'gene_diagram' ~ list(create_gene_diagram(gene_data, 
                                                                                 region_start, region_end, 
                                                                                 plot_start, plot_end)),
                          plot_type == 'p' ~ list(make_es_p_plot('p', locus_data, 
                                                                 region_start, region_end,
                                                                 plot_start, plot_end,
                                                                 plot_region = FALSE)),
                          plot_type == 'nonSyn_all' ~ list(make_es_p_plot('functional', locus_data, 
                                                                          region_start, region_end, 
                                                                          plot_start, plot_end, 
                                                                          plot_region = FALSE,
                                                                      clump_inclusion = NULL)),
                          plot_type == 'nonSyn_contributePGS' ~ list(make_es_p_plot('functional', locus_data, 
                                                                      region_start, region_end, 
                                                                      plot_start, plot_end, 
                                                                      plot_region = FALSE,
                                                                      clump_inclusion = clumps)),
                          # plot_type == 'clump_locations' ~ list(make_clump_plot(locus_data, 
                          #                                                           region_start, region_end, 
                          #                                                           plot_start, plot_end,
                          #                                                           clump_inclusion = clumps)),
                          TRUE ~ list(make_es_p_plot('es', locus_data, 
                                                     region_start, region_end, 
                                                     plot_start, plot_end,
                                                     plot_region = FALSE)))) %>%
  filter(nrow(locus_data) > 0) %>%
  ungroup %>%
  
  group_by(chromosome, plot_group) %>%
  mutate(position = row_number(),
         max_pos = position != max(position)) %>%
  rowwise %>%
  mutate(plot = if_else(max_pos,
                        list(plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                                          axis.ticks.x = element_blank())),
                        list(plot))) %>%
  
  group_by(chromosome, plot_group) %>%
  summarise(linkage_plot = list(wrap_plots(plot, guides = 'collect', ncol = 1, 
                                           height = c(0.05, rep(0.95 / (n() - 1), (n() - 1))))),
            .groups = 'drop_last') %>%
  
  mutate(keep_legend = row_number() == max(row_number())) %>%
  rowwise %>%
  mutate(linkage_plot = if_else(keep_legend, 
                                list(linkage_plot), 
                                list(linkage_plot & theme(legend.position = 'none')))) %>%
  
  group_by(chromosome, plot_group) %>%
  summarise(chromosome_plot = list(wrap_plots(linkage_plot, guides = 'keep', nrow = 1) & plot_annotation(title = str_c(chromosome, plot_group, sep = '_'))),
            .groups = 'drop') %>%
  mutate(plot_name = str_c(if_else(only_fdr_sig_clumps,
                                   '../Results/linkage_group_plots_functional_onlySigLinkage/',
                                   '../Results/linkage_group_plots_functional/'), 
                           chromosome, '_', plot_group, '.png'))


walk2(chromosome_plots$plot_name, chromosome_plots$chromosome_plot,
      ~ggsave(.x, plot = .y, height = 15, width = 15))

chromosome_plots %>%
  filter(chromosome == 'Acerv_scaffold_142') %>%
  pull(chromosome_plot)

chromosome_plots %>%
  filter(chromosome == 'Acerv_scaffold_8') %>%
  pull(chromosome_plot)

chromosome_plots$chromosome_plot[[6]]

#### Make Overall/Clump Plots ####
chromosome_plots <- snp_info_joined %>%
  
  # filter(chromosome == 'Acerv_scaffold_115') %>%
  
  expand_grid(tibble(plot_type = c('gene_diagram', 'p', 'es', 'clump'))) %>%
  unnest(clumps) %>%
  mutate(plot_type = if_else(plot_type == 'clump', clumps, plot_type)) %>%
  select(-clumps) %>%
  distinct %>%
  rowwise %>%
  mutate(locus_data = case_when(str_detect(plot_type, 'clump') ~ list(filter(locus_data, clump_id == plot_type)),
                                TRUE ~ list(locus_data)),
         plot = case_when(plot_type == 'gene_diagram' ~ list(create_gene_diagram(gene_data, 
                                                                                 region_start, region_end, 
                                                                                 plot_start, plot_end)),
                          plot_type == 'p' ~ list(make_es_p_plot('p', locus_data, 
                                                                 region_start, region_end,
                                                                 plot_start, plot_end)),
                          TRUE ~ list(make_es_p_plot('es', locus_data, 
                                                     region_start, region_end, 
                                                     plot_start, plot_end)))) %>%
  filter(nrow(locus_data) > 0) %>%
  ungroup %>%
  
  group_by(chromosome, plot_group) %>%
  mutate(position = row_number(),
         max_pos = position != max(position)) %>%
  rowwise %>%
  mutate(plot = if_else(max_pos,
                        list(plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                                          axis.ticks.x = element_blank())),
                        list(plot))) %>%
  
  group_by(chromosome, plot_group) %>%
  summarise(linkage_plot = list(wrap_plots(plot, guides = 'collect', ncol = 1, 
                                           height = c(0.05, rep(0.95 / (n() - 1), (n() - 1))))),
            .groups = 'drop_last') %>%
  
  mutate(keep_legend = row_number() == max(row_number())) %>%
  rowwise %>%
  mutate(linkage_plot = if_else(keep_legend, 
                                list(linkage_plot), 
                                list(linkage_plot & theme(legend.position = 'none')))) %>%
  
  group_by(chromosome) %>%
  summarise(chromosome_plot = list(wrap_plots(linkage_plot, guides = 'keep', nrow = 1) & plot_annotation(title = chromosome)),
            .groups = 'drop') %>%
  mutate(plot_name = str_c('../Results/linkage_group_plots_onlySigLinkage/', 
                           chromosome, '.png'))

chromosome_plots$chromosome_plot[[6]]

walk2(chromosome_plots$plot_name, chromosome_plots$chromosome_plot,
      ~ggsave(.x, plot = .y, height = 15, width = 15))

# for(i in 1:nrow(chromosome_plots)){
#   readline(prompt = str_c("Press [enter] to plot ", chromosome_plots$chromosome[i]))
#   print(chromosome_plots$chromosome_plot[[i]])
# }

#### Output linkage block summary ####


#### Output Gene info for genes in linkage blocks ####