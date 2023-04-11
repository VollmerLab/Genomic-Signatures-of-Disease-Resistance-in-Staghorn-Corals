#### Plot Full Genome & LFMM/PGS functional effects

#### Libraries ####
library(Biostrings)
library(IRanges)
library(bioseq)
library(tidyverse)
library(magrittr)
library(patchwork)

alpha <- 0.05
pgs_p_cutoff <- 0.0001

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


#### Data ####
genome <- readBStringSet('../../Bioinformatics/genome_assembly/k2HybridAssembly.fasta')
new_genome <- readBStringSet('../../Bioinformatics/genome_assembly/k2_genome.fasta')
genome <- genome[names(genome) %in% names(new_genome)]

gene_keep <- microseq::readGFF('../../Bioinformatics/genome_annotation/k2_structuralAnnotations.gff3') %>%
  filter(Type == 'gene') %>%
  mutate(gene_id = str_extract(Attributes, 'Acer_[0-9]+')) %>%
  pull(gene_id)


index_snp_summary <- read_csv('../Results/index_snp_summary.csv', 
                              show_col_types = FALSE)

lfmm_results <- read_csv('../Results/lfmm_results.csv', show_col_types = FALSE) %>%
  filter(pop_group == 'all',
         model_family == 'norm') %>%
  select(-pop_group, -model_family) %>%
  mutate(locus = str_c(chromosome, position, sep = '-'), .after = 'position')  %>%
  filter(chromosome %in% names(genome))

gene_locations <- read_csv('../../Bioinformatics/genome_annotation/r5_gene_id_locations.csv', show_col_types = FALSE) %>%
  filter(gene_type != 'tRNA') %>%
  filter(gene_id %in% gene_keep)
structural_annotations <- read_csv('../../Bioinformatics/genome_annotation/r5_structural_annotation.csv', show_col_types = FALSE) %>%
  filter(gene_id %in% gene_keep)
functional_annotations <- read_rds('../../Bioinformatics/genome_annotation/r5_functional_annotations.rds.gz') %>%
  dplyr::rename(gene_id = qseqid) %>%
  filter(gene_id %in% gene_keep)
codon_sequences <- read_csv_chunked('../../Bioinformatics/genome_annotation/r5_codon_positions.csv.gz', 
                                    callback = DataFrameCallback$new(preprocess_codons),
                                    chunk_size = 100000, col_types = 'cccciiciic') %>%
  filter(gene_id %in% gene_keep)

linkage_blocks <- read_csv('../intermediate_files/lfmm_snp_clumps.csv', show_col_types = FALSE) %>%
  filter(chromosome %in% names(genome))
threshold_linkage_inclusion <- read_csv('../intermediate_files/clumping_thersholding_results.csv', show_col_types = FALSE) %>%
  filter(threshold_p == pgs_p_cutoff) %>%
  select(chromosome, clump_id) %>%
  filter(chromosome %in% names(genome))

snp_data <- read_delim('../../Bioinformatics/variant_calling/18October2022/genolike.mafs.gz', delim = '\t', col_types = 'cicccnnnni') %>%
  select(chromo, position, major, minor, ref) %>%
  dplyr::rename(chromosome = chromo) %>%
  mutate(locus = str_c(chromosome, position, sep = '-'),
         alt = if_else(ref == major, minor, major)) %>%
  filter(chromosome %in% names(genome))

chromosome_groups <- tibble(chromosome = names(genome),
                            length = width(genome)) %>%
  mutate(chrom_group = row_number() %% 2 == 1) %>%
  mutate(concat_start_position = cumsum(lag(length, default = 0L)))

#### Merge all Data ####
fully_joined_snp_data <- full_join(lfmm_results, linkage_blocks,
          by = c('chromosome', 'locus')) %>%
  left_join(mutate(threshold_linkage_inclusion, include_pgs = TRUE), 
            by = c('chromosome', 'clump_id')) %>%
  mutate(include_pgs = if_else(is.na(include_pgs), FALSE, include_pgs)) %>%
  group_by(chromosome) %>%
  # filter(any(include_pgs)) %>%
  ungroup %>%
  
  left_join(chromosome_groups, by = 'chromosome') %>%
  select(-chrom_group) %>%
  nest(locus_info = -c(chromosome, length, concat_start_position)) %>%
  arrange(concat_start_position) %>%
  
  mutate(chrom_group = row_number() %% 2 == 1) %>%
  mutate(concat_start_position = cumsum(lag(length, default = 0L))) %>%
  unnest(locus_info) %>%
  
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
  mutate(mutation_type = if_else(str_detect(mutation_type, 'Non-Coding'), 
                                 'Non-Coding', mutation_type) %>%
           factor(levels = c('Nonsense', 'Non-Synonomous', 'Synonomous', 'Non-Coding')) %>%
           fct_rev,
         coding = mutation_type %in% c('Non-Synonomous', 'Synonomous'))


all_loci_info <- fully_joined_snp_data %>% 
  select(chromosome, position, locus, es, z, p, fdr, clump_id, is_index, include_pgs, ref, alt, gene_id, mutation_type) %>%
  left_join(select(functional_annotations, gene_id, starts_with('swissprot'), blast_details, kegg_id, entap_results) %>%
              unnest(c(blast_details, entap_results), keep_empty = TRUE), 
            by = 'gene_id') 

write_csv(all_loci_info, xzfile('../Results/Table_S2.csv.xz', compression = 9))
# write_csv(all_loci_info, '../Results/Table_S2.csv.gz')

#### Make Plots ####
smoosh_whitespace <- function(data){
  data %>%
    select(-chrom_group, -concat_start_position) %>%
    nest_by(chromosome) %>%
    arrange(str_extract(chromosome, '[0-9]+') %>% as.integer) %>%
    rowwise %>%
    mutate(start_position = min(data$position),
           length = max(data$position) - min(data$position)) %>%
    ungroup %>%
    mutate(chrom_group = row_number() %% 2 == 1) %>%
    mutate(concat_start_position = cumsum(lag(length, default = 0L))) %>%
    select(-length) %>%
    unnest(data) %>%
    mutate(position = position - start_position) 
}

plot_data <- fully_joined_snp_data %>%
  
  group_by(chromosome, clump_id) %>%
  # filter(any(fdr < alpha)) %>%
  ungroup %>%
  
  # inner_join(threshold_linkage_inclusion, by = c('chromosome', 'clump_id')) %>%
  
  smoosh_whitespace %>%
  
  select(chromosome, chrom_group, concat_start_position, 
         position, locus, es, p, fdr, clump_id, is_index, include_pgs, 
         mutation_type) %>%
  
  mutate(position_concat = position + concat_start_position) %>%
  mutate(coloration = case_when(#include_pgs & is_index & (!fdr < alpha) ~ 'linked',
                                fdr < alpha ~ 'significant',
                                chrom_group ~ 'A',
                                !chrom_group ~ 'B'),
         sization = case_when(include_pgs & is_index ~ 'interesting',
                              fdr < alpha ~ 'interesting_2',
                              TRUE ~ 'not')) %>%
  arrange(coloration) %>%
  mutate(chrom_number = str_extract(chromosome, '[0-9]+') %>% as.integer) %>%
  filter(chrom_number <= max(chrom_number[include_pgs & is_index])) %>%
  mutate(position_concat2 = position_concat + (1e6 * chrom_number))


make_chromosome_linkage_diagram <- function(data){
  chromosome_data <- plot_data %>%
    group_by(chromosome) %>%
    filter(position_concat == min(position_concat) | position_concat == max(position_concat)) %>%
    select(chromosome, position_concat) %>%
    distinct %>%
    mutate(which_end = if_else(position_concat == min(position_concat), 'min', 'max')) %>%
    pivot_wider(names_from = which_end, values_from = position_concat) %>%
    ungroup %>%
    
    mutate(chrom_mid = (max + min) / 2,
           chromosome = str_extract(chromosome, '[0-9]+')) %>%
    arrange(as.integer(chromosome))
  
  linkage_data <- plot_data %>%
    group_by(chromosome, clump_id) %>%
    filter(position_concat == min(position_concat) | position_concat == max(position_concat)) %>%
    select(chromosome, clump_id, position_concat) %>%
    mutate(which_end = if_else(position_concat == min(position_concat), 'min', 'max')) %>%
    distinct %>%
    pivot_wider(names_from = which_end, values_from = position_concat) %>%
    ungroup %>%
    
    group_by(chromosome) %>%
    mutate(n_linkage = n())  %>%
    
    
    mutate(n_row_within = 2 + n_linkage + (n_linkage - 1),
           which_group = row_number(),
           row_within = 1 + which_group + lag(which_group, default = 0L)) %>%
    mutate(y_width = 0.5 / n_row_within) %>%
    mutate(y_start = 0.5 - (y_width * (row_within - 1)),
           y_end = y_start - y_width) %>%
    ungroup
    
  
  
  chromosome_data %>%
    
    ggplot(aes(xmin = (min + 100) / 1e6, 
               xmax = (max - 100) / 1e6, 
               ymin = 0, ymax = 0.5)) +
    
    geom_rect(fill = 'grey50') +
    geom_rect(data = linkage_data, inherit.aes = FALSE,
              aes(xmin = (min + 100) / 1e6, xmax = (max - 100) / 1e6, 
                  ymin = y_start, ymax = y_end), fill = 'grey5') +
    geom_label(aes(x = chrom_mid / 1e6, y = 0.25, label = chromosome),
              hjust = 0.5, vjust = 0.5, colour = 'white', fill = 'grey50') +
    theme_void() +
    theme(axis.text = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text())
  
}

chromosome_diagram <- make_chromosome_linkage_diagram(chromosome_data)


max_first  <- max(plot_data$position_concat2)   # Specify max of first y axis
max_second <- max(plot_data$position_concat) # Specify max of second y axis
min_first  <- min(plot_data$position_concat2)   # Specify min of first y axis
min_second <- min(plot_data$position_concat) # Specify min of second y axis

# scale and shift variables calculated based on desired mins and maxes
scale = (max_second - min_second)/(max_first - min_first)
shift = min_first - min_second

# Function to scale secondary axis
scale_function <- function(x, scale, shift){
  return ((x)*scale - shift)
}

# Function to scale secondary variable values
inv_scale_function <- function(x, scale, shift){
  return ((x + shift)/scale)
}


#https://github.com/tidyverse/ggplot2/issues/2075

colour_choice <- c(wesanderson::wes_palette("Zissou1", 9, type = "continuous")[c(9,1)],
                   'gray25', 'gray75') %>%
  set_names(c('significant', 'linked', 'A', 'B'))

p_values_genome <- plot_data %>%
  # sample_frac(0.05) %>%
  
  ggplot(aes(x = position_concat2 / 1e6, y = -log(p, base = 10), 
             colour = coloration, size = sization,
             shape = include_pgs & is_index)) +
  geom_point() +
  guides(x = guide_axis(position = 'none')) +
  scale_y_continuous(expand = c(0, 0), limits = c(NA, 8)) +
  scale_x_continuous(labels = scales::comma_format(), 
                     sec.axis = sec_axis(~scale_function(., scale, shift), 
                                         guide = guide_axis(position = 'bottom')),
                     expand = c(0, 0)) +
  scale_colour_manual(values = c('significant' = 'red', 'A' = 'grey25', 'B' = 'grey75')) +
  scale_size_manual(values = c('interesting' = 2.5, 'interesting_2' = 1, 'not' = 0.25)) +
  scale_shape_manual(values = c('TRUE' = 'cross', 'FALSE' = 'circle')) +
  
  # facet_grid(~chrom_number, scales = 'free_x', space = 'free_x') +
  coord_cartesian(clip = "off") +
  guides(size = 'none',
         colour = 'none',
         shape = 'none') +
  labs(x = 'Position (MB)',
       # y = '-log10(p)',
       y = bquote(-log[10](p)),
       shape = 'Signficant FDR') +
  theme_classic() +
  theme(axis.text = element_text(face = 'plain', colour = 'black', size = 24),
        axis.title = element_text(face = 'plain', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(5, 15, 5, 5))


genome_es <- plot_data %>%
  
  # filter(coloration == 'index') %>%
  
  # sample_frac(0.01) %>%
  ggplot(aes(x = position_concat2 / 1e6, y = es, 
             colour = coloration, size = sization,
             shape = include_pgs & is_index)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  guides(x = guide_axis(position = 'none')) +
  scale_x_continuous(labels = scales::comma_format(), 
                     sec.axis = sec_axis(~scale_function(., scale, shift), 
                                         guide = guide_axis(position = 'bottom')),
                     expand = c(0, 0)) +
  scale_colour_manual(values = c('significant' = 'red', 'A' = 'grey25', 'B' = 'grey75')) +
  scale_size_manual(values = c('interesting' = 2.5, 'interesting_2' = 1, 'not' = 0.25)) +
  scale_shape_manual(values = c('TRUE' = 'cross', 'FALSE' = 'circle')) +
  
  guides(size = 'none',
         colour = 'none',
         shape = 'none') +
  labs(x = 'Position (MB)',
       y = 'Effect Size',
       shape = 'Signficant FDR') +
  theme_classic() +
  theme(axis.text = element_text(face = 'plain', colour = 'black', size = 24),
        axis.title = element_text(face = 'plain', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 1),
        plot.margin = margin(5, 15, 5, 5))

genome_plot <- (p_values_genome / genome_es) &
  plot_layout(ncol = 1) &
  plot_annotation(tag_levels = 'A')  &
  theme(plot.tag = element_text(size = 24, colour = 'black'))
ggsave('../Results/genome_lfmm_plot.png', plot = genome_plot, height = 10, width = 15)

#Effect Size all N across genome color by if in PGS
genome_functional <- plot_data %>%
  
  group_by(chromosome, clump_id) %>%
  mutate(fdr = min(fdr)) %>%
  ungroup %>%
  
  mutate(coloration = case_when(fdr < alpha ~ 'significant',
                                include_pgs ~ 'pgs_clump',
                                chrom_group ~ 'A',
                                !chrom_group ~ 'B'),
         sization = case_when(include_pgs & is_index ~ 'interesting',
                              fdr < alpha ~ 'interesting',
                              include_pgs ~ 'interesting',
                              TRUE ~ 'not')) %>%
  arrange(coloration) %>%
  filter(mutation_type == 'Non-Synonomous') %>%
  
  # sample_frac(0.01) %>%
  ggplot(aes(x = position_concat / 1e6, y = es, 
             colour = coloration, size = sization,
             shape = include_pgs & is_index)) +
  geom_point() +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_colour_manual(values = c('significant' = 'red', 'pgs_clump' = 'orange', 
                                 'A' = 'grey25', 'B' = 'grey75')) +
  scale_size_manual(values = c('interesting' = 2.5, 'not' = 0.25)) +
  scale_shape_manual(values = c('TRUE' = 'cross', 'FALSE' = 'circle')) +
  
  guides(size = 'none',
         colour = 'none',
         shape = 'none') +
  labs(x = 'Position (MB)',
       y = 'ES',
       shape = 'Signficant FDR') +
  theme_classic() +
  theme(axis.text = element_text(face = 'plain', colour = 'black', size = 24),
        axis.title = element_text(face = 'plain', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA))


#### Other Plots ####


interesting_plot <- fully_joined_snp_data %>%
  select(chromosome, chrom_group, concat_start_position, position, locus, es, p, fdr, clump_id, is_index, include_pgs, 
         mutation_type) %>%
  
  # filter(include_pgs) %>%
  # filter(mutation_type == 'Non-Synonomous') %>%
  
  mutate(position_concat = position + concat_start_position) %>%
  mutate(coloration = case_when(fdr < alpha ~ 'significant',
                                chrom_group ~ 'A',
                                !chrom_group ~ 'B'),
         sization = case_when(include_pgs & is_index ~ 'interesting',
                              fdr < alpha ~ 'interesting',
                              TRUE ~ 'not')) %>%
  arrange(coloration) %>%
  # filter(include_pgs & is_index) %>%
  filter(sization == 'interesting') %>%
  
  # filter(coloration == 'index') %>%
  
  # sample_frac(0.01) %>%
  ggplot(aes(x = position_concat / 1e6, y = es, 
             colour = coloration, size = sization,
             shape = include_pgs & is_index)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point() +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_colour_manual(values = c('significant' = 'red', 'A' = 'grey25', 'B' = 'grey75')) +
  scale_size_manual(values = c('interesting' = 2.5, 'not' = 0.25)) +
  scale_shape_manual(values = c('TRUE' = 'cross', 'FALSE' = 'circle')) +
  
  guides(size = 'none',
         colour = 'none',
         shape = 'none') +
  labs(x = 'Position (MB)',
       y = 'Effect Size',
       shape = 'Signficant FDR') +
  theme_classic() +
  theme(axis.text = element_text(face = 'plain', colour = 'black', size = 24),
        axis.title = element_text(face = 'plain', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA))


functional_plot <- tmp %>%
  select(chromosome, chrom_group, concat_start_position, position, locus, es, p, fdr, clump_id, is_index, include_pgs, 
         mutation_type) %>%
  
  # filter(include_pgs) %>%
  # filter(mutation_type == 'Non-Synonomous') %>%
  
  mutate(position_concat = position + concat_start_position) %>%
  mutate(coloration = case_when(fdr < alpha ~ 'significant',
                                chrom_group ~ 'A',
                                !chrom_group ~ 'B'),
         sization = case_when(include_pgs & is_index ~ 'interesting',
                              fdr < alpha ~ 'interesting',
                              TRUE ~ 'not')) %>%
  arrange(coloration) %>%
  # filter(include_pgs & is_index) %>%
  group_by(chromosome, clump_id) %>%
  filter(any(sization == 'interesting')) %>%
  ungroup %>%
  filter(mutation_type == 'Non-Synonomous' | (include_pgs & is_index)) %>%
  
  # filter(coloration == 'index') %>%
  
  # sample_frac(0.01) %>%
  ggplot(aes(x = position_concat / 1e6, y = es, 
             colour = coloration, 
             shape = include_pgs & is_index)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point() +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_colour_manual(values = c('significant' = 'red', 'A' = 'grey25', 'B' = 'grey75')) +
  # scale_size_manual(values = c('interesting' = 2.5, 'not' = 0.25)) +
  scale_shape_manual(values = c('TRUE' = 'cross', 'FALSE' = 'circle')) +
  
  guides(size = 'none',
         colour = 'none',
         shape = 'none') +
  labs(x = 'Position (MB)',
       y = 'Effect Size',
       shape = 'Signficant FDR') +
  theme_classic() +
  theme(axis.text = element_text(face = 'plain', colour = 'black', size = 24),
        axis.title = element_text(face = 'plain', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA))



(p_values_genome / interesting_plot / functional_plot) &
  plot_layout(ncol = 1)


tmp %>%
  select(chromosome, locus, es, clump_id, is_index, include_pgs, mutation_type) %>%
  distinct %>%
  
  filter(include_pgs) %>%
  
  group_by(chromosome, clump_id) %>%
  mutate(index_effect = es[is_index]) %>%
  filter(!is_index) %>%
  ungroup %>%
  # filter(chromosome == 'Acerv_scaffold_109', clump_id == 'clump_1') %>%
  summarise(cor(index_effect, es))


tmp %>%
  select(chromosome, locus, es, clump_id, is_index, include_pgs, mutation_type) %>%
  distinct %>%
  
  filter(include_pgs) %>%
  
  group_by(chromosome, clump_id) %>%
  mutate(index_effect = es[is_index]) %>%
  filter(!is_index) %>%
  ungroup %>%
  filter(mutation_type == 'Non-Synonomous') %>%
  
  ggplot(aes(x = index_effect, y = es)) +
  geom_point()


tst_dat <- tmp %>%
  select(chromosome, chrom_group, concat_start_position, position, locus, es, p, fdr, clump_id, is_index, include_pgs, 
         mutation_type) %>%
  
  # filter(include_pgs) %>%
  # filter(mutation_type == 'Non-Synonomous') %>%
  
  mutate(position_concat = position + concat_start_position) %>%
  mutate(coloration = case_when(fdr < alpha ~ 'significant',
                                chrom_group ~ 'A',
                                !chrom_group ~ 'B'),
         sization = case_when(include_pgs & is_index ~ 'interesting',
                              fdr < alpha ~ 'interesting',
                              TRUE ~ 'not')) %>%
  arrange(coloration) %>%
  # filter(include_pgs & is_index) %>%
  group_by(chromosome, clump_id) %>%
  filter(any(sization == 'interesting')) %>%
  mutate(dist_index = abs(es[is_index]) - abs(es)) %>%
  ungroup 


tst_dat %>%
  # filter(mutation_type == 'Non-Synonomous') %>%
  ggplot(aes(x = dist_index, y = interaction(chromosome, clump_id))) +
  stat_summary()



tmp %>%
  select(chromosome, chrom_group, concat_start_position, position, locus, es, p, fdr, clump_id, is_index, include_pgs, 
         mutation_type) %>%
  
  # filter(include_pgs) %>%
  # filter(mutation_type == 'Non-Synonomous') %>%
  
  mutate(position_concat = position + concat_start_position) %>%
  mutate(coloration = case_when(fdr < alpha ~ 'significant',
                                chrom_group ~ 'A',
                                !chrom_group ~ 'B'),
         sization = case_when(include_pgs & is_index ~ 'interesting',
                              fdr < alpha ~ 'interesting',
                              TRUE ~ 'not')) %>%
  arrange(coloration) %>%
  filter(include_pgs & is_index) %>%
  
  sample_frac(0.01) %>%
  ggplot(aes(x = position_concat / 1e6, y = -log(p, base = 10), 
             colour = coloration, size = sization,
             shape = include_pgs & is_index)) +
  geom_point() +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_colour_manual(values = c('significant' = 'red', 'A' = 'grey25', 'B' = 'grey75')) +
  scale_size_manual(values = c('interesting' = 2.5, 'not' = 0.25)) +
  scale_shape_manual(values = c('TRUE' = 'cross', 'FALSE' = 'circle')) +
  
  guides(size = 'none',
         colour = 'none',
         shape = 'none') +
  labs(x = 'Position (MB)',
       y = 'Effect Size',
       shape = 'Signficant FDR') +
  theme_classic() +
  theme(axis.text = element_text(face = 'plain', colour = 'black', size = 24),
        axis.title = element_text(face = 'plain', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA))






  
  
snp_info_joined %>%
  select(chromosome, plot_group, plot_start, plot_end, locus_data) %>%
  mutate(length = plot_end - plot_start) %>%
  arrange(str_extract(chromosome, '[0-9]+') %>% as.integer,
          plot_start) %>%
  
  mutate(chrom_group = row_number() %% 2 == 1) %>%
  mutate(concat_start_position = cumsum(lag(length, default = 0L))) %>%
  select(-plot_start, -plot_end) %>%
  
  rowwise %>%
  mutate(locus_data = list(mutate(locus_data, position = position - min(position)))) %>%
  ungroup %>%
  
  unnest(locus_data) %>%
  select(chromosome, plot_group, locus, chrom_group, concat_start_position, 
         position, es, clump_id, mutation_type) %>%
  mutate(position_concat = position + concat_start_position) %>%
  
  left_join(select(tmp, chromosome, clump_id, include_pgs) %>%
               filter(include_pgs) %>%
               distinct,
             by = c('chromosome', 'clump_id')) %>%
  filter(mutation_type == 'Non-Synonomous') %>%
  mutate(include_pgs = if_else(is.na(include_pgs), FALSE, TRUE)) %>%
  mutate(coloration = case_when(include_pgs ~ 'aaaindex',
                                chrom_group ~ 'A',
                                !chrom_group ~ 'B')) %>%
  arrange(coloration) %>%
  # filter(position_concat / 1e6 < 3) %>%
  
  ggplot(aes(x = position_concat / 1e6, y = es, colour = coloration)) +
  geom_point() +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_colour_manual(values = c('aaaindex' = 'red', 'A' = 'grey25', 'B' = 'grey75')) +
  guides(size = 'none',
         colour = 'none',
         shape = guide_legend(override.aes = list(size = 2))) +
  labs(x = 'Position (MB)',
       y = 'Effect Size') +
  theme_classic() +
  theme(axis.text = element_text(face = 'plain', colour = 'black', size = 24),
        axis.title = element_text(face = 'plain', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA))
  



snp_info_joined %>%
  rowwise %>%
  mutate(locus_data = list(group_by(locus_data, clump_id) %>%
                             filter(any(fdr < alpha)) %>%
                             ungroup)) %>%
  filter(nrow(locus_data) != 0) %>%
  ungroup %>%
  select(chromosome, plot_group, plot_start, plot_end, locus_data) %>%
  mutate(length = plot_end - plot_start) %>%
  arrange(str_extract(chromosome, '[0-9]+') %>% as.integer,
          plot_start) %>%
  
  mutate(chrom_group = row_number() %% 2 == 1) %>%
  mutate(concat_start_position = cumsum(lag(length, default = 0L))) %>%
  select(-plot_start, -plot_end) %>%
  
  rowwise %>%
  mutate(locus_data = list(mutate(locus_data, position = position - min(position)))) %>%
  ungroup %>%
  
  unnest(locus_data) %>%
  select(chromosome, plot_group, locus, chrom_group, concat_start_position, 
         position, es, clump_id, mutation_type) %>%
  mutate(position_concat = position + concat_start_position) %>%
  
  left_join(select(tmp, chromosome, clump_id, include_pgs) %>%
              filter(include_pgs) %>%
              distinct,
            by = c('chromosome', 'clump_id')) %>%
  filter(mutation_type == 'Non-Synonomous') %>%
  mutate(include_pgs = if_else(is.na(include_pgs), FALSE, TRUE)) %>%
  mutate(coloration = case_when(include_pgs ~ 'aaaindex',
                                chrom_group ~ 'A',
                                !chrom_group ~ 'B')) %>%
  arrange(coloration) %>%
  # filter(position_concat / 1e6 < 3) %>%
  
  ggplot(aes(x = position_concat / 1e6, y = es, colour = chromosome)) +
  geom_point() +
  scale_x_continuous(labels = scales::comma_format()) +
  # scale_colour_manual(values = c('aaaindex' = 'red', 'A' = 'grey25', 'B' = 'grey75')) +
  guides(size = 'none',
         # colour = 'none',
         shape = guide_legend(override.aes = list(size = 2))) +
  labs(x = 'Position (MB)',
       y = 'Effect Size') +
  theme_classic() +
  theme(axis.text = element_text(face = 'plain', colour = 'black', size = 24),
        axis.title = element_text(face = 'plain', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA))







snp_info_joined %>%
  select(chromosome, plot_group, clumps, locus_data) %>% 
  # filter(chromosome %in% c('Acerv_scaffold_8', 'Acerv_scaffold_142')) %>%
  unnest(locus_data) %>%
  filter(mutation_type == 'Non-Synonomous') %>%
  
  rowwise %>%
  mutate(in_pgs = clump_id %in% clumps) %>%
  ungroup %>%
  
  arrange(-1 * log(p, base = 10)) %>%
  
  ggplot(aes(x = position / position_divisor, y = es)) +
  
  geom_hline(yintercept = 0, linetype = 'dashed') +
  
  geom_point(aes(colour = in_pgs,
                 size = coding,
                 shape = mutation_type)) +
  
  # scale_x_continuous(limits = range_of_interest / position_divisor) +
  scale_size_manual(values = c('TRUE' = 2, 'FALSE' = 0.5), drop = FALSE) +
  scale_shape_manual(values = c('Non-Coding' = 'circle',
                                'Synonomous' = 'square',
                                'Non-Synonomous' = 'triangle',
                                'Nonsense' = 'cross'), 
                     drop = FALSE) +
  guides(size = 'none',
         shape = guide_legend(override.aes = list(size = 3))) +
  facet_grid(~chromosome + plot_group, scales = 'free_x', space = 'free_x') +
  labs(x = 'Position (MB)',
       y = 'Effect Size',
       colour = 'In PGS',
       # colour = bquote(-log[10](p)),
       shape = 'Mutation Type') +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 16),
        panel.border = element_rect(colour = 'black', fill = NA),
        axis.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 16),
        legend.title = element_text(colour = 'black', size = 14))
  
  


