#### After running ANGSD & PCAngsd output plots & csv #### 

if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  working_dir <- args[1]
  message(args)
} else {
  working_dir <- '/scratch/j.selwyn/variant_calling/k2_18-October-2022/genotyping/subset'
}

setwd(working_dir)
# options(bitmapType="cairo")

#suppressMessages(suppressWarnings(library(tidyverse)))
# suppressMessages(suppressWarnings(library(Cairo)))
suppressMessages(suppressWarnings(library(magrittr)))
# suppressMessages(suppressWarnings(library(ggtree)))
# suppressMessages(suppressWarnings(library(treeio)))
# suppressMessages(suppressWarnings(library(tidytree)))

suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(tibble)))
suppressMessages(suppressWarnings(library(ggplot2)))

#### Read In Data ####

sample_info <- read_delim("bam.list", col_names = 'ID', delim = '\t', show_col_types = FALSE) %>%
  mutate(ID =  str_remove(ID, dirname(ID)) %>% str_remove('/') %>% str_remove('-fp-.*bam$')) %>% 
  mutate(species = str_extract(ID, 'A[a-z]'),
         region = str_extract(ID, '[FP][A-Z]'),
         genotype = str_extract(ID, '[a-zA-Z0-9]+$'),
         
         species = if_else(is.na(species), str_extract(ID, 'acerv|apalm|aprol') %>% 
                             str_remove_all('erv|lm|ol') %>% str_to_sentence(),
                           species)) %>% 
  mutate(reef = if_else(region == 'PA', str_extract(genotype, '[a-zA-Z]+'), 'FL')) %>%
  mutate(id_number = row_number())

pca_results <- read_delim("output.pcangsd.cov", col_names = FALSE, delim = ' ', show_col_types = FALSE) %>%
  as.matrix %>% 
  eigen

admix_files <- list.files(pattern = 'output.pcangsd.admix.*Q')
admix_info <- read_delim(admix_files, col_names = FALSE, delim = ' ', show_col_types = FALSE) %>%
  rename_with(~str_replace(., 'X', 'Cluster')) %>%
  rowwise() %>%
  mutate(admix_assignment = which.max(c_across(starts_with('Cluster'))) %>% 
           str_c('Cluster', .)) %>%
  mutate(max = max(c_across(starts_with('Cluster'))),
         match = which.max(c_across(starts_with('Cluster')))) %>%
  ungroup

#### Convert Genotype File to pair of CSV files with sample names ####
geno_probs <- read_delim('genolike.geno.gz', 
                         col_names = c('contig', 'position', 'ref', 'alt', 
                                       c(rbind(str_c(sample_info$ID, 'geno', sep = '_'), 
                                               str_c(sample_info$ID, 'prob', sep = '_')))), 
                         col_types = str_c(c('cicc', 
                                             rep('id', length(sample_info$ID))), 
                                           collapse = ''),
                         show_col_types = FALSE,
                         delim = '\t',
                         progress = show_progress())

#### If the data is subset then filter genotypes/bamlist ####
if(file.exists('subset_data.txt')){
  kept_samples <- read_delim('subset_data.txt', delim = '\t', col_names = 'keep', col_types = 'l') %>%
    pull(keep) %>%
    which
  
  sample_info <- slice(sample_info, kept_samples) %>%
    mutate(id_number = row_number())
  geno_probs <- select(geno_probs, c(contig, position, ref, alt, contains(sample_info$ID)))
}

if(file.exists('subset_loci.txt')){
  kept_loci <- read_delim('subset_loci.txt', delim = '\t', col_names = 'keep', col_types = 'l') %>%
    pull(keep) %>%
    which
  
  geno_probs <- slice(geno_probs, kept_loci)
}

genotypes <- geno_probs %>%
  select(contig, position, ref, alt, ends_with('geno')) %>%
  rename_with(~str_remove(., '_geno')) %>%
  mutate(across(c(where(is.integer), -position), ~na_if(., -1))) %T>%
  write_csv('genotypes.csv')

genotype_probabilities <- geno_probs %>%
  select(contig, position, ref, alt, ends_with('prob')) %>%
  rename_with(~str_remove(., '_prob')) %>%
  mutate(across(c(where(is.numeric)), ~na_if(., 0))) %T>%
  write_csv('genotype_probabilities.csv')

message('Number of Loci = ', scales::comma(nrow(genotypes)))
message('Number of Individuals = ', scales::comma(nrow(sample_info)))

#### Basic Individual/locus Summaries of missing data ####
percent_missing_data_ind <- select(genotypes, where(is.integer), -position) %>%
  is.na %>%
  colSums() %>%
  divide_by(nrow(genotypes)) %>%
  enframe(name = 'ID',
          value = 'percent_missing_loci') %>%
  arrange(-percent_missing_loci) %T>%
  write_csv('individual_missing_data.csv')

missing_inds <- filter(percent_missing_data_ind, percent_missing_loci > 0.25) %>%
  pull(ID) 

missing_inds_txt <- missing_inds %>%
  c('Individuals with more\nthan 25% missing loci', length(missing_inds), '\n', .) %>%
  str_c(collapse = '\n')

individual_missing_plot <- percent_missing_data_ind %>%
  ggplot(aes(x = percent_missing_loci, fill = percent_missing_loci > 0.25)) +
  geom_histogram(bins = floor(nrow(sample_info) / 2), show.legend = FALSE) +
  scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::comma_format(accuracy = 1)) +
  labs(y = 'Number of Individuals',
       x = 'Percent of Loci Missing in Individual') +
  annotate(geom = 'text', x = 0.5, y = 30, label = missing_inds_txt) +
  theme_classic()
ggsave('individual_missing_plot.png', plot = individual_missing_plot, height = 7, width = 7)

loci_missing_plot <- select(genotypes, where(is.integer), -position) %>%
  is.na %>%
  rowSums() %>%
  divide_by(ncol(genotypes) - 4) %>%
  tibble(percent_missing = .) %>%
  ggplot(aes(x = percent_missing)) +
  geom_histogram(bins = 100) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::comma_format(accuracy = 1)) +
  labs(y = 'Number of Loci',
       x = 'Percent of Individuals Missing Locus') +
  theme_classic()
ggsave('loci_missing_plot.png', plot = loci_missing_plot, height = 7, width = 7)


loci_missing_plot_without_poor_individuals <- select(genotypes, where(is.integer), -position) %>%
  select(-all_of(missing_inds)) %>%
  is.na %>%
  rowSums() %>%
  divide_by(ncol(genotypes) - 4) %>%
  tibble(percent_missing = .) %>%
  ggplot(aes(x = percent_missing)) +
  geom_histogram(bins = 100) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::comma_format(accuracy = 1)) +
  labs(y = 'Number of Loci',
       x = 'Percent of Individuals Missing Locus') +
  theme_classic()
ggsave('loci_missing_plot_without_poor_individuals.png', plot = loci_missing_plot_without_poor_individuals, 
       height = 7, width = 7)

#### Join Data ####
full_data <- bind_cols(
  sample_info,
  admix_info,
  as_tibble(pca_results$vectors, .name_repair = ~str_c('PC', 1:ncol(pca_results$vectors)))
)
write_csv(full_data, 'sample_structure_pca_data.csv')

#### Plot PCA plot ####
var_explained <- pca_results$values / sum(pca_results$values)
varExplained_plot <- tibble(pct_explained = var_explained,
       pc = 1:length(pca_results$values)) %>%
  ggplot(aes(x = pc, y = pct_explained)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = 'Percent Genetic Variance Explained',
       y = 'Principle Component') +
  theme_classic()
ggsave('varExplained_plot.png', plot = varExplained_plot, height = 7, width = 7)

pca_plot <- full_data %>%
  mutate(region = if_else(is.na(region), 'baum', region),
         reef = if_else(is.na(reef), 'baum', reef)) %>%
  
  ggplot(aes(x = PC1, y = PC2, colour = species, shape = region, fill = reef)) +
  geom_point() +
  scale_shape_manual(values = c('baum' = 20, 'FL' = 21, 'PA' = 22)) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = str_c('PC1 (', round(var_explained[1]*100, 1),'%)'),
       y = str_c('PC2 (', round(var_explained[2]*100, 1),'%)')) +
  theme_classic()
ggsave('pca_plot.png', plot = pca_plot, height = 7, width = 7)

pca_cluster <- full_data %>%
  ggplot(aes(x = PC1, y = PC2, colour = admix_assignment)) +
  geom_point() +
  labs(x = str_c('PC1 (', round(var_explained[1]*100, 1),'%)'),
       y = str_c('PC2 (', round(var_explained[2]*100, 1),'%)')) +
  theme_classic()
ggsave('pca_plotClustering.png', plot = pca_cluster, height = 7, width = 7)

#### Investigate Admixture ####
structure_plot <- full_data %>%
  select(-starts_with('PC')) %>%
  mutate(region = if_else(is.na(region), 'Baum', region)) %>%
  group_by(region) %>%
  arrange(match, 
          -max) %>%
  mutate(number = row_number()) %>%
  ungroup %>%
  pivot_longer(cols = starts_with('Cluster')) %>%
  ggplot(aes(x = number, y = value, fill = name)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid( ~ region, scales = 'free_x') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL,
       y = 'Assignment Probability',
       fill = NULL) +
  theme_classic() +
  theme(legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "vertical",
        axis.text = element_text(colour = 'black'),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour = 'black', fill = 'transparent'))
ggsave('structure_plot.png', plot = structure_plot, height = 7, width = 7)

#### Plot Tree ####
tree_data <- read.newick('output.pcangsd.tree') %>%
  as.treedata() %>%
  as_tibble() %>%
  left_join(mutate(sample_info, id_number = as.character(id_number)),
            by = c('label' = 'id_number')) %>%
  select(-label) %>%
  rename(label = ID) %>%
  as.treedata() 

tree_plot <- tree_data %>%
  ggtree() +
  # layout_dendrogram() +
  geom_tippoint(aes(colour = reef)) +
  geom_tiplab(size = 2) +
  geom_treescale(fontsize = 6, linesize = 1, offset = 1) +
  theme(legend.position = 'bottom') +
  scale_x_continuous(expand = c(0,0.1))
ggsave('tree_plot.png', plot = tree_plot, height = 7, width = 7)
