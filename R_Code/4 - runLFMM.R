#### Libraries ####
library(tidyverse)
library(magrittr)
library(adegenet)
library(poppr)
library(LEA)
library(patchwork)

#### Function ####
write_phenotype <- function(dat, in_name){
  out <- write.env(select(dat, disease_resistance), str_replace(in_name, '\\.snmfProject$', '.env'))
}

impute_missing <- function(in_lfmm, in_structure, K_choice){
  if(!file.exists(str_c(in_lfmm, '_imputed.lfmm'))){
    structure_out <- load.snmfProject(in_structure)
    imputation_out <- impute(structure_out, in_lfmm, method = 'mode', 
                             K = K_choice, run = 1)
  }
  
  imputation_out <- str_c(in_lfmm, '_imputed.lfmm')
}

#### Read in Data ####
meta_data <- read_csv('../intermediate_files/preprocessed_metadata.csv', 
                      show_col_types = FALSE) %>%
  filter(!is.na(disease_resistance)) %>%
  mutate(resist_norm = qqnorm(disease_resistance, plot.it = FALSE)$x)

intial_data <- expand_grid(pop_group = c('acerv'),
            model_family = c('pois_norm_binom'),
            only_unlinked = c(TRUE, FALSE)) %>%
  
  filter(!only_unlinked) %>%
  
  mutate(genomic_file = str_c('../intermediate_files/', if_else(only_unlinked, 'unlinked_', ''), 'preprocessed_', pop_group, '.rds'),
         structure_file = str_c('../lea_projects/', if_else(only_unlinked, 'unlinked_', ''), 'preprocessed_', pop_group, '.snmfProject'),
         lfmm_file = str_c('../lea_projects/', if_else(only_unlinked, 'unlinked_', ''), 'preprocessed_', pop_group, '.lfmm')) %>%
  mutate(K = c(2)) %>%
  rowwise %>%
  mutate(metadata = case_when(pop_group == 'acerv' ~ list(meta_data),
                               pop_group == 'acerv_panama' ~ list(meta_data %>% filter(location == 'Panama')),
                               pop_group == 'acerv_florida' ~ list(meta_data %>% filter(location == 'Florida')),
                               TRUE ~ list(NA))) %>%
  mutate(phenotype_file = write_phenotype(metadata, structure_file)) %>%
  select(-metadata) %>%
  ungroup

#### Impute Missing Data for LFMM ####
setwd('../lea_projects')

imputed_data <- intial_data %>%
  mutate(second_model = 'binom',
         model_family = str_remove(model_family, '_binom')) %>%
  pivot_longer(cols = c(model_family, second_model),
               names_to = 'tmp',
               values_to = 'model_family') %>%
  select(-tmp) %>%
  filter(str_detect(model_family, 'norm')) %>%
  mutate(structure_file = if_else(model_family == 'binom', str_replace(structure_file, '\\.snmfProject$', '_binom.snmfProject'), structure_file),
         lfmm_file = if_else(model_family == 'binom', str_replace(lfmm_file, '\\.lfmm$', '_binom.lfmm'), lfmm_file)) %>%
  rowwise %>%
  mutate(impute_file = impute_missing(lfmm_file, structure_file, K))

#### Run LFMM ####
if(file.exists('../Results/lfmm_results.csv')){
  lfmm_results <- read_csv('../Results/lfmm_results.csv')
} else {
  library(multidplyr)
  cluster <- new_cluster(parallel::detectCores() - 1)
  cluster_library(cluster, c('LEA', 'readr', 'adegenet',
                             'dplyr', 'magrittr', 'purrr'))
  
  lfmm_results <- imputed_data %>%
    select(pop_group, model_family, impute_file, phenotype_file, K, genomic_file) %>%
    rowwise %>%
    partition(cluster) %>%
    mutate(initial_lfmm = list(lfmm2(impute_file, phenotype_file, K = K, effect.sizes = TRUE))) %>%
    collect %>%
    ungroup() %>%
    separate(model_family, into = c("A", 'B'), sep = '_', fill = 'right') %>%
    pivot_longer(cols = c('A', 'B'),
                 names_to = 'tmp',
                 values_to = 'model_family') %>%
    select(-tmp) %>%
    filter(!is.na(model_family)) %>%
    rowwise(pop_group, model_family) %>%
    filter(model_family == 'norm') %>%
    
    partition(cluster) %>%
    mutate(genomic_data = list(read_rds(genomic_file))) %>%
    summarise(chromosome = genomic_data@chromosome,
              position = genomic_data@position,
              effect_size = as.numeric(initial_lfmm@B),
              lfmm2.test(initial_lfmm, 
                         impute_file,
                         phenotype_file,
                         full = FALSE,
                         genomic.control = TRUE,
                         linear = FALSE,
                         family = case_when(model_family == 'pois' ~ 'poisson',
                                            model_family == 'norm' ~ 'gaussian',
                                            model_family == 'binom' ~ 'binomial')) %>%
                extract(2:1) %>%
                map(as.numeric) %>%
                bind_cols() %>%
                rename(z = zscores,
                       p = pvalues) %>%
                mutate(fdr = p.adjust(p, method = 'fdr'))) %>%
    collect %>%
    ungroup %>%
    rename(es = effect_size) %>%
    mutate(pop_group = str_replace_all(pop_group, c('^acerv$' = 'all', 'acerv_panama' = 'panama', 'acerv_florida' = 'florida')))
  write_csv(lfmm_results, '../Results/lfmm_results.csv')
}

#### Diagnostic Plots ####
ks_tests <- lfmm_results %>%
  group_by(pop_group, model_family) %>%
  summarise(broom::tidy(ks.test(p, "punif", 0, 1)),
            .groups = 'drop') %>%
  mutate(annotation_text = str_c('D = ', round(statistic, 3), '\np = ', scales::pvalue(p.value)),
         x = 0.5, y = -Inf)

lfmm_results %>%
  ggplot(aes(x = p)) +
  geom_histogram(bins = 50) +
  geom_text(data = ks_tests, aes(x = x, y = y, label = annotation_text), colour = 'white', vjust = -1) +
  labs(x = 'p-value', 
       y = 'Count') +
  facet_grid(pop_group ~ model_family) +
  theme_classic() 
ggsave('../Results/lfmm_diagnostic_histograms.png', height = 8, width = 8)

lfmm_results %>%
  ggplot(aes(sample = p)) +
  stat_qq_line(distribution = qunif, dparams = list(min = 0, max = 1), size = 2, colour = 'red') +
  stat_qq(distribution = qunif, dparams = list(min = 0, max = 1), size = 0.01, alpha = 0.3) +
  facet_grid(pop_group ~ model_family) +
  theme_classic()
ggsave('../Results/lfmm_diagnostic_qqPlot.png', height = 8, width = 8)
