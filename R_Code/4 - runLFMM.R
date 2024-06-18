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

run_lfmm <- function(lfmm_file, environment_file, ...){
  
  geno_mat <- data.table::fread(lfmm_file, header = FALSE, verbose = FALSE, 
                                showProgress = FALSE, 
                                nThread = parallel::detectCores() / 2, 
                                logical01 = FALSE) %>%
    as.matrix
  
  env_mat <- read_delim(environment_file, delim = ' ', 
                        col_names = FALSE, show_col_types = FALSE) %>%
    as.matrix()
  
  lfmm2(geno_mat, env_mat, ...)
}

process_out <- function(x, NAMES){
  if(dim(x)[1] < dim(x)[2]){
    t(x) %>%
      set_colnames(NAMES) %>%
      as_tibble()
  } else {
    x %>%
      set_colnames(NAMES) %>%
      as_tibble()
  }
}


post_process_lfmm <- function(lfmm_model, lfmm_file, environment_file){
  
  geno_mat <- data.table::fread(lfmm_file, header = FALSE, verbose = FALSE, 
                                showProgress = FALSE, 
                                nThread = parallel::detectCores() / 2, 
                                logical01 = FALSE) %>%
    as.matrix
  
  env_mat <- read_delim(environment_file, delim = ' ', 
                        col_names = FALSE, show_col_types = FALSE) %>%
    as.matrix()
  
  effect_sizes <- lfmm_model@B %>%
    set_colnames('es') %>%
    as_tibble %>%
    mutate(locus_id = row_number())
  
  constant_loci <- apply(geno_mat, 2, function(a) all(a[1] == a))
  
  var_lfmm <- lfmm2.test(lfmm_model, 
                         geno_mat[,!constant_loci],
                         env_mat,
                         full = FALSE,
                         genomic.control = TRUE,
                         linear = TRUE)
  
  out_pz <- var_lfmm %>%
    extract(2:1) %>%
    map(process_out, NAMES = 'A') %>%
    map2(., names(.), function(x, y) rename_with(x, ~str_c(y, ., sep = '_'))) %>%
    bind_cols() %>%
    rename_with(~str_remove(., '_A')) %>%
    mutate(locus_id = which(!constant_loci))
  
  out <- full_join(effect_sizes, 
                   out_pz,
                   by = 'locus_id') %>%
    arrange(locus_id) %>%
    select(-locus_id)
  
  out
}


#### Read in Data ####
meta_data <- read_csv('../intermediate_files/preprocessed_metadata.csv', 
                      show_col_types = FALSE) %>%
  filter(!is.na(disease_resistance)) %>%
  mutate(resist_norm = qqnorm(disease_resistance, plot.it = FALSE)$x)

intial_data <- expand_grid(pop_group = c('acerv'),
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
  rowwise %>%
  mutate(impute_file = impute_missing(lfmm_file, structure_file, K)) %>%
  ungroup

#### Run LFMM ####
if(file.exists('../Results/lfmm_results.csv')){
  lfmm_results <- read_csv('../Results/lfmm_results.csv')
} else {
  lfmm_results <- imputed_data %>%
    select(pop_group, impute_file, phenotype_file, K, genomic_file) %>%
    rowwise(pop_group) %>%
    mutate(initial_lfmm = list(run_lfmm(lfmm_file = impute_file, environment_file = phenotype_file, 
                                        K = K, effect.sizes = TRUE)),
           genomic_data = list(read_rds(genomic_file))) %>%
    reframe(chromosome = genomic_data@chromosome,
              position = genomic_data@position,
              post_process_lfmm(initial_lfmm, 
                                impute_file,
                                phenotype_file)) %>%
    rename(z = zscores,
           p = pvalues) %>%
    mutate(fdr = p.adjust(p, method = 'fdr'),
           pop_group = str_replace_all(pop_group, c('^acerv$' = 'all', 'acerv_panama' = 'panama', 'acerv_florida' = 'florida')))
  write_csv(lfmm_results, '../Results/lfmm_results.csv')
}

#### Diagnostic Plots ####
ks_tests <- lfmm_results %>%
  group_by(pop_group) %>%
  sample_n(5000) %>%
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
  facet_grid(pop_group ~ .) +
  theme_classic() 
ggsave('../Results/lfmm_diagnostic_histograms.png', height = 8, width = 8)

lfmm_results %>%
  ggplot(aes(sample = p)) +
  stat_qq_line(distribution = qunif, dparams = list(min = 0, max = 1), size = 2, colour = 'red') +
  stat_qq(distribution = qunif, dparams = list(min = 0, max = 1), size = 0.01, alpha = 0.3) +
  facet_grid(pop_group ~ .) +
  theme_classic()
ggsave('../Results/lfmm_diagnostic_qqPlot.png', height = 8, width = 8)
