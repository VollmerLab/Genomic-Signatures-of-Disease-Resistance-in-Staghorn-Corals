#### Libraries ####
library(tidyverse)
library(magrittr)
library(LEA)
library(mgcv)

meta_data <- read_csv('../intermediate_files/preprocessed_metadata.csv', 
                      show_col_types = FALSE) %>%
  filter(!is.na(disease_resistance))

lfmm_results <- read_csv('../Results/lfmm_results.csv', show_col_types = FALSE) %>%
  filter(pop_group == 'all',
         model_family == 'norm') %>%
  select(-pop_group, -model_family) %>%
  mutate(locus = str_c(chromosome, position, sep = '-'), .after = 'position') %>%
  select(-es) %>% 
  left_join(read_csv('../intermediate_files/lfmm_snp_clumps.csv', 
                     show_col_types = FALSE),
            by = c('chromosome', 'locus'))

updated_dr <- read_csv('../intermediate_files/disease_resistance.csv', show_col_types = FALSE) %>%
  filter(treatment == 'D') %>%
  select(gen_id, disease_resistance_se) %>%
  inner_join(select(meta_data, ID, gen_id, disease_resistance),
             by = c('gen_id')) %>%
  arrange(match(ID, meta_data$ID))

lfmm_linkage_loci <- lfmm_results %>%
  mutate(clump_id = str_c(chromosome, clump_id, sep = '-')) %>%
  filter(fdr < 0.05) %>%
  group_by(clump_id) %>%
  filter(p == min(p),
         is_index) %>%
  # sample_n(1) %>%
  ungroup 

calculate_pgs <- function(locus_effects, genotypes){
  # Genotypes: Ind x Locus
  # Effect Size: Locus x 1 
  # Product: Ind x 1 (sum across loci)
  
  pgs_numerator <- genotypes[,locus_effects$locus] %*% matrix(locus_effects$es)
  
  (pgs_numerator / (2 * nrow(locus_effects))) %>%
    set_colnames('pgs') %>%
    set_rownames(rownames(genotypes)) %>%
    as_tibble(rownames = 'ID')
}

#### Bootstrap LFMM - PGS ####
# https://patrick-rockenschaub.com/post/2022-01-26-tidymodels-optimism-bootstrap/
# https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-020-01201-w


if(file.exists('../intermediate_files/optimism_values.csv')){
  optimisim_values <- read_csv('../intermediate_files/optimism_values.csv', show_col_types = FALSE)
} else {
  boot_resample <- function(lfmm_file, phenotype_file, sample_size, loci_names, observed){
    
    if(!observed){
      boot_rows <- sample(sample_size, sample_size, replace = TRUE)
    } else {
      boot_rows <- 1:sample_size
    }
    
    
    out <- tibble(split = if_else(observed, 'obs', 'boot'),
                  which_rows = list(boot_rows)) %>%
      rowwise %>%
      mutate(geno_mat = list(data.table::fread(lfmm_file, header = FALSE, verbose = FALSE, 
                                               showProgress = FALSE, 
                                               nThread = parallel::detectCores() / 2, 
                                               logical01 = FALSE) %>%
                               as.matrix %>%
                               extract(which_rows, ) %>%
                               set_colnames(loci_names)),
             
             pheno_mat = list(read_delim(phenotype_file, delim = ' ', 
                                         col_names = FALSE, show_col_types = FALSE) %>%
                                as.matrix() %>%
                                extract(which_rows, )))
    out
  }
  
  observed <- boot_resample("../lea_projects/preprocessed_acerv.lfmm_imputed.lfmm", 
                            "../lea_projects/preprocessed_acerv.env",
                            sample_size = nrow(meta_data),
                            loci_names = lfmm_results$locus,
                            observed = TRUE) %>%
    ungroup
  
  #Fit observed model
  observed_pgs_fit <- observed %>%
    mutate(K = 2) %>%
    rowwise %>%
    mutate(lfmm = list(lfmm2(geno_mat, pheno_mat, K = K, effect.sizes = TRUE)),
           loci = list(lfmm@B %>%
                         set_colnames(str_c('es')) %>%
                         as_tibble %>%
                         bind_cols(lfmm_results))) %>%
    select(-lfmm) %>%
    expand_grid(threshold_p = c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 
                                5e-3, 1e-2, 0.05, 0.1, 0.2, 0.3, NA_real_)) %>%
    rowwise %>%
    mutate(obs_loci = case_when(!is.na(threshold_p) ~ list(filter(loci, is_index, p < threshold_p)),
                                TRUE ~ list(filter(loci, locus %in% lfmm_linkage_loci$locus))),
           n_loci = nrow(obs_loci)) %>%
    
    mutate(obs_pgs = list(calculate_pgs(obs_loci, geno_mat) %>%
                            mutate(pheno = pheno_mat,
                                   pheno_var = updated_dr$disease_resistance_se[which_rows] ^ 2)),
           # obs_model = list(lm(pheno ~ pgs, data = obs_pgs)),
           obs_model = list(gam(pheno ~ pgs, 
                                family=betar(link="logit"), 
                                data = obs_pgs, 
                                weights = pheno_var)),
           obs_pgs = list(bind_cols(obs_pgs, 
                                    predict(obs_model, 
                                            newdata = obs_pgs, 
                                            se.fit = TRUE, 
                                            type = 'response')))) %>%
    mutate(yardstick::rsq_trad(obs_pgs, truth = pheno, estimate = fit)) %>%
    rename(obs_r2 = .estimate) %>%
    ungroup
  
  #Fit bootstrap models
  bootstrap_resamples <- expand_grid(boot_resamp = 1:100,
                                     K = 2) %>%
    rowwise(boot_resamp, K) %>%
    reframe(boot_resample("../lea_projects/preprocessed_acerv.lfmm_imputed.lfmm", 
                          "../lea_projects/preprocessed_acerv.env",
                          sample_size = nrow(meta_data),
                          loci_names = lfmm_results$locus,
                          observed = FALSE)) %>%
    rowwise %>%
    mutate(lfmm = list(lfmm2(geno_mat, pheno_mat, K = K, effect.sizes = TRUE)),
           loci = list(lfmm@B %>%
                         set_colnames(str_c('es')) %>%
                         as_tibble %>%
                         bind_cols(lfmm_results))) %>%
    select(-lfmm) %>%
    expand_grid(threshold_p = c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 
                                5e-3, 1e-2, 0.05, 0.1, 0.2, 0.3, NA_real_)) %>%
    rowwise %>%
    # mutate(boot_loci = list(filter(loci, is_index, p < threshold_p))) %>%
    mutate(boot_loci = case_when(!is.na(threshold_p) ~ list(filter(loci, is_index, p < threshold_p)),
                                TRUE ~ list(filter(loci, locus %in% lfmm_linkage_loci$locus))),
           n_loci = nrow(boot_loci)) %>%
    
    mutate(boot_pgs = list(calculate_pgs(boot_loci, geno_mat) %>%
                             mutate(pheno = pheno_mat,
                                    pheno_var = updated_dr$disease_resistance_se[which_rows] ^ 2)),
           # boot_model = list(lm(pheno ~ pgs, data = boot_pgs)),
           boot_model = list(gam(pheno ~ pgs, 
                                 family=betar(link="logit"), 
                                 data = boot_pgs, 
                                 weights = pheno_var)),
           
           boot_pgs = list(bind_cols(boot_pgs, predict(boot_model, newdata = boot_pgs, se.fit = TRUE, 
                                                       type = 'response')))) %>%
    mutate(yardstick::rsq_trad(boot_pgs, truth = pheno, estimate = fit)) %>%
    rename(boot_r2 = .estimate) %>%
    ungroup %>%
    select(threshold_p, n_loci, starts_with('boot'))
  
  
  #Calculate optimism
  optimisim_values <- bind_cols(bootstrap_resamples,
                                observed) %>%
    rowwise %>%
    mutate(obs_pgs = list(calculate_pgs(boot_loci, geno_mat) %>%
                            mutate(pheno = pheno_mat,
                                   pheno_var = updated_dr$disease_resistance_se[which_rows] ^ 2)),
           obs_pgs = list(bind_cols(obs_pgs, predict(boot_model, newdata = obs_pgs, se.fit = TRUE, 
                                                     type = 'response')))) %>%
    mutate(yardstick::rsq_trad(obs_pgs, truth = pheno, estimate = fit)) %>%
    rename(obs_r2 = .estimate) %>%
    mutate(optimism = boot_r2 - obs_r2) %>%
    select(threshold_p, n_loci, ends_with('r2'), optimism) %>%
    group_by(threshold_p, n_loci) %>%
    summarise(across(where(is.numeric), list(mean = mean, sd = sd)),
              n_boot = n(), .groups = 'drop') %>%
    left_join(select(observed_pgs_fit, threshold_p, n_loci, obs_r2),
              by = c('threshold_p', 'n_loci')) %>%
    mutate(adj_obs_r2 = obs_r2 - optimism_mean,
           adj_obs_lwr_95 = adj_obs_r2 + c(-2) * optimism_sd / sqrt(n_boot),
           adj_obs_upr_95 = adj_obs_r2 + c(2) * optimism_sd / sqrt(n_boot)) 
  write_csv(optimisim_values, '../intermediate_files/optimism_values.csv')
}

optimisim_values %>%
  select(n_loci, adj_obs_r2, optimism_mean, optimism_sd)

filter(optimisim_values, n_loci %in% c(10, 63)) %>% 
  arrange(n_loci) %>%
  select(n_loci, obs_r2, starts_with('adj_obs'), optimism_sd)

thresh_p_plot <- optimisim_values %>%
  mutate(adj_obs_r2 = obs_r2 - optimism_mean,
         adj_obs_lwr_95 = adj_obs_r2 - optimism_sd,
         adj_obs_upr_95 = adj_obs_r2 + optimism_sd) %>%
  
  ggplot(aes(x = log(threshold_p, base = 10), y = adj_obs_r2, 
             ymin = adj_obs_lwr_95, ymax = adj_obs_upr_95)) +
  # geom_path(aes(y = obs_r2), colour = 'black') +
  # geom_pointrange(aes(y = obs_r2), colour = 'black') +
  
  geom_path(colour = 'black') +
  geom_pointrange(colour = 'black') +
  
  # geom_text_repel(aes(label = scales::pvalue(threshold_p, accuracy = 1e-6))) +
  scale_x_continuous(labels = scales::math_format()) +
  annotation_logticks(sides = 'b') +
  labs(y = expression(r^2),
       x = 'P-Value Threshold',
       colour = NULL,
       fill = NULL) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 16),
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 16))
ggsave('../Results/optimism_adjusted_r2.png', height = 5, width = 5)

n_loci_plot <- optimisim_values %>%
  arrange(n_loci) %>%
  
  ggplot(aes(x = log(n_loci, base = 10), y = adj_obs_r2, 
             ymin = adj_obs_lwr_95, ymax = adj_obs_upr_95)) +
  # geom_path(aes(y = obs_r2), colour = 'black') +
  # geom_pointrange(aes(y = obs_r2), colour = 'black') +
  
  geom_path(colour = 'black') +
  geom_pointrange(colour = 'black') +
  
  # geom_text_repel(aes(label = scales::pvalue(threshold_p, accuracy = 1e-6))) +
  scale_x_continuous(labels = scales::math_format()) +
  annotation_logticks(sides = 'b') +
  labs(y = expression(r^2),
       x = 'Linkage Groups',
       colour = NULL,
       fill = NULL) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 16),
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 16))
