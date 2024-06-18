library(tidyverse)
library(magrittr)
library(adegenet)
library(hierfstat)
library(dartR)
library(furrr)
library(progressr)

NPERM <- 1000
plan('multisession', workers = 7L)

#### Functions ####
permute_fst <- function(gi, p){
  p()
  pop(gi) <- sample(gi@pop)
  wc(gi)$FST
}

find_loci_all_missing_pop <- function(gi){
  loci_remove <- seppop(gi) %>%
    map(~as.matrix(.x) %>% is.na %>% colSums() %>% divide_by(nInd(.x))) %>%
    map(~.x[.x == 1]) %>%
    unlist %>%
    names %>%
    str_remove('Ac\\.|Apa\\.|Apr\\.') %>%
    str_remove('\\.[AGCT]$') %>%
    unique
  
  loc_keep <- as.character(gi@loc.fac)[!as.character(gi@loc.fac) %in% loci_remove]
  gi[loc = loc_keep]
}

#### Read in Data ####
full_data <- read_rds('../intermediate_files/initial_full_genomic.rds')
unlinked_loci_acerv <- read_rds('../intermediate_files/unlinked_preprocessed_acerv.rds')
meta_data <- read_csv('../intermediate_files/preprocessed_metadata.csv', 
                      show_col_types = FALSE) %>%
  filter(!is.na(disease_resistance))

#### Subset to just our cervicornis & prolifera/palmata from baums ####
subset_data <- full_data[full_data@ind.names %in% unlinked_loci_acerv@ind.names | str_detect(full_data@ind.names, 'baum_a(palm|prol)'),
                         full_data@loc.names %in% unlinked_loci_acerv@loc.names]

#### Fst between Florida & Panama within Acerv #### 
pop(unlinked_loci_acerv) <- factor(as_tibble(strata(unlinked_loci_acerv))$location)
acerv_gi <- gl2gi(unlinked_loci_acerv)

(within_acerv_fst <- wc(acerv_gi))

with_progress({
  p <- progressor(steps = length(1:(NPERM - 1)))
  null_samples_acerv <- future_map_dbl(1:(NPERM - 1), ~permute_fst(gi = acerv_gi, p = p),
                                       .options = furrr_options(seed = TRUE))
})

#p-value
mean(c(null_samples_acerv, within_acerv_fst$FST) >= within_acerv_fst$FST)

#### Fst between Acerv and Apalm ####
acerv_apalm <- subset_data[str_detect(subset_data@ind.names, 'baum_aprol', negate = TRUE)]
pop(acerv_apalm) <- as_tibble(strata(acerv_apalm))$species

acerv_apalm_gi <- gl2gi(acerv_apalm)

(btw_species_fst <- wc(acerv_apalm_gi))

with_progress({
  p <- progressor(steps = length(1:(NPERM - 1)))
  null_samples_btw <- future_map_dbl(1:(NPERM - 1), ~permute_fst(gi = acerv_apalm_gi, p = p),
                                       .options = furrr_options(seed = TRUE))
})

mean(c(null_samples_btw, btw_species_fst$FST) >= btw_species_fst$FST)

#### Hybridization Analysis ####
pop(subset_data) <- as_tibble(strata(subset_data)) %>%
  mutate(species = if_else(species == 'Ac', str_c(species, str_sub(location, 1, 1)), 
                           as.character(species))) %>%
  pull(species) %>%
  factor

acropora_gi <- gl2gi(subset_data) %>%
  find_loci_all_missing_pop
acropora_split <-seppop(acropora_gi)


hybridization_classes <- tribble(
  ~'name', ~'P1', ~'P2', ~'gen', ~'P1_gen', ~'P2_gen', ~'label', ~'simple_label',
  'AcF', 'AcF', 'AcF', 'zero', 'zero', 'zero', 'AcF', 'Ac',
  'AcP', 'AcP', 'AcP', 'zero', 'zero', 'zero', 'AcP', 'Ac',
  'Apa', 'Apa', 'Apa', 'zero', 'zero', 'zero', 'Ap', 'Ap',
  
  'AcF_AcP_F1', 'AcF', 'AcP', 'first', 'zero', 'zero', 'AcF x AcP', 'Ac',
  'AcF_Apa_F1', 'AcF', 'Apa', 'first', 'zero', 'zero', 'AcF x Ap', 'F1',
  'AcP_Apa_F1', 'AcP', 'Apa', 'first', 'zero', 'zero', 'AcP x Ap', 'F1',
  
  'AcF_AcP_F2', 'AcF_AcP_F1', 'AcF_AcP_F1', 'second', 'first', 'first', '(AcF x AcP) x (AcF x AcP)', 'Ac',
  'AcF_Apa_F2', 'AcF_Apa_F1', 'AcF_Apa_F1', 'second', 'first', 'first', '(AcF x Ap) x (AcF x Ap)', 'F2',
  'AcP_Apa_F2', 'AcP_Apa_F1', 'AcP_Apa_F1', 'second', 'first', 'first', '(AcP x Ap) x (AcP x Ap)', 'F2',
  
  'AcF_AcP_F1_AcF_bx', 'AcF_AcP_F1', 'AcF', 'second', 'first', 'zero', '(AcF x AcP) x AcF', 'Ac',
  'AcF_AcP_F1_AcP_bx', 'AcF_AcP_F1', 'AcP', 'second', 'first', 'zero', '(AcF x AcP) x AcP', 'Ac',
  'AcF_AcP_F1_Apa_F1', 'AcF_AcP_F1', 'Apa', 'second', 'first', 'zero', '(AcF x AcP) x Ap', 'F1',
  
  'AcF_Apa_F1_AcF_bx', 'AcF_Apa_F1', 'AcF', 'second', 'first', 'zero', '(AcF x Ap) x AcF', 'Ac_bx',
  'AcF_Apa_F1_Apa_bx', 'AcF_Apa_F1', 'Apa', 'second', 'first', 'zero', '(AcF x Ap) x Ap', 'Ap_bx',
  'AcF_Apa_F1_AcP_x', 'AcF_Apa_F1', 'AcP', 'second', 'first', 'zero', '(AcF x Ap) x AcP', 'Ac_bx',
  
  'AcP_Apa_F1_AcP_bx', 'AcP_Apa_F1', 'AcP', 'second', 'first', 'zero', '(AcP x Ap) x AcP', 'Ac_bx',
  'AcP_Apa_F1_Apa_bx', 'AcP_Apa_F1', 'Apa', 'second', 'first', 'zero', '(AcP x Ap) x Ap', 'Ap_bx',
  'AcP_Apa_F1_AcF_x', 'AcP_Apa_F1', 'AcF', 'second', 'first', 'zero', '(AcP x Ap) x AcF', 'Ac_bx',
  
  'AcF_AcP_F1_x_AcF_Apa_F1', 'AcF_AcP_F1', 'AcF_Apa_F1', 'second', 'first', 'first', '(AcF x AcP) x (AcF x Ap)', 'Ac_bx',
  'AcF_AcP_F1_x_AcP_Apa_F1', 'AcF_AcP_F1', 'AcP_Apa_F1', 'second', 'first', 'first', '(AcF x AcP) x (AcP x Ap)', 'Ac_bx',
  'AcF_Apa_F1_x_AcP_Apa_F1', 'AcF_Apa_F1', 'AcP_Apa_F1', 'second', 'first', 'first', '(AcF x Ap) x (AcP x Ap)', 'F2'
)


if(file.exists('../intermediate_files/hybrid_pca.rds')){
  acropora_pca <- read_rds('../intermediate_files/hybrid_pca.rds')
} else { 
  hybrid_groups <- vector('list', nrow(hybridization_classes))
  names(hybrid_groups) <- hybridization_classes$name
  N_hybrid <- 50
  
  for(i in 1:length(hybrid_groups)){
    message('Start simulating hybrids for class ', hybridization_classes$name[i])
    message(Sys.time())
    #Get real data unless first/second gen hybrid is one of parents
    if(hybridization_classes$P1[i] %in% names(acropora_split)){
      p1 <- acropora_split[[hybridization_classes$P1[i]]]
    } else {
      p1 <- hybrid_groups[[hybridization_classes$P1[i]]]
    }
    
    if(hybridization_classes$P2[i] %in% names(acropora_split)){
      p2 <- acropora_split[[hybridization_classes$P2[i]]]
    } else {
      p2 <- hybrid_groups[[hybridization_classes$P2[i]]]
    }
    #Make Hybrids
    hybrid_groups[[i]] <- hybridize(p1, p2, n = N_hybrid, 
                                    pop = hybridization_classes$name[i], 
                                    hyb.label = str_c(hybridization_classes$name[i], 
                                                      1:N_hybrid, sep = '_'))
    message('Finished simulating hybrids for class ', hybridization_classes$name[i])
    message(Sys.time())
  }
  
  #### Plot PCA with observed & simulated hybrids ####
  acropora_hybrids <- repool(c(acropora_gi, hybrid_groups)) 
  
  acropora_pca <- acropora_hybrids %>%
    scaleGen(., center = FALSE, scale = FALSE, NA.method = "mean") %>%
    dudi.pca(center = TRUE, scale = TRUE, scannf = FALSE, nf = nLoc(acropora_gi))
  write_rds(acropora_pca, '../intermediate_files/hybrid_pca.rds')
  
}

pct_explained <- acropora_pca$eig / sum(acropora_pca$eig)


wesanderson::wes_palette("Zissou1", 5, type = "continuous")
location_colour <- set_names(wesanderson::wes_palette("Zissou1", 9, type = "continuous")[c(2, 8)],
                             c('AcF', 'AcP'))

colour_scale <- wesanderson::wes_palette("Zissou1", 5, type = "continuous")[c(1,2,3,3,4,5)] %>%
  set_names(c('Ac', 'Ac_bx', 'F1', 'F2', 'Ap_bx', 'Ap'))

plot_data <- acropora_pca$li %>%
  as_tibble(rownames = 'ID') %>%
  select(ID, Axis1, Axis2) %>%
  mutate(is_simulated = if_else(row_number() %in% 1:nInd(acropora_gi), 'Observed', 'Simulated')) %>%
  
  mutate(hybrid_type = if_else(is_simulated == 'Simulated',
                               str_extract(ID, str_c(rev(unique(hybridization_classes$name)), 
                                             collapse = '|')),
                               str_extract(ID, 'Apa|apalm|aprol|F1|F2|AcBx|ApBx|Ac'))) %>%
  left_join(select(hybridization_classes, name, label, simple_label),
            by = c('hybrid_type' = 'name')) %>%
  mutate(label = if_else(is_simulated == 'Observed',
                         case_when(hybrid_type == 'Ac' ~ if_else(str_detect(ID, 'FL'), 'AcF', 'AcP'),
                                   hybrid_type == 'apalm' ~ 'Ap',
                                   hybrid_type == 'aprol' ~ '(Ac x Ap)'),
                         label)) %>%
  
  mutate(simple_label = if_else(is_simulated == 'Observed',
                                case_when(hybrid_type == 'Ac' ~ if_else(str_detect(ID, 'FL'), 'Ac', 'Ac'),
                                          hybrid_type == 'apalm' ~ 'Ap',
                                          hybrid_type == 'aprol' ~ 'F1'),
                                simple_label)) %>%
  mutate(location = case_when(str_detect(ID, 'FL') & is_simulated == 'Observed' ~ 'FL',
                              str_detect(ID, 'PA') & is_simulated == 'Observed' ~ 'PA',
                              is_simulated == 'Observed' & hybrid_type == 'apalm' ~ 'palmata',
                              is_simulated == 'Observed' & hybrid_type == 'aprol' ~ 'prolifera',
                              simple_label == 'F2' ~ 'F2',
                              TRUE ~ 'other')) %>%
  mutate(simple_label = factor(simple_label, levels = c('Ac', 'Ac_bx', 'F1', 'F2', 'Ap_bx', 'Ap')),
         location = factor(location, levels = c('FL', 'PA', 'prolifera', 'palmata', 'other', 'F2')))

library(ggtext)
filter(plot_data, is_simulated == 'Simulated') %>%
  ggplot(aes(x = Axis1, y = Axis2, fill = simple_label)) +
  geom_hex(bins = 100, alpha = 0.8) +
  geom_point(data = filter(plot_data, is_simulated != 'Simulated'), size = 3,
             colour = 'black', aes(shape = location)) +
  scale_fill_manual(values = colour_scale) +
  scale_shape_manual(values = c('FL' = 21, 'PA' = 24, 'palmata' = 22, 'prolifera' = 23),
                     labels = c('*A. cervicornis* FL', '*A. cervicornis* PA', 
                                '*A. prolifera* (hybrid)', '*A. palmata*')) +
  labs(x = str_c('PC1 (', scales::percent(pct_explained[1]), ')'),
       y = str_c('PC1 (', scales::percent(pct_explained[2]), ')'),
       fill = 'Type',
       shape = 'Species') +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5))) +
  theme(aspect.ratio = 0.5,
        axis.text = element_text(colour = 'black', size = 24),
        axis.title = element_text(colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 24),
        legend.title = element_text(colour = 'black', size = 18),
        # legend.text = element_text(face = 'bold', colour = 'black', size = 18),
        legend.text = element_markdown(colour = 'black', size = 18))
ggsave('../Results/hybridization_pca.png', height = 10, width = 15)




plot_data %>%
  
  ggplot(aes(x = Axis1, y = Axis2, colour = simple_label, shape = location)) +
  geom_point(size = 3) +
  # geom_text(aes(label = simple_label)) +
  facet_wrap(~ is_simulated, ncol = 1) +
  scale_shape_manual(values = c('FL' = 'circle', 'PA' = 'triangle', 'palmata' = 'square', 'prolifera' = 'diamond',
                                'other' = 'bullet', 'F2' = 'cross')) +
  scale_colour_manual(values = colour_scale) +
  labs(x = str_c('PC1 (', scales::percent(pct_explained[1]), ')'),
       y = str_c('PC1 (', scales::percent(pct_explained[2]), ')')) +
  theme_classic() +
  labs(colour = 'Species') +
  guides(colour = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5))) +
  theme(aspect.ratio = 0.5,
        axis.text = element_text(face = 'bold', colour = 'black', size = 24),
        axis.title = element_text(face = 'bold', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold', colour = 'black', size = 24),
        legend.title = element_text(face = 'bold', colour = 'black', size = 18),
        legend.text = element_text(face = 'bold', colour = 'black', size = 18))

acropora_pca$li %>%
  as_tibble(rownames = 'ID') %>%
  select(ID, Axis1, Axis2) %>%
  mutate(is_simulated = if_else(row_number() %in% 1:nInd(acropora_gi), 'Observed', 'Simulated')) %>%
  
  mutate(hybrid_type = if_else(is_simulated == 'Simulated',
                               str_extract(ID, str_c(rev(unique(hybridization_classes$name)), 
                                                     collapse = '|')),
                               str_extract(ID, 'Apa|apalm|aprol|F1|F2|AcBx|ApBx|Ac'))) %>%
  left_join(select(hybridization_classes, name, label),
            by = c('hybrid_type' = 'name')) %>%
  mutate(label = if_else(is_simulated == 'Observed',
                         case_when(hybrid_type == 'Ac' ~ if_else(str_detect(ID, 'FL'), 'AcF', 'AcP'),
                                   hybrid_type == 'apalm' ~ 'Ap',
                                   hybrid_type == 'aprol' ~ '(Ac x Ap)'),
                         label)) %>%
  
  ggplot(aes(x = Axis1, y = Axis2, colour = label)) +
  geom_point(size = 3) +
  facet_wrap(label ~ is_simulated) +
  labs(x = str_c('PC1 (', scales::percent(pct_explained[1]), ')'),
       y = str_c('PC1 (', scales::percent(pct_explained[2]), ')')) +
  theme_classic() +
  labs(colour = 'Species') +
  guides(colour = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5))) +
  theme(aspect.ratio = 0.5,
        # axis.text = element_text(face = 'bold', colour = 'black', size = 24),
        # axis.title = element_text(face = 'bold', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA),
        strip.background = element_blank(),
        # strip.text = element_text(face = 'bold', colour = 'black', size = 24),
        legend.title = element_text(face = 'bold', colour = 'black', size = 18),
        legend.text = element_text(face = 'bold', colour = 'black', size = 18),
        legend.position = 'none')


#### PC-component vs Disease Resistance ####
resistance_v_hybrid <- acropora_pca$li %>%
  as_tibble(rownames = 'ID') %>%
  select(ID, Axis1) %>%
  inner_join(meta_data, by = 'ID') 

with(resistance_v_hybrid, cor.test(disease_resistance, Axis1, method = 'kendall'))

resist_hybrid_lm <- lm(disease_resistance ~ Axis1 * location, data = resistance_v_hybrid) 
summary(resist_hybrid_lm)
car::Anova(resist_hybrid_lm, type = 2) 

performance::check_model(resist_hybrid_lm)

ggplot(resistance_v_hybrid, aes(x = Axis1, y = disease_resistance, colour = location)) +
  geom_point()
