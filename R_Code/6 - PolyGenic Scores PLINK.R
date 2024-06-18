#Clump & threshold snps

#### Libraries ####
library(tidyverse)
library(magrittr)
library(adegenet)
library(broom)
library(patchwork)
library(ggdist)
library(gghalves)
library(ggrepel)
library(ggprism)
library(ggside)
library(mgcv)
library(RColorBrewer)

alpha <- 0.05

#### Functions ####
replace_na <- function(mat){
  mat[is.na(mat)] <- -1
  mat
}

impute_mean <- function(x){
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}

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

#### Data ####
full_metadata <- read_csv('../intermediate_files/clone_metadata.csv', show_col_types = FALSE) %>%
  group_by(clone_group) %>%
  filter(pct_missing == min(pct_missing)) %>% 
  ungroup %>%
  filter((species == 'Ac' & data_origin == 'vollmer') | species == 'Apa' | species == 'Apr') %>%
  filter((species == 'Ac' & pct_missing  < 0.3) | species == 'Apa' | species == 'Apr') %>%
  mutate(species_location = if_else(data_origin == 'vollmer', location, species)) 


lfmm_results <- read_csv('../Results/lfmm_results.csv', show_col_types = FALSE) %>%
  filter(pop_group == 'all') %>%
  select(-pop_group) %>%
  mutate(locus = str_c(chromosome, position, sep = '-'), .after = 'position') 

full_snp_info <- read_csv('../intermediate_files/snps_with_mutationType.csv', show_col_types = FALSE)

unfiltered_snps <- read_rds('../intermediate_files/initial_full_genomic.rds')

functional_annotations <- read_csv('../genome/genomic_structural_annotation.csv', show_col_types = FALSE) %>%
  filter(structure_type == 'CDS') %>%
  select(gene_id, structure_attributes) %>%
  distinct

#### Calculate linkage disequilibrium within each contig ####
if(!file.exists('../variant_calling/linkage/lfmm_snp_clumps.csv.gz')){
  dir.create('../variant_calling/linkage')
  
  #Make Position File
  subset_loci <- unfiltered_snps[unfiltered_snps@ind.names %in% full_metadata$ID[full_metadata$data_origin == 'vollmer'], 
                                 unfiltered_snps@loc.names %in% lfmm_results$locus]
  posFile <- tibble(locus = subset_loci@loc.names) %>%
    separate(locus, into = c('contig', 'position'), sep = '-') %>%
    nest_by(contig) %>%
    mutate(data = list(mutate(data, contig = contig) %>%
                         select(contig, position))) %>%
    mutate(file_out = str_c('../variant_calling/linkage/', contig, '_posFile.pos')) 
  
   walk2(posFile$data, posFile$file_out, ~write_tsv(.x, .y, col_names = FALSE))

  #Make Genotype File
  genoFile <- as.matrix(subset_loci) %>%
    t %>%
    replace_na %>%
    as_tibble(rownames = "locus") %>%
    left_join(select(lfmm_results, chromosome, locus),
              by = 'locus') %>%
    nest_by(chromosome) %>%
    mutate(data = list(select(data, -locus)),
           file_out = str_c('../variant_calling/linkage/', chromosome, '_genotypes.tsv.gz')) %>%
    ungroup
  
  walk2(genoFile$data, genoFile$file_out, ~write_tsv(.x, .y, col_names = FALSE))
  
  # Run below block on SLURM
  # sbatch \
  #   --output=$(pwd)/variant_calling/linkage/linkage_%A_%a.output \
  #   --array=0-$((355-1))%30 \
  #   $(pwd)/bash_code/runNGSLD_array.slurm \
  #     $(pwd)/variant_calling/linkage
  
}

#### Clump SNPs into linkage blocks ####
if(!file.exists('../variant_calling/linkage/lfmm_snp_clumps.csv.gz')){
  # Run below block on SLURM
  # sbatch \
  #   --output=$(pwd)/variant_calling/linkage/linkage_clumping_%j.output \
  #   $(pwd)/bash_code/runRscript.slurm \
  #     $(pwd)/bash_code/r_utils/clumpify_lfmm_snps.R \
  #       $(pwd)/variant_calling/linkage \
  #       $(pwd)/Results/lfmm_results.csv
} else {
  snp_clumps <- read_csv('../variant_calling/linkage/lfmm_snp_clumps.csv.gz', 
                         show_col_types = FALSE)
}

#### Impute Missing Genotypes ####
unfiltered_snps <- unfiltered_snps[unfiltered_snps@ind.names %in% full_metadata$ID,
                                   unfiltered_snps@loc.names %in% lfmm_results$locus] 
setPop(unfiltered_snps) <- ~species

imputed_genotypes <- seppop(unfiltered_snps) %>%
  map(tab, NA.method = 'mean') %>%
  do.call(rbind, .) %>%
  
  #if all missing in pop need to impute from everyone
  apply(2, impute_mean)

#### Threshold Linkage PGS Linkage Blocks ####
threshold_pgs <- expand_grid(threshold_p = c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 
                                             5e-3, 1e-2, 0.05, 0.1, 0.2, 0.3)) %>%
  rowwise(threshold_p) %>%
  
  mutate(data = list(full_join(lfmm_results, snp_clumps,
                               by = c('chromosome', 'locus')) %>%
                       filter(is_index) %>%
                       mutate(clump_id = str_c(chromosome, clump_id, sep = '-')))) %>%
  mutate(data = list(filter(data, p < threshold_p))) %>%
  reframe(n_clumps = nrow(data),
            chromosomes = list(data$clump_id), 
            calculate_pgs(data, imputed_genotypes)) %>%
  left_join(full_metadata, 
            by = 'ID') %>%
  nest(pgs_data = c(ID:species_location)) %>%
  rowwise %>%
  # mutate(glance(lm(disease_resistance ~ pgs, data = pgs_data))) %>%
  mutate(model = list(gam(disease_resistance ~ pgs, family=betar(link="logit"), data = pgs_data, 
                  weights = (disease_resistance_se^2))),
         r.squared = summary(model)$r.sq,
         pred_data = list(tibble(pgs = modelr::seq_range(pgs_data$pgs, n = 100)) %>%
                            bind_cols(., predict(model, newdata = ., se.fit = TRUE, 
                                                 type = 'response')))) %>%
  ungroup

threshold_pgs %>%
  select(threshold_p, n_clumps, chromosomes, r.squared) %>%
  unnest(chromosomes) %>%
  separate(chromosomes, into = c('chromosome', 'clump_id'), sep = '-') %>%
  write_csv('../intermediate_files/clumping_thersholding_results.csv')

n_clump_plot <- threshold_pgs %>%
  ggplot(aes(x = log(n_clumps, base = 10), y = r.squared)) +
  geom_path() +
  geom_point() +
  # geom_text_repel(aes(label = scales::comma(n_clumps))) +
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

thres_p_plot <- threshold_pgs %>%
  ggplot(aes(x = log(threshold_p, base = 10), y = r.squared)) +
  geom_path() +
  geom_point() +
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

n_clump_plot + thres_p_plot & plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 24))
ggsave('../Results/pgs_thresholding_r2.png', height = 6, width = 10)

#### Pick PGS threshold ####
p_cutoff <- 0.0001

polygenic_scores <- threshold_pgs %>%
  filter(threshold_p == p_cutoff) %>%
  select(pgs_data) %>%
  unnest(pgs_data) 

lm(disease_resistance ~ pgs, polygenic_scores) %>% car::Anova(type = '2')
lm(disease_resistance ~ pgs, polygenic_scores) %>% summary

#### PGS of linkage blocks but only those which are significant in LFMM ####
lfmm_linkage_loci <- full_join(lfmm_results, snp_clumps,
          by = c('chromosome', 'locus')) %>%
  mutate(clump_id = str_c(chromosome, clump_id, sep = '-')) %>%
  filter(fdr < alpha) %>%
  group_by(clump_id) %>%
  filter(p == min(p),
         is_index) %>%
  # sample_n(1) %>%
  ungroup 

lfmm_linkage_loci_pgs <- lfmm_linkage_loci %>%
  calculate_pgs(imputed_genotypes) %>%
  left_join(full_metadata, 
            by = 'ID') 

ten_snp_model <- gam(disease_resistance ~ pgs, family=betar(link="logit"), 
    data = filter(lfmm_linkage_loci_pgs, !is.na(disease_resistance)), 
    weights = (disease_resistance_se^2)) 

summary(ten_snp_model); anova(ten_snp_model)

#### Make Combined Figure ####
plot_data <- bind_rows(all = polygenic_scores %>%
                         filter(!is.na(disease_resistance)) %>%
                         mutate(colour_group = case_when(disease_resistance >= quantile(disease_resistance, 0.75) &
                                                           pgs >= quantile(pgs, 0.75) ~ 'A',
                                                         disease_resistance > quantile(disease_resistance, 0.5) &
                                                           pgs > quantile(pgs, 0.5) ~ 'B',
                                                         disease_resistance <= quantile(disease_resistance, 0.25) &
                                                           pgs <= quantile(pgs, 0.25) ~ 'D',
                                                         disease_resistance < quantile(disease_resistance, 0.5) &
                                                           pgs < quantile(pgs, 0.5) ~ 'C',
                                                         
                                                         TRUE ~ 'error')),
                       lfmm_sig = lfmm_linkage_loci_pgs %>%
                         filter(!is.na(disease_resistance)) %>%
                         mutate(colour_group = case_when(disease_resistance >= quantile(disease_resistance, 0.75) &
                                                           pgs >= quantile(pgs, 0.75) ~ 'A',
                                                         disease_resistance > quantile(disease_resistance, 0.5) &
                                                           pgs > quantile(pgs, 0.5) ~ 'B',
                                                         disease_resistance <= quantile(disease_resistance, 0.25) &
                                                           pgs <= quantile(pgs, 0.25) ~ 'D',
                                                         disease_resistance < quantile(disease_resistance, 0.5) &
                                                           pgs < quantile(pgs, 0.5) ~ 'C',
                                                         
                                                         TRUE ~ 'error')),
                       .id = 'type') 

models_preds <- plot_data %>%
  group_by(type) %>%
  summarise(the_model = list(gam(disease_resistance ~ pgs, family=betar(link="logit"),  
                                 weights = (disease_resistance_se^2))),
            pred_data = list(tibble(pgs = modelr::seq_range(pgs, n = 100))),
            .groups = 'rowwise') %>%
  mutate(pred_data = list(bind_cols(pred_data, predict(the_model, newdata = pred_data, se.fit = TRUE, 
                                          type = 'response'))),
         r.squared = summary(the_model)$r.sq) %>%
  ungroup %>%
  mutate(n_clump = if_else(type == 'all', threshold_pgs %>%
                             filter(threshold_p == p_cutoff) %>%
                             pull(n_clumps), 
                           nrow(lfmm_linkage_loci))) %>%
  mutate(type = factor(type, levels = c('lfmm_sig', 'all'))) %>%
  mutate(annotation = c("B", 'A'))

facet_labels <- models_preds %$%
  set_names(annotation, type)

models_preds$the_model %>%
  map(summary)


colour_options <- set_names(rev(c(brewer.pal(3, 'Reds')[c(3,2)], 
                              brewer.pal(3, 'Blues')[c(2,3)], 'gray50')),
                            c('error', LETTERS[1:4]))

wesanderson::wes_palette("Zissou1", 9, type = "continuous")
colour_options <- set_names(c(wesanderson::wes_palette("Zissou1", 9, type = "continuous")[c(1,3)],
                              'gray50', 
                              wesanderson::wes_palette("Zissou1", 9, type = "continuous")[c(5,9)]),
          c('A', 'B', 'error', 'C', 'D'))

plot_data %>%
  group_by(ID) %>%
  # mutate(colour_group = colour_group[type == 'lfmm_sig']) %>%
  ungroup %>%
  mutate(type = factor(type, levels = c('lfmm_sig', 'all'))) %>%
  mutate(colour_group = case_when(colour_group == 'A' ~ '(0,1]',
                                  colour_group == 'B' ~ '(1,2]',
                                  colour_group == 'error' ~ '(2,3]',
                                  colour_group == 'C' ~ '(3,4]',
                                  colour_group == 'D' ~ '(4,5]')) %>%
  
  ggplot(aes(x = pgs, y = disease_resistance)) +
  geom_vline(data = . %>% group_by(type) %>%
               reframe(pgs = quantile(pgs, c(0.25, 0.75))),
             aes(xintercept = pgs), linetype = 'solid', linewidth = 0.1,
             colour = 'gray50') +

  geom_vline(data = . %>% group_by(type) %>% 
               reframe(pgs = quantile(pgs, c(0.5))),
             aes(xintercept = pgs), linetype = 'solid', linewidth = 0.25,
             colour = 'black') +
  
  geom_hline(data = . %>% group_by(type) %>%
               reframe(disease_resistance = quantile(disease_resistance,
                                                       c(0.25, 0.75))),
             aes(yintercept = disease_resistance), linetype = 'solid', linewidth = 0.1,
             colour = 'gray50') +
  
  geom_hline(data = . %>% group_by(type) %>% 
               reframe(disease_resistance = quantile(disease_resistance,
                                                       c(0.5))),
             aes(yintercept = disease_resistance), linetype = 'solid', linewidth = 0.25,
             colour = 'black') +
  
  # geom_smooth(method = 'lm', formula = y ~ x, colour = 'black') +
  geom_ribbon(data = select(models_preds, type, pred_data) %>%
                unnest(pred_data), #inherit.aes = FALSE,
              aes(x = pgs, ymin = fit - se.fit, ymax = fit + se.fit, y = NULL),
              alpha = 0.5) +
  geom_line(data = select(models_preds, type, pred_data) %>%
              unnest(pred_data), #inherit.aes = FALSE,
            aes(x = pgs, y = fit)) +
  geom_point(aes(shape = location, colour = colour_group), size = 2) +
  geom_xsideboxplot(aes(group = 1),
                    show.legend = FALSE,
                    orientation = "y") +
  geom_ysideboxplot(aes(group = 1),
                    show.legend = FALSE,
                    orientation = "x") +
  scale_colour_manual(values = unname(colour_options),
                      labels = c('Highly\nResistant', '', '', 'Highly\nSusceptible')) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  guides(colour = guide_colorsteps(reverse = TRUE),
         shape = 'none') +
  facet_wrap(~ type, labeller = as_labeller(facet_labels)) +
  ggside(collapse = "y") +
  labs(y = 'Disease Resistance',
       x = 'Polygenic Score',
       colour = NULL,
       fill = NULL,
       shape = NULL) +
  # coord_equal(ratio = 4) +
  theme_classic() +
  theme(aspect.ratio=1,
        panel.background = element_rect(colour = 'black', linewidth = 1),
        axis.text = element_text(colour = 'black', size = 10),
        axis.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 10),
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 16, hjust = 0),
        ggside.panel.border = element_blank(),
        ggside.panel.background = element_blank(),
        ggside.axis.text = element_blank(),
        ggside.axis.line = element_blank(),
        ggside.axis.ticks = element_blank(),
        ggside.panel.scale = 0.01)
ggsave('../Results/pgs_v_diseaseResistance_combined.png', height = 7, width = 10)



plot_data %>%
  mutate(group = case_when(colour_group == 'A' ~ 'winner',
                           colour_group == 'D' ~ 'loser',
                           TRUE ~ 'middle')) %>%
  count(type, species_location, group) %>%
  pivot_wider(names_from = group, values_from = n)

plot_data %>%
  group_by(type) %>%
  mutate(quadrant = case_when(disease_resistance >= quantile(disease_resistance, 0.5) &
                                pgs >= quantile(pgs, 0.5) ~ 'topright',
                              disease_resistance >= quantile(disease_resistance, 0.5) &
                                pgs < quantile(pgs, 0.5) ~ 'topleft',
                              disease_resistance < quantile(disease_resistance, 0.5) &
                                pgs >= quantile(pgs, 0.5) ~ 'bottomright',
                              disease_resistance < quantile(disease_resistance, 0.5) &
                                pgs < quantile(pgs, 0.5) ~ 'bottomleft')) %>%
  # ggplot(aes(x = pgs, y = disease_resistance, colour = quadrant)) +
  # geom_point() +
  # facet_wrap(~type)
  
  count(type, quadrant, species_location) %>%
  pivot_wider(names_from = species_location,
              values_from = n, values_fill = 0L)

