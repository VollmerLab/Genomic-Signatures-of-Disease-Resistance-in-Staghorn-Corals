#### Libraries ####
library(tidyverse)
library(magrittr)
library(adegenet)
library(LEA)
library(patchwork)

#### Functions ####
write_lea_genetics <- function(genomic_data, outname, model_family){
  
  base_name <- str_remove(outname, str_c(dirname(outname), '/')) %>%
    str_remove('\\.rds$') %>%
    str_c('../lea_projects/', ., if_else(model_family == 'binomial', '_binom', ''))
  
  if(file.exists(str_c(base_name, '.lfmm')) & file.exists(str_c(base_name, '.geno'))){
    lfmm_file <- str_c(base_name, '.lfmm')
    geno_file <- str_c(base_name, '.geno')
  } else {
    if(model_family == 'binomial'){
      genomic_data <- as.matrix(genomic_data) %>% 
        is_greater_than(0) %>%
        apply(2, as.integer) %>%
        set_rownames(rownames(genomic_data))
    }
    
    lfmm_file <- write.lfmm(as.matrix(genomic_data), str_c(base_name, '.lfmm'))
    geno_file <- write.geno(as.matrix(genomic_data), str_c(base_name, '.geno'))
  }
  
  tibble(lfmm_file = lfmm_file, geno_file = geno_file)
}

run_pca <- function(lfmm_file){
  
  if(file.exists(str_replace(lfmm_file, '\\.lfmm$', '.pca'))){
    pca_out <- load.pcaProject(str_replace(lfmm_file, '\\.lfmm$', '.pcaProject'))
  } else {
    pca_out <- pca(lfmm_file, center = TRUE, scale = TRUE)
  }
  pca_out
}

plot_pca <- function(pca_file, region, point_size = 1){
  if(is.na(region)){region <- c('Panama', 'Florida')}
  pct_explained <- pca_file$eigenvalue[,1] / sum(pca_file$eigenvalue[,1])
  
  as_tibble(pca_file$projections) %>%
    select(V1:V5) %>%
    # bind_cols(tibble(sample = colnames(snp_data)[-1:-4])) %>%
    bind_cols(filter(meta_data, location %in% region)) %>%
    ggplot(aes(x = V1, y = V2, shape = location, colour = disease_resistance)) +
    geom_point(size = point_size) +
    labs(x = str_c('PC1 (', round(pct_explained[1]*100, 2), '%)'),
         y = str_c('PC2 (', round(pct_explained[2]*100, 2), '%)'),
         colour = 'Disease\nResistance',
         shape = NULL) +
    theme_classic()
}

run_structure <- function(geno_file, k_choices, N, calc_entropy){
  if(file.exists(str_replace(geno_file, '\\.geno$', '.snmf'))){
    unlinked_structure <- load.snmfProject(str_replace(geno_file, '\\.geno$', '.snmfProject'))
  } else {
    unlinked_structure <- snmf(geno_file,
                               K = k_choices, # normally we would run many values of K
                               entropy = calc_entropy,
                               repetitions = N, #normally we would do 10 at least
                               project = "new")
  }
  unlinked_structure
}

cross_ent_plot <- function(genetic_structure){
  map(unique(genetic_structure@K), cross.entropy, object = genetic_structure) %>%
    set_names(unique(genetic_structure@K)) %>%
    map(as_tibble) %>%
    map(~rename(.x, 'cross_entropy' = 1)) %>%
    bind_rows(.id = 'K') %>%
    mutate(K = as.integer(K)) %>%
    ggplot(aes(x = K, y = cross_entropy)) +
    geom_violin(aes(group = K)) +
    stat_summary(fun.data = mean_cl_boot) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
    theme_classic()
}

structure_plot <- function(genetic_structure, choose_K = 2, region){
  if(is.na(region)){region <- c('Panama', 'Florida')}
  run_choice <- ifelse(is.null(cross.entropy(genetic_structure, K = choose_K)), 
                       1L, 
                       which.min(cross.entropy(genetic_structure, K = choose_K)))
  
  Q(genetic_structure, K = choose_K, run = run_choice) %>%
    as_tibble() %>%
    rename_with(~str_replace(., '^V', 'Cluster_')) %>%
    bind_cols(filter(meta_data, location %in% region)) %>%
    rowwise() %>%
    mutate(admix_assignment = which.max(c_across(starts_with('Cluster'))) %>% 
             str_c('Cluster', .)) %>%
    mutate(max = max(c_across(starts_with('Cluster'))),
           match = which.max(c_across(starts_with('Cluster')))) %>%
    ungroup %>%
    arrange(location, match, 
            -max) %>%
    mutate(number = row_number()) %>%
    pivot_longer(cols = starts_with('Cluster')) %>%
    ggplot(aes(x = number, y = value, fill = name)) +
    geom_bar(position="stack", stat="identity") +
    facet_grid( ~ location, scales = 'free_x', space = 'free_x') +
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
}

#### Read in Data ####
meta_data <- read_csv('../intermediate_files/preprocessed_metadata.csv', 
                      show_col_types = FALSE) %>%
  filter(!is.na(disease_resistance))

genetic_data <- expand_grid(pop_group = c('acerv'),
            only_unlinked = c(TRUE, FALSE)) %>%
  mutate(filename = str_c(if_else(only_unlinked, 
                                  '../intermediate_files/unlinked_preprocessed_',
                                  '../intermediate_files/preprocessed_'), 
                          pop_group, '.rds')) %>%
  
  rowwise() %>%
  mutate(genomic_data = list(read_rds(filename))) %>%
  ungroup 

#### Run interspecific ANGSD on used subset and unlinked loci ####
if(file.exists('../variant_calling/subset_interspecific/sample_structure_pca_data.csv')){
  angsd_results_interspecific <- read_csv('../variant_calling/subset_interspecific/sample_structure_pca_data.csv', show_col_types = FALSE)
} else {
  dir.create('../variant_calling/subset_interspecific')
  kept_samples <- read_csv('../intermediate_files/clone_metadata.csv', show_col_types = FALSE) %>%
    group_by(clone_group) %>%
    mutate(keep = pct_missing == min(pct_missing)) %>%
    # filter(pct_missing == min(pct_missing)) %>% 
    ungroup %>%
    mutate(keep = keep & ((species == 'Ac' & data_origin == 'vollmer' & ID %in% meta_data$ID) | species == 'Apa' | species == 'Apr')) %>%
    mutate(keep = keep & ((species == 'Ac' & pct_missing  < 0.3) | species == 'Apa' | species == 'Apr')) %>%
    select(keep) %>%
    mutate(keep = as.integer(keep)) %T>%
    write_csv('../variant_calling/subset_interspecific/subset_data.txt', col_names = FALSE)
  
  kept_loci <- read_csv('../variant_calling/genotyping/genotypes.csv', 
                        col_types = cols(.default = col_integer(),
                                         contig = col_character(),
                                         ref = col_character(),
                                         alt = col_character())) %>%
    select(contig, position) %>%
    mutate(locus = str_c(contig, position, sep = '-'), .keep = 'unused') %>%
    full_join(genetic_data %>%
                filter(only_unlinked) %>%
                select(genomic_data) %>%
                pull(genomic_data) %>%
                pluck(1) %>%
                as.matrix %>%
                colnames %>%
                tibble(locus = .) %>%
                mutate(keep = TRUE),
              by = 'locus') %>%
    mutate(keep = if_else(is.na(keep), FALSE, keep)) %>%
    select(keep) %>%
    mutate(keep = as.integer(keep)) %T>%
    write_csv('../variant_calling/subset_interspecific/subset_loci.txt', col_names = FALSE)

  #Copy subset_loci and subset_data to Discovery
  # srun -t 24:00:00 --nodes=1 --cpus-per-task=40 --pty /bin/bash
  # source activate pcangsd
  # scriptDir=$(pwd)/bash_code
  # outdir=$(pwd)/variant_calling/subset_interspecific
  # cd ${outdir}
  
  # ln -s ${outdir}/../genotyping/genolike.beagle.gz
  # ln -s ${outdir}/../genotyping/bam.list
  # ln -s ${outdir}/../genotyping/genolike.geno.gz
  
  # pcangsd \
  #   -b genolike.beagle.gz \
  #   -t ${SLURM_CPUS_PER_TASK} \
  #   -o output.pcangsd \
  #   --filter subset_data.txt \
  #   --filterSites subset_loci.txt \
  #   --minMaf 0.05 \
  #   --tree \
  #   --maf_save \
  #   --pi_save \
  #   --admix \
  #   --admix_K 2 \
  #   --admix_auto 10
  
  # bash ${scriptDir}/runRscript.slurm \
  #   ${scriptDir}/r_utils/post_genotyping.R \
  #     $(pwd)
  
  angsd_results_interspecific <- read_csv('../variant_calling/subset_interspecific/sample_structure_pca_data.csv', show_col_types = FALSE)
}

angsd_results_interspecific %>%
  mutate(species_location = if_else(species == 'Ac', region, species)) %>%
  ggplot(aes(x = PC1, y = PC2, colour = species_location)) +
  geom_point(size = 2) +
  scale_colour_discrete(labels = c(expression(italic("A. palmata")),
                                   expression(italic("A. prolifera")),
                                   'Florida',
                                   'Panama')) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(x = 'PC1 (12.0%)',
       y = 'PC2 (4.3%)',
       colour = NULL) +
  theme_classic() +
  theme(legend.title = element_text(face = 'plain', colour = 'black', size = 20),
        legend.text = element_text(face = 'plain', colour = 'black', size = 20),
        axis.text = element_text(face = 'plain', colour = 'black', size = 24),
        axis.title = element_text(face = 'plain', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA)) 
ggsave('../Results/unlinked_interspecific_pca.png', height = 7, width = 9)



angsd_results_interspecific %>%
  select(-starts_with('PC')) %>%
  mutate(species_location = if_else(species == 'Ac', region, species),
         species_location = case_when(species_location == 'Apa' ~ "italic('A. palmata')",
                                      species_location == 'Apr' ~ "italic('A. prolifera')",
                                      species_location == 'FL' ~ 'Florida',
                                      species_location == 'PA' ~ 'Panama'),
         species_location = factor(species_location, 
                                   levels = c("italic('A. palmata')",
                                              "italic('A. prolifera')",
                                              "Florida",
                                              "Panama"))) %>%
  group_by(species_location) %>%
  arrange(match, 
          -max) %>%
  mutate(number = row_number()) %>%
  ungroup %>%
  pivot_longer(cols = starts_with('Cluster')) %>%
  ggplot(aes(x = number, y = value, fill = name)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid( ~ species_location, scales = 'free_x', labeller = label_parsed) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL,
       y = 'Assignment Probability',
       fill = NULL) +
  theme_classic() +
  theme(legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "vertical",
        legend.title = element_text(face = 'plain', colour = 'black', size = 20),
        legend.text = element_text(face = 'plain', colour = 'black', size = 20),
        axis.text = element_text(face = 'plain', colour = 'black', size = 24),
        axis.text.x = element_blank(),
        axis.title = element_text(face = 'plain', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = 'transparent'),
        strip.text = element_text(face = 'plain', colour = 'black', size = 24),
        strip.background = element_blank())
ggsave('../Results/unlinked_interspecific_structure_plot.png', height = 7, width = 9)


angsd_results_interspecific %>%
  left_join(meta_data, by = 'ID') %>%
  select(ID, region, disease_resistance, PC1) %>%
  filter(!is.na(disease_resistance)) %>%
  ggplot(aes(x = PC1, y = disease_resistance, colour = region)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y ~ x, colour = 'black') +
  scale_color_discrete(labels = c('Florida', 'Panama')) +
  guides(colour = list(size = 3)) +
  labs(x = 'PC1 (12.0%)',
       y = 'Disease Resistance',
       colour = NULL) +
  theme_classic() +
  theme(legend.title = element_text(face = 'plain', colour = 'black', size = 20),
        legend.text = element_text(face = 'plain', colour = 'black', size = 20),
        axis.text = element_text(face = 'plain', colour = 'black', size = 24),
        axis.title = element_text(face = 'plain', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA)) 
ggsave('../Results/unlinked_interspecific_PC1_vs_disease_resistance.png', height = 7, width = 9)

angsd_results_interspecific %>%
  left_join(meta_data, by = 'ID') %>%
  select(ID, region, disease_resistance, PC1) %>%
  filter(!is.na(disease_resistance)) %>%
  lm(disease_resistance ~ PC1, data = .) %T>%
  {summary(.) %>% print} %>%
  anova

#### Run intraspecific ANGSD on used subset and unlinked loci ####
if(file.exists('../variant_calling/subset_intraspecific/sample_structure_pca_data.csv')){
  angsd_results_intraspecific <- read_csv('../variant_calling/subset_intraspecific/sample_structure_pca_data.csv', show_col_types = FALSE)
} else {
  dir.create('../variant_calling/subset_intraspecific')
  kept_samples <- read_csv('../intermediate_files/clone_metadata.csv', show_col_types = FALSE) %>%
    group_by(clone_group) %>%
    mutate(keep = pct_missing == min(pct_missing)) %>%
    # filter(pct_missing == min(pct_missing)) %>% 
    ungroup %>%
    mutate(keep = keep & ((species == 'Ac' & data_origin == 'vollmer' & ID %in% meta_data$ID))) %>%
    mutate(keep = keep & ((species == 'Ac' & pct_missing  < 0.3))) %>%
    select(keep) %>%
    mutate(keep = as.integer(keep)) %T>%
    write_csv('../variant_calling/subset_intraspecific/subset_data.txt', col_names = FALSE)
  
  kept_loci <- read_csv('../variant_calling/genotyping/genotypes.csv', 
                        col_types = cols(.default = col_integer(),
                                         contig = col_character(),
                                         ref = col_character(),
                                         alt = col_character())) %>%
    select(contig, position) %>%
    mutate(locus = str_c(contig, position, sep = '-'), .keep = 'unused') %>%
    full_join(genetic_data %>%
                filter(only_unlinked) %>%
                select(genomic_data) %>%
                pull(genomic_data) %>%
                pluck(1) %>%
                as.matrix %>%
                colnames %>%
                tibble(locus = .) %>%
                mutate(keep = TRUE),
              by = 'locus') %>%
    mutate(keep = if_else(is.na(keep), FALSE, keep)) %>%
    select(keep) %>%
    mutate(keep = as.integer(keep)) %T>%
    write_csv('../variant_calling/subset_intraspecific/subset_loci.txt', col_names = FALSE)
  
  #Copy subset_loci and subset_data to Discovery
  # srun -t 24:00:00 --nodes=1 --cpus-per-task=40 --pty /bin/bash
  # source activate pcangsd
  # scriptDir=$(pwd)/bash_code
  # outdir=$(pwd)/variant_calling/subset_intraspecific
  # cd ${outdir}
  
  # ln -s ${outdir}/../genotyping/genolike.beagle.gz
  # ln -s ${outdir}/../genotyping/bam.list
  # ln -s ${outdir}/../genotyping/genolike.geno.gz
  
  
  # pcangsd \
  #   -b genolike.beagle.gz \
  #   -t ${SLURM_CPUS_PER_TASK} \
  #   -o output.pcangsd \
  #   --filter subset_data.txt \
  #   --filterSites subset_loci.txt \
  #   --minMaf 0.05 \
  #   --tree \
  #   --maf_save \
  #   --pi_save \
  #   --admix \
  #   --admix_K 2 \
  #   --admix_auto 10
  
  # bash ${scriptDir}/runRscript.slurm \
  #   ${scriptDir}/r_utils/post_genotyping.R \
  #     $(pwd)
  
  angsd_results_intraspecific <- read_csv('../variant_calling/subset_intraspecific/sample_structure_pca_data.csv', show_col_types = FALSE)
}

wesanderson::wes_palette("Zissou1", 9, type = "continuous")

location_colour <- set_names(wesanderson::wes_palette("Zissou1", 9, type = "continuous")[c(2, 8)],
                             c('FL', 'PA'))

pca_plot <- angsd_results_intraspecific %>%
  left_join(meta_data, by = 'ID') %>%
  mutate(panel = 'A') %>%
  ggplot(aes(x = PC1, y = PC2, colour = region, shape = region)) +
  geom_point(size = 2) +
  scale_shape_discrete(labels = c('Florida', 'Panama')) +
  scale_colour_manual(values = location_colour,
                        labels = c('Florida', 'Panama')) +
  # guides(shape = guide_legend(override.aes = list(size = 3))) +
  facet_grid(~panel, switch = 'x') +
  guides(shape = 'none',
         colour = 'none') +
  labs(x = 'PC1 (6.7%)',
       y = 'PC2 (2.2%)',
       colour = NULL,
       shape = NULL) +
  theme_classic() +
  theme(#aspect.ratio = 1,
        legend.title = element_text(face = 'plain', colour = 'black', size = 20),
        legend.text = element_text(face = 'plain', colour = 'black', size = 20),
        axis.text = element_text(face = 'plain', colour = 'black', size = 10),
        axis.title = element_text(face = 'plain', colour = 'black', size = 16),
        panel.border = element_rect(colour = 'black', fill = NA),
        strip.placement = 'outside',
        strip.text = element_blank(),
        strip.background = element_blank()) 
ggsave('../Results/unlinked_intraspecific_pca.png', plot = pca_plot, height = 7, width = 9)

structure_plot <- angsd_results_intraspecific %>%
  select(-starts_with('PC')) %>%
  mutate(species_location = if_else(species == 'Ac', region, species),
         species_location = case_when(species_location == 'Apa' ~ "italic('A. palmata')",
                                      species_location == 'Apr' ~ "italic('A. prolifera')",
                                      species_location == 'FL' ~ 'Florida',
                                      species_location == 'PA' ~ 'Panama'),
         species_location = factor(species_location, 
                                   levels = c("italic('A. palmata')",
                                              "italic('A. prolifera')",
                                              "Florida",
                                              "Panama"))) %>%
  
  mutate(match = if_else(region == 'FL', 2, 1)) %>%
  
  group_by(species_location) %>%
  arrange(match, 
          -max) %>%
  mutate(number = row_number()) %>%
  mutate(number = if_else(region == 'PA', n():1, number)) %>%
  # arrange(-number) %>%
  # mutate(number = row_number()) %>%
  ungroup %>%
  pivot_longer(cols = starts_with('Cluster')) %>%
  mutate(name = case_when(name == 'Cluster1' ~ 'PA',
                          name == 'Cluster2' ~ 'FL')) %>%
  ggplot(aes(x = number, y = value, fill = name)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid( ~ species_location, scales = 'free_x', labeller = label_parsed,
              switch = 'x') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = location_colour) +
  guides(fill = 'none') +
  labs(x = NULL,
       y = 'Assignment Probability',
       fill = NULL) +
  theme_classic() +
  theme(#aspect.ratio = 1,
        legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "vertical",
        legend.title = element_text(face = 'plain', colour = 'black', size = 20),
        legend.text = element_text(face = 'plain', colour = 'black', size = 20),
        axis.text = element_text(face = 'plain', colour = 'black', size = 10),
        axis.text.x = element_blank(),
        # axis.title = element_text(face = 'plain', colour = 'black', size = 16),
        axis.title = element_blank(),
        panel.border = element_rect(colour = 'black', fill = 'transparent'),
        strip.text = element_text(face = 'plain', colour = 'black', size = 16),
        strip.placement = "outside", 
        # strip.text = element_blank(),
        strip.background = element_blank())
ggsave('../Results/unlinked_intraspecific_structure_plot.png', plot = structure_plot, height = 7, width = 9)

pca_plot / structure_plot + 
  plot_layout(heights = c(1, 0.3)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18))
ggsave('../Results/unlinked_intraspecific_plots.png', height = 7, width = 9)

angsd_results_intraspecific %>%
  left_join(meta_data, by = 'ID') %$%
  cor.test((1 - max), disease_resistance, method = 'kendall')

angsd_results_intraspecific %>%
  left_join(meta_data, by = 'ID') %$%
  cor.test(PC1, disease_resistance, method = 'kendall')

angsd_results_intraspecific %>%
  left_join(meta_data, by = 'ID') %>%
  ggplot(aes(x = rank(1 - max), y = disease_resistance, colour = region, shape = region)) +
  geom_point() +
  scale_colour_manual(values = location_colour, labels = c('Florida', 'Panama')) +
  scale_shape_manual(values = c('FL' = 'circle', 'PA' = 'triangle'), labels = c('Florida', 'Panama')) +
  labs(x = 'Rank-Order Admixture',
       y = 'Disease Resistance',
       colour = NULL,
       shape = NULL) +
  theme_classic() +
  theme(legend.text = element_text(face = 'plain', colour = 'black', size = 20),
    axis.text = element_text(face = 'plain', colour = 'black', size = 10),
    axis.title = element_text(face = 'plain', colour = 'black', size = 16),
    panel.border = element_rect(colour = 'black', fill = 'transparent'),
    strip.text = element_text(face = 'plain', colour = 'black', size = 16),
    strip.placement = "outside", 
    # strip.text = element_blank(),
    strip.background = element_blank())
ggsave('../Results/admixture_v_disease.png', height = 7, width = 7)

#### Write LFMM Necessary Files ####
mkdir <- map_lgl(c('../lea_projects'), ~dir.create(.x, showWarnings = FALSE))
setwd('../lea_projects')

#write genotype files
lea_data <- genetic_data %>%
  rowwise(pop_group, only_unlinked) %>%
  summarise(write_lea_genetics(genomic_data, filename, model_family = 'norm'), 
            .groups = 'drop')

#### Run Genetic PCA ####
pca_out <- lea_data %>%
  ungroup %>%
  filter(only_unlinked) %>%
  mutate(location = str_extract(pop_group, 'panama|florida') %>% str_to_sentence) %>%
  rowwise(pop_group) %>%
  summarise(pca_out = list(run_pca(lfmm_file)),
            plot = list(plot_pca(pca_out, location, 3)),
            .groups = 'drop')
  
wrap_plots(pca_out$plot)
pca_out$plot[[1]] +
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  theme(aspect.ratio = 1,
        legend.title = element_text(face = 'bold', colour = 'black', size = 20),
        legend.text = element_text(face = 'bold', colour = 'black', size = 20),
        axis.text = element_text(face = 'bold', colour = 'black', size = 24),
        axis.title = element_text(face = 'bold', colour = 'black', size = 24),
        panel.border = element_rect(colour = 'black', fill = NA)) 
ggsave('../Results/unlinked_pca.png', height = 10, width = 10)

#### Run Structure ####
structure_run <- lea_data %>%
  ungroup %>%
  filter(only_unlinked) %>%
  select(pop_group, geno_file) %>%
  rowwise(pop_group) %>%
  summarise(structure_out = list(run_structure(geno_file, k_choices = c(1:10), N = 10, calc_entropy = TRUE)),
            cross_plot = list(cross_ent_plot(structure_out)),
            .groups = 'drop')

structure_run %>%
  pull(cross_plot) %>%
  wrap_plots() +
  theme(panel.background = element_rect(colour = 'black', linewidth = 1),
        axis.text = element_text(colour = 'black', size = 12),
        axis.title = element_text(colour = 'black', size = 16)) +
  labs(x = 'K',
       y = 'Cross Entropy')
ggsave('../Results/unlinked_Choosing_K.png', height = 5, width = 5)

structure_plots <- structure_run %>%
  mutate(location = str_extract(pop_group, 'panama|florida') %>% str_to_sentence,
         K = c(2)) %>%
  rowwise(pop_group) %>%
  summarise(plots = list(structure_plot(structure_out, K, location)),
            .groups = 'drop')

wrap_plots(structure_plots$plots)
ggsave('../Results/unlinked_structure_plot.png', height = 5, width = 10)

#### Run best K once for linked (binomial and nonbinomial) for imputation purposes
linked_model_out <- lea_data %>%
  ungroup %>%
  filter(!only_unlinked) %>%
  select(pop_group, geno_file) %>%
  mutate(location = str_extract(pop_group, 'panama|florida') %>% str_to_sentence,
         K = c(2)) %>%
  rowwise(pop_group) %>%
  summarise(structure_out = list(run_structure(geno_file, k_choices = K, N = 1, calc_entropy = FALSE)),
            structure_plot = list(structure_plot(structure_out, K, location) + labs(title = pop_group)),
            .groups = 'drop')

linked_model_out %>%
  pull(structure_plot) %>%
  wrap_plots()
ggsave('../Results/linked_structure_plot.png', height = 5, width = 10)

          