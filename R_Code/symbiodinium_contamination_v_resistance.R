library(tidyverse)
library(patchwork)
library(janitor)
library(lme4)
library(lmerTest)

location_colour <- set_names(c(wesanderson::wes_palette("Zissou1", 9, type = "continuous")[c(2, 8)]),
                             c('Florida', 'Panama'))

data <- read_csv('../../intermediate_files/preprocessed_metadata.csv', show_col_types = FALSE) %>%
  
  inner_join(read_csv('../../../Bioinformatics/variant_calling/18October2022/fastq_screen_results_summary.csv',
                      show_col_types = FALSE) %>%
               clean_names %>%
               select(-starts_with('percent')) %>%
               group_by(id, genome, symbiodinium_species) %>%
               summarise(across(where(is.numeric), mean), .groups = 'drop') %>%
               mutate(across(where(is.numeric), floor)) %>%
               ungroup %>%
               nest(contam_data = -c(id)),
            
            by = c('ID' = 'id')) %>%
  mutate(species_location = if_else(data_origin == 'vollmer', location, species))


#### Plot & Analysis for Overall amount of Symbiodinium Contamination in Genome ####
percent_plot <- data %>%
  unnest(contam_data) %>%
  filter(genome == 'Symbiodinium') %>%
  mutate(percent_mapped = 1 - number_unmapped / number_reads_processed) %>%
  filter(!is.na(disease_resistance)) %>%
  ggplot(aes(x = disease_resistance, y = percent_mapped, colour = location, shape = location)) +
  geom_point() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = location_colour, 
                     labels = c('Florida' = 'Florida', 'Panama' = 'Panama')) +
  scale_shape_manual(values = c('Florida' = 'circle', 'Panama' = 'triangle') , 
                     labels = c('Florida' = 'Florida', 'Panama' = 'Panama')) +
  labs(x = 'Disease Resistance',
       y = 'Symbiodinium sp. Content (%)',
       colour = NULL,
       shape = NULL) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 16))

data %>%
  unnest(contam_data) %>%
  filter(genome == 'Symbiodinium') %>%
  mutate(across(c(number_unmapped, number_reads_processed), round)) %>%
  mutate(percent_mapped = 1 - number_unmapped / number_reads_processed,
         ID = factor(ID)) %>%
  filter(!is.na(disease_resistance)) %>% 
  glmer(percent_mapped ~ disease_resistance * location + (1 | ID), data = ., family = 'binomial', 
        weights = number_reads_processed) %>% summary

#### Plot and Analysis for major Clade Composition - B & C ####

symb_diversity <- data %>%
  unnest(contam_data) %>%
  filter(!is.na(symbiodinium_species)) %>%
  mutate(n_mapped = number_one_hit_one_genome + number_multiple_hits_one_genome + number_one_hit_multiple_genomes + multiple_hits_multiple_genomes) %>%
  select(ID:disease_resistance, species_location, symbiodinium_species, number_reads_processed, n_mapped) %>%
  pivot_wider(names_from = symbiodinium_species,
              values_from = n_mapped,
              values_fn = mean) %>%
  filter(number_reads_processed > 1000)



symb_diversity %>%
  pivot_longer(cols = -c(ID:number_reads_processed)) %>%
  mutate(name = case_when(str_detect(name, 'microadriaticum') ~ 'A',
                          str_detect(name, 'pilosum') ~ 'A',
                          str_detect(name, 'natans') ~ 'A',
                          str_detect(name, 'necroappetens') ~ 'A',
                          str_detect(name, 'Breviolum') ~ 'B',
                          str_detect(name, 'clade A') ~ 'A',
                          str_detect(name, 'clade C') ~ 'C',
                          str_detect(name, 'kawagutii') ~ 'F',
                          TRUE ~ name)) %>%
  group_by(ID, species_location, name) %>%
  summarise(value = mean(value),
            .groups = 'drop') %>%
  # filter(str_detect(name, 'Breviolum|microadriaticum')) %>%
  ggplot(aes(x = ID, y = value, fill = name)) +
  geom_col(position="fill") +
  facet_grid(~species_location, scales = 'free_x', space = 'free_x')

symb_diversity %>%
  pivot_longer(cols = -c(ID:number_reads_processed))  %>%
  filter(str_detect(name, 'Breviolum|microadriaticum'),
         data_origin == 'vollmer') %>%
  select(ID, location, name, value) %>%
  group_by(ID, location) %>%
  mutate(percent = value / sum(value)) %>%
  ungroup %>%
  filter(str_detect(name, 'microadriaticum')) %>%
  t.test(percent ~ location, data = .)
  
  
library(afex)
diversity_data <- symb_diversity %>%
  pivot_longer(cols = -c(ID:number_reads_processed))  %>%
  filter(str_detect(name, 'Breviolum|microadriaticum'),
         data_origin == 'vollmer') %>%
  select(ID, location, disease_resistance, name, value) %>%
  group_by(ID, location, disease_resistance) %>%
  mutate(percent = value / sum(value),
         value = floor(value),
         total = sum(value)) %>%
  ungroup %>%
  filter(str_detect(name, 'microadriaticum')) %>%
  mutate(across(c(ID, location), factor))

diversity_data$name

diversity_model <- mixed(cbind(value, total - value) ~ location + (1 | ID),
        data = diversity_data, 
        family = 'binomial',
        method = 'LRT')

library(emmeans)
library(ggdist)
library(gghalves)
library(patchwork)
library(ggtext)

composition_plot <- emmeans(diversity_model, ~location, type = 'response') %>%
  as_tibble() %>%
  ggplot(aes(x = location, y = prob, ymin = asymp.LCL, ymax = asymp.UCL, colour = location, shape = location)) +
  stat_halfeye(data = diversity_data, inherit.aes = FALSE,
               aes(x = location, y = value / total, colour = location, fill = location),
               adjust = 0.5, width = 0.6, .width = 0, alpha = 0.5, show.legend = FALSE,
               justification = -0.155, show_point = FALSE) +
  geom_half_point(data = diversity_data, inherit.aes = FALSE,
                  aes(x = location, y = value / total, colour = location, fill = location, shape = location),
                  side = 'l', range_scale = 0.1, alpha = 1,
                  transformation = position_jitter(height = 0, width = 0.1),
                  show.legend = FALSE) +
  
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = location_colour , 
                     labels = c('Florida' = 'Florida', 'Panama' = 'Panama')) +
  scale_fill_manual(values = location_colour , 
                     labels = c('Florida' = 'Florida', 'Panama' = 'Panama')) +
  scale_shape_manual(values = c('Florida' = 'circle', 'Panama' = 'triangle') , 
                    labels = c('Florida' = 'Florida', 'Panama' = 'Panama')) +
  labs(y = "*Symbiodinium*/*Breviolum* (%)",
       x = NULL,
       colour = NULL,
       fill = NULL) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        axis.title.y = element_markdown(),
        legend.text = element_text(colour = 'black', size = 16))



percent_plot + composition_plot + 
  plot_layout(guides = 'collect', widths = c(0.5, 0.5)) + 
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 24, colour = 'black'))
ggsave('../../Results/symbiodinium content.png', height = 5, width = 10)

ggsave('../../Results/symbiodinium content.png', plot = composition_plot, height = 5, width = 5)



diversity_model <- mixed(cbind(value, total - value) ~ location * disease_resistance + (1 | ID),
                         data = diversity_data, 
                         family = 'binomial',
                         method = 'LRT')

diversity_model
summary(diversity_model)

diversity_model_noID <- glm(cbind(value, total - value) ~ location * disease_resistance,
                         data = diversity_data, 
                         family = 'binomial')
anova(diversity_model_noID, diversity_model, test = 'LRT')

logLik(diversity_model_noID)
logLik(diversity_model$full_model)
lmtest::lrtest(diversity_model_noID, diversity_model$full_model)

(d0 <- deviance(diversity_model_noID))
(d1 <- deviance(diversity_model$full_model))
(LR <- d0 - d1)
pchisq(LR, 1, lower = FALSE)

diversity_v_resistance_plot <- emmeans(diversity_model, ~location * disease_resistance, 
        at = list(disease_resistance = modelr::seq_range(diversity_data$disease_resistance, 25)),
        type = 'response') %>%
  as_tibble() %>%
  ggplot(aes(x = disease_resistance, y = prob, 
             colour = location, fill = location,
             shape = location,
             ymin = asymp.LCL, ymax = asymp.UCL)) +
  geom_ribbon(colour = NA, alpha = 0.5, show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  geom_point(data = diversity_data, aes(x = disease_resistance, y = value / total, colour = location, shape = location),
             inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = location_colour , 
                     labels = c('Florida' = 'Florida', 'Panama' = 'Panama')) +
  scale_fill_manual(values = location_colour , 
                    labels = c('Florida' = 'Florida', 'Panama' = 'Panama')) +
  scale_shape_manual(values = c('Florida' = 'circle', 'Panama' = 'triangle') , 
                    labels = c('Florida' = 'Florida', 'Panama' = 'Panama')) +
  # guides(shape = 'none') +
  labs(colour = NULL,
       fill = NULL,
       x = 'Disease Resistance', 
       y = "*Symbiodinium*/*Breviolum* (%)",
       shape = NULL) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        axis.title.y = element_markdown(),
        legend.text = element_text(colour = 'black', size = 16))

composition_plot + diversity_v_resistance_plot + 
  plot_layout(guides = 'collect', widths = c(0.5, 0.5)) + 
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 24, colour = 'black'))
ggsave('../../Results/symbiodinium content.png', height = 5, width = 10)



#### Extra post review ####

symb_comp <- data %>%
  unnest(contam_data) %>%
  filter(!is.na(symbiodinium_species)) %>%
  mutate(unique_map = number_one_hit_one_genome + number_multiple_hits_one_genome) %>%
  select(ID, gen_id, species_location, disease_resistance, symbiodinium_species, number_reads_processed, unique_map) %>%
  arrange(ID) %>%
  mutate(symbiodinium_clade = case_when(str_detect(symbiodinium_species, 'microadriaticum') ~ 'A',
                                        str_detect(symbiodinium_species, 'pilosum') ~ 'A',
                                        str_detect(symbiodinium_species, 'natans') ~ 'A',
                                        str_detect(symbiodinium_species, 'necroappetens') ~ 'A',
                                        str_detect(symbiodinium_species, 'Breviolum') ~ 'B',
                                        str_detect(symbiodinium_species, 'clade A') ~ 'A',
                                        str_detect(symbiodinium_species, 'clade C') ~ 'C',
                                        str_detect(symbiodinium_species, 'kawagutii') ~ 'F',
                                        str_detect(symbiodinium_species, 'CCMP') ~ 'A',
                                        str_detect(symbiodinium_species, 'KB8') ~ 'A',
                                        TRUE ~ symbiodinium_species)) %>%
  mutate(symbiodinium_clade = str_replace_all(symbiodinium_clade,
                                              c('A' = 'Symbiodinium',
                                                'B' = 'Breviolum',
                                                'C' = 'Cladocopium', 
                                                'F' = 'Fugacium')),
         symbiodinium_clade = factor(symbiodinium_clade, 
                                     levels = c('Symbiodinium', 'Breviolum', 
                                                'Cladocopium', 'Fugacium'))) %>% 
  group_by(ID, gen_id, species_location, disease_resistance, 
           number_reads_processed, symbiodinium_clade) %>%
  summarise(unique_map = sum(unique_map),
            .groups = 'drop') %>%
  group_by(ID) %>%
  mutate(total = sum(unique_map),
         percent_symbiont = unique_map / total) %>%
  ungroup %>%
  filter(species_location %in% c('Florida', 'Panama')) %>%
  
  ggplot(aes(x = ID, y = percent_symbiont, fill = symbiodinium_clade)) +
  geom_col(position = position_fill(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(x = NULL,
       y = 'Read Composition (%)',
       fill = 'Symbiont Clade') +
  guides(fill = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  facet_grid(~species_location, scales = 'free_x', space = 'free_x',
             switch = "x") +
  theme_classic() +
  theme(legend.position = 'bottom',
        panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 14),
        strip.background = element_blank(),
        strip.placement = 'outside',
        strip.text = element_text(colour = 'black', size = 16),
        axis.text.x = element_blank())


symb_comp



nmds_data <- data %>%
  unnest(contam_data) %>%
  filter(!is.na(symbiodinium_species)) %>%
  mutate(unique_map = number_one_hit_one_genome + number_multiple_hits_one_genome) %>%
  select(ID, gen_id, species_location, disease_resistance, symbiodinium_species, number_reads_processed, unique_map) %>%
  arrange(ID) %>%
  mutate(symbiodinium_clade = case_when(str_detect(symbiodinium_species, 'microadriaticum') ~ 'A',
                                        str_detect(symbiodinium_species, 'pilosum') ~ 'A',
                                        str_detect(symbiodinium_species, 'natans') ~ 'A',
                                        str_detect(symbiodinium_species, 'necroappetens') ~ 'A',
                                        str_detect(symbiodinium_species, 'Breviolum') ~ 'B',
                                        str_detect(symbiodinium_species, 'clade A') ~ 'A',
                                        str_detect(symbiodinium_species, 'clade C') ~ 'C',
                                        str_detect(symbiodinium_species, 'kawagutii') ~ 'F',
                                        str_detect(symbiodinium_species, 'CCMP') ~ 'A',
                                        str_detect(symbiodinium_species, 'KB8') ~ 'A',
                                        TRUE ~ symbiodinium_species)) %>%
  mutate(symbiodinium_clade = str_replace_all(symbiodinium_clade,
                                              c('A' = 'Symbiodinium',
                                                'B' = 'Breviolum',
                                                'C' = 'Cladocopium', 
                                                'F' = 'Fugacium'))) %>% 
  group_by(ID, gen_id, species_location, disease_resistance, number_reads_processed, symbiodinium_clade) %>%
  summarise(unique_map = sum(unique_map),
            .groups = 'drop') %>%
  group_by(ID) %>%
  mutate(total = sum(unique_map),
         percent_symbiont = unique_map / total) %>%
  ungroup %>%
  
  filter(species_location %in% c('Florida', 'Panama')) %>%
  select(-number_reads_processed, -total, -percent_symbiont) %>%
  pivot_wider(names_from = symbiodinium_clade, values_from = unique_map)


library(vegan)
the_nmds <- select(nmds_data, -gen_id:-disease_resistance) %>%
  column_to_rownames('ID') %>%
  metaMDS()
plot(the_nmds)

adonis2(select(nmds_data, -gen_id:-disease_resistance) %>%
          column_to_rownames('ID') ~ species_location * disease_resistance, 
          data = select(nmds_data, gen_id:disease_resistance),
        by = 'margin')

adonis2(select(nmds_data, -gen_id:-disease_resistance) %>%
          column_to_rownames('ID') ~ species_location + disease_resistance, 
        data = select(nmds_data, gen_id:disease_resistance),
        by = 'margin')



location_colour <- set_names(c(wesanderson::wes_palette("Zissou1", 9, type = "continuous")[c(2, 8)]),
                             c('Florida', 'Panama'))

nmds_plot <- scores(the_nmds, 'sites') %>%
  as_tibble(rownames = 'ID') %>%
  left_join(nmds_data,
            by = 'ID') %>%
  mutate(panel = 'A') %>%
  ggplot(aes(x = NMDS1, y = NMDS2, colour = species_location, shape = species_location)) +
  geom_point() +
  geom_text(data = scores(the_nmds, 'species') %>%
              as_tibble(rownames = 'Symbiont Clade'),
            aes(x = NMDS1, y = NMDS2, label = `Symbiont Clade`), 
            inherit.aes = FALSE, size = 6, hjust = 'inward', fontface = 'italic') +
  scale_colour_manual(values = location_colour) +
  scale_shape_manual(values = c('Florida' = 'circle', 'Panama' = 'triangle')) +
  facet_grid(~panel, switch = 'x') +
  theme_classic() +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 16),
        strip.placement = 'outside',
        strip.text = element_blank(),
        strip.background = element_blank())


nmds_plot / symb_comp + 
  plot_layout(heights = c(1, 0.3)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18))
ggsave('../../Results/symb_composition.png', height = 7, width = 9)



library(brms)
library(tidybayes)
brms_data <- nmds_data %>%
  mutate(Y = cbind(Symbiodinium, Breviolum, Cladocopium, Fugacium)) %>%
  select(-c(Symbiodinium, Breviolum, Cladocopium, Fugacium)) %>%
  mutate(total = rowSums(Y))

brms_model <- brm(Y | trials(total) ~ species_location * disease_resistance,
                  family = multinomial(),
                  data = brms_data,
                  backend = 'cmdstanr',
                  chains = 4,
                  cores = 4)
summary(brms_model)


brms_data %>%
  modelr::data_grid(species_location = unique(species_location),
                    disease_resistance = modelr::seq_range(disease_resistance, 10)) %>%
  mutate(total = 1) %>%
  add_epred_draws(brms_model, re_formula = NA) %>%
  point_interval() %>%
  ungroup %>%
  
  ggplot(aes(x = disease_resistance, y = .epred, ymin = .lower, ymax = .upper, 
             colour = species_location, fill = species_location)) +
  geom_ribbon(alpha = 0.5, colour = NA) +
  geom_line() +
  facet_wrap(~.category, scales = 'free_y') +
  scale_y_continuous(labels = scales::percent_format()) +
  guides(colour = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  theme_classic() +
  theme(axis.title = element_text(colour = 'black', size = 16),
        axis.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 12),
        panel.background = element_rect(colour = 'black'),
        legend.position = 'bottom')


brms_data %>%
  modelr::data_grid(species_location = unique(species_location),
                    disease_resistance = modelr::seq_range(disease_resistance, 10)) %>%
  mutate(total = 1) %>%
  add_epred_draws(brms_model, re_formula = NA) %>%
  point_interval() %>%
  ungroup %>%
  
  ggplot(aes(x = disease_resistance, y = .epred, ymin = .lower, ymax = .upper, 
             colour = .category, fill = .category)) +
  geom_ribbon(alpha = 0.5, colour = NA) +
  geom_line() +
  facet_wrap(~species_location, scales = 'free_y') +
  scale_y_continuous(labels = scales::percent_format()) +
  guides(colour = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  theme_classic() +
  theme(axis.title = element_text(colour = 'black', size = 16),
        axis.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 12),
        panel.background = element_rect(colour = 'black'),
        legend.position = 'bottom')


brms_data %>%
  modelr::data_grid(species_location = unique(species_location),
                    disease_resistance = median(disease_resistance)) %>%
  mutate(total = 1) %>%
  add_epred_draws(brms_model, re_formula = NA) %>%
  point_interval() %>%
  ungroup %>%
  
  ggplot(aes(x = .category, y = .epred, ymin = .lower, ymax = .upper, 
             colour = species_location)) +
  geom_pointrange(position = position_dodge(0.5)) +
  facet_wrap(~.category, scales = 'free_y') +
  scale_y_continuous(labels = scales::percent_format()) +
  guides(colour = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  theme_classic() +
  theme(axis.title = element_text(colour = 'black', size = 16),
        axis.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 12),
        panel.background = element_rect(colour = 'black'),
        legend.position = 'bottom')


hypothesis(brms_model, 'muBreviolum_species_locationPanama > 0')
hypothesis(brms_model, 'muCladocopium_species_locationPanama > 0')
hypothesis(brms_model, 'muFugacium_species_locationPanama > 0')


library(emmeans)