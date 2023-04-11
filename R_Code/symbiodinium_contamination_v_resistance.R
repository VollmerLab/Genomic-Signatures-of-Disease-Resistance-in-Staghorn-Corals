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
