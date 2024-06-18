
#### Libraries ####
library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)

dir.create('../Results')

#### Read in Data ####
full_tank_data <- read_csv('../Data/Combo_Census.csv', show_col_types = FALSE) %>% 
  pivot_longer(cols = where(is.numeric),
               names_to = 'day',
               values_to = 'infected', 
               names_transform = as.numeric,
               values_transform = as.integer, 
               values_drop_na = TRUE) %>%
  janitor::clean_names() %>%
  
  #Keep only diseased treatment since healthy treatment not involved in calculation of disease resistance
  # filter(treatment == 'D') %>%
  
  #Remove Tanks which crashed due to anoxia
  filter(treatment != 'Treatment') %>%
  filter(!((tank == 'H1' & location == 'PA') | (tank == 'H2' & location == 'FL') |
           (tank == 'H3' & location == 'FL'))) %>%
  
  #Recode to match with 
  mutate(location = str_replace_all(location, c('FL' = 'Florida', 'PA' = 'Panama')),
         across(c(site, genotype), ~str_replace_all(., c('Holyshit' = 'HS', 'Sebastian' = 'SR', 'Tetas' = 'Tet'))),
         site = if_else(location == 'Florida', NA_character_, site),
         genotype = str_remove(genotype, '_[0-9]+_')) %>%
  dplyr::rename(reef = site,
         gen_id = genotype) %>%
  select(-reef) %>%
  left_join(read_csv('../intermediate_files/collected_sample_metadata.csv', show_col_types = FALSE) %>%
              select(-location),
             by = 'gen_id') %>%
  mutate(dummy1 = 1,
         dummy2 = 1,
         dummy3 = 1,
         across(where(is.character), as.factor)) %>%
  select(-species, -data_origin, -ID)

#### Estimate Disease Resistance - assume location effect is purely the result of experimental differences #### 
survival_model_locationExperimental <- gam(day ~ treatment + s(location, bs = 're', by = dummy1) +
                        s(location, gen_id, bs = 're', by = dummy3) +
                        s(location, tank, bs = 're', by = dummy2), 
                      weights = infected,
                      family = cox.ph(),
                      data = full_tank_data,
                      method = 'REML',
                      control = gam.control(nthreads = 4, trace = TRUE))
summary(survival_model_locationExperimental)

survival_model_locationExperimental2 <- gam(day ~ s(location, bs = 're', by = dummy1) +
                                             s(location, gen_id, bs = 're', by = dummy3) +
                                             s(location, tank, bs = 're', by = dummy2), 
                                           weights = infected,
                                           family = cox.ph(),
                                           data = filter(full_tank_data, treatment == 'D'),
                                           method = 'REML',
                                           control = gam.control(nthreads = 4, trace = TRUE))
summary(survival_model_locationExperimental2)

#### Get Deviance Explained by Individual Components ####
survival_model_locationExperimental_noID <- gam(day ~ s(location, bs = 're', by = dummy1) +
                                                  s(location, tank, bs = 're', by = dummy2), 
                                                weights = infected,
                                                family = cox.ph(),
                                                data = filter(full_tank_data, treatment == 'D'),
                                                method = 'REML',
                                                control = gam.control(nthreads = 4, trace = TRUE))
summary(survival_model_locationExperimental_noID)



survival_model_locationExperimental_noTank <- gam(day ~ s(location, bs = 're', by = dummy1) +
                                                    s(location, gen_id, bs = 're', by = dummy3), 
                                                  weights = infected,
                                                  family = cox.ph(),
                                                  data = filter(full_tank_data, treatment == 'D'),
                                                  method = 'REML',
                                                  control = gam.control(nthreads = 4, trace = TRUE))
summary(survival_model_locationExperimental_noTank)

survival_model_locationExperimental_onlyLocation <- gam(day ~ s(location, bs = 're', by = dummy1), 
                                                        weights = infected,
                                                        family = cox.ph(),
                                                        data = filter(full_tank_data, treatment == 'D'),
                                                        method = 'REML',
                                                        control = gam.control(nthreads = 4, trace = TRUE))
summary(survival_model_locationExperimental_onlyLocation)

#### Percent var explained by genotype in each location alone ####
just_florida_full <- gam(day ~ s(gen_id, bs = 're', by = dummy3) +
                           s(tank, bs = 're', by = dummy2), 
                         weights = infected,
                         family = cox.ph(),
                         data = filter(full_tank_data, treatment == 'D',
                                       location == 'Florida'),
                         method = 'REML',
                         control = gam.control(nthreads = 4, trace = TRUE))
summary(just_florida_full)

# 14.6% variance in florida only explained by genotype

just_florida_full_noTank <- gam(day ~ s(gen_id, bs = 're', by = dummy3), 
                                weights = infected,
                                family = cox.ph(),
                                data = filter(full_tank_data, treatment == 'D',
                                              location == 'Florida'),
                                method = 'REML',
                                control = gam.control(nthreads = 4, trace = TRUE))
summary(just_florida_full_noTank)


just_florida_full_noGen <- gam(day ~ s(tank, bs = 're', by = dummy2), 
                               weights = infected,
                               family = cox.ph(),
                               data = filter(full_tank_data, treatment == 'D',
                                             location == 'Florida'),
                               method = 'REML',
                               control = gam.control(nthreads = 4, trace = TRUE))
summary(just_florida_full_noGen)




just_panama_full <- gam(day ~ s(gen_id, bs = 're', by = dummy3) +
                          s(tank, bs = 're', by = dummy2), 
                        weights = infected,
                        family = cox.ph(),
                        data = filter(full_tank_data, treatment == 'D',
                                      location == 'Panama'),
                        method = 'REML',
                        control = gam.control(nthreads = 4, trace = TRUE))
summary(just_panama_full)


just_panama_full_noTank <- gam(day ~ s(gen_id, bs = 're', by = dummy3), 
                               weights = infected,
                               family = cox.ph(),
                               data = filter(full_tank_data, treatment == 'D',
                                             location == 'Panama'),
                               method = 'REML',
                               control = gam.control(nthreads = 4, trace = TRUE))
summary(just_panama_full_noTank)

# 3.6% variance in panama only explained by genotype

just_panama_full_noGen <- gam(day ~ s(tank, bs = 're', by = dummy2), 
                              weights = infected,
                              family = cox.ph(),
                              data = filter(full_tank_data, treatment == 'D',
                                            location == 'Panama'),
                              method = 'REML',
                              control = gam.control(nthreads = 4, trace = TRUE))
summary(just_panama_full_noGen)


#### Estimate Disease Resistance ####
disease_resistance <- full_tank_data %>%
  select(gen_id, location, treatment) %>%
  distinct %>%
  mutate(dummy1 = 0,
         dummy2 = 0,
         dummy3 = 1,
         tank = 'sim',
         day = 6) %>%
  bind_cols(., predict(survival_model_locationExperimental, newdata = ., se.fit = TRUE, type = 'response')) %>%
  dplyr::rename(disease_resistance = fit) %>%
  select(-se.fit) %>%
  mutate(dummy1 = 1) %>%
  # mutate(day = 5) %>%
  bind_cols(., predict(survival_model_locationExperimental, newdata = ., se.fit = TRUE, type = 'response')) %>%
  dplyr::rename(disease_resistance_with_location = fit) %>%
  select(-se.fit, -tank, -dummy1, -dummy2, -day) %>%
  mutate(disease_resistance = as.numeric(disease_resistance),
         disease_resistance_with_location = as.numeric(disease_resistance_with_location)) %>%
  select(-disease_resistance_with_location, -dummy3)
write_csv(disease_resistance, '../intermediate_files/disease_resistance.csv')


#### Paper Plots ####
location_colour <- set_names(c('black', wesanderson::wes_palette("Zissou1", 9, type = "continuous")[c(2, 8)]),
                             c('overall', 'Florida', 'Panama'))

survival_curve <- full_tank_data %>%
  filter(treatment == 'D') %>%
  mutate(gen_id = str_c('G', location)) %>%
  select(gen_id, location, treatment) %>%
  distinct %>%
  mutate(dummy1 = 1,
         dummy2 = 0,
         dummy3 = 0,
         tank = 'sim') %>%
  add_row(gen_id = 'overall', location = 'overall', treatment = 'D', dummy1 = 0, dummy2 = 0, dummy3 = 0, tank = 'sim') %>%
  expand_grid(day = seq(0, 7, length.out = 1000)) %>%
  bind_cols(., predict(survival_model_locationExperimental, newdata = ., se.fit = TRUE, type = 'response')) %>%
  
  ggplot(aes(x = day, y = fit, colour = location,
             ymin = fit - se.fit, ymax = fit + se.fit,
             fill = location)) +
  geom_ribbon(data = . %>% filter(location != 'overall'), alpha = 0.5, colour = NA, show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  scale_color_manual(values = location_colour, 
                     labels = c('overall' = 'Standardized', 'Florida' = 'Florida', 'Panama' = 'Panama')) +
  scale_fill_manual(values = location_colour,
                    labels = c('overall' = 'Standardized', 'Florida' = 'Florida', 'Panama' = 'Panama')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(x = 'Experiment Day',
       y = 'Fragment Survival (%)',
       colour = NULL,
       fill = NULL) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 10),
        axis.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 16),
        legend.position = c(0.1, 0.1))


individual_resistance <- disease_resistance %>%
  filter(treatment == 'D') %>%
  mutate(gen_id = fct_reorder(gen_id, disease_resistance)) %>%
  ggplot(aes(x = disease_resistance, y = gen_id, colour = location)) +
  geom_point() +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_color_manual(values = location_colour, 
                     labels = c('Florida' = 'Florida', 'Panama' = 'Panama')) +
  labs(y = 'Genotype',
       x = 'Disease Resistance',
       colour = NULL) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 16),
        axis.text.y = element_blank(),
        axis.title = element_text(colour = 'black', size = 20),
        legend.text = element_text(colour = 'black', size = 16),
        legend.position = c(0.1, 0.9))


library(ggdist)
library(gghalves)
location_disease_resitance <- disease_resistance %>%
  filter(treatment == 'D') %>%
  ggplot(aes(x = location, y = disease_resistance, colour = location, fill = location, shape = location)) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, alpha = 0.5, show.legend = FALSE,
               justification = -0.15) +
  geom_boxplot(width = 0.05, outlier.shape = NA, fill = 'white', show.legend = FALSE) +
  # stat_dots(side = 'left', dotsize = 0.05, justification = 1.05, binwidth = 0.05)
  geom_half_point(side = 'l', range_scale = 0.1, alpha = 1,
                  transformation = position_jitter(height = 0, width = 0.1)) +
  scale_color_manual(values = location_colour, 
                     labels = c('overall' = 'Normalized', 'Florida' = 'Florida', 'Panama' = 'Panama')) +
  scale_fill_manual(values = location_colour, 
                    labels = c('overall' = 'Normalized', 'Florida' = 'Florida', 'Panama' = 'Panama')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(y = 'Normalized Disease Resistance',
       x = NULL,
       colour = NULL,
       fill = NULL) +
  theme_classic() +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 16),
        axis.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 16))

survival_curve + location_disease_resitance + 
  plot_layout(guides = 'collect', widths = c(0.75, 0.25)) + 
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 24, colour = 'black'))
ggsave('../Results/disease_resistance.png', height = 5, width = 10, scale = 1.1)
