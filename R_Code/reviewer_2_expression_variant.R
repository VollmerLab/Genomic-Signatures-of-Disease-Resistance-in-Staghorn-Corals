#### Libraries ####
library(tidyverse)
library(lmerTest)

#### Data ####
full_data <- read_csv('../../Gene_Expression/intermediate_files/normalized_transcript_counts.csv.gz')

gene_names <- tibble(gene_id = c('Acer_00005575', 'Acer_00005576', 'Acer_00005577', 'Acer_00005578', 'Acer_00005579',
                                 'Acer_00031547', 'Acer_00031549', 'Acer_00031552', 'Acer_00031553', 'Acer_00031554'),
                     gene_name = c('HMCN1', 'PTPRD', 'SECG', 'NPAL2', 'AP3D1',
                                   'CFA61', 'Unc', 'Unc', 'Unc', 'LRP2')) %>%
  filter(gene_name != 'Unc') %>%
  filter(gene_name %in% c('PTPRD', 'AP3D1', 'SECG', 'LRP2', 'CFA61'))

full_metadata <- read_csv('../../Combo_Tank_Survival/intermediate_files/clone_metadata.csv', 
                          show_col_types = FALSE) %>%
  group_by(clone_group) %>%
  filter(pct_missing == min(pct_missing)) %>% 
  ungroup %>%
  filter((species == 'Ac' & data_origin == 'vollmer') | species == 'Apa' | species == 'Apr') %>%
  filter((species == 'Ac' & pct_missing  < 0.3) | species == 'Apa' | species == 'Apr') %>%
  mutate(species_location = if_else(data_origin == 'vollmer', location, species))  %>%
  filter(data_origin == 'vollmer') %>%
  select(gen_id, disease_resistance) %>%
  rename(genotype = gen_id)


analysis_data <- full_data %>%
  inner_join(gene_names, by = 'gene_id') %>%
  left_join(full_metadata, by = 'genotype') %>%
  mutate(tank = str_extract(sequence_id, '[DH][0-9]')) %>%
  select(gene_name, gene_id, log2_cpm, genotype, 
         time, exposure, disease_resistance, tank) %>%
  filter(!is.na(disease_resistance))

#### Summary Stats ####
analysis_data %>%
  count(genotype, tank)

analysis_data %>%
  filter(gene_name == 'AP3D1') %>%
  count(genotype, time)

analysis_data %>%
  select(genotype, time, exposure, disease_resistance) %>%
  distinct %>%
  # group_by(genotype, disease_resistance, time) %>%
  # summarise()
  pivot_wider(names_from = time, 
              values_from = exposure,
              values_fn = ~str_c(., collapse = ' + '))

analysis_data %>%
  count(gene_name)

#### Analysis ####
gene_models <- analysis_data %>%
  filter(!genotype %in% c('U46')) %>% #, 'U37'
  nest_by(gene_name, gene_id) %>%
  mutate(model = list(lmer(log2_cpm ~ time * exposure * disease_resistance + (1 | genotype),
                           data = data)))


map(gene_models$model, ~anova(.x, ddf = "Kenward-Roger")) %>%
  set_names(gene_models$gene_name) %>%
  map(~as_tibble(.x, rownames = 'Term')) %>%
  bind_rows(.id = 'Gene Name') %>%
  mutate(across(where(is.double), ~round(., digits = 3)),
         across(where(is.double), ~sprintf("%.3f", .))) %>%
  mutate(df = str_c('F(', NumDF, ', ', DenDF, ') = ', `F value`),
         .keep = 'unused') %>%
  select(-`Mean Sq`) %>%
  rename(SS = `Sum Sq`,
         p = `Pr(>F)`) %>%
  relocate(df, .before = 'p') %>%
  mutate(Term = str_replace_all(Term, ':', ' x '),
         Term = str_replace_all(Term, '_', ' '),
         Term = str_to_title(Term)) %>% View
  write_csv('gene_expression_anova.csv')

analysis_data %>%
  ggplot(aes(x = disease_resistance, y = log2_cpm, colour = interaction(time, exposure))) +
  stat_summary_bin(bins = 3) +
  facet_wrap(~gene_name)


analysis_data %>%
  ggplot(aes(x = time, y = log2_cpm, shape = exposure, colour = disease_resistance)) +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~gene_name)


analysis_data %>%
  group_by(gene_name) %>%
  mutate(dr_cat = c('low', 'high')[(disease_resistance > median(disease_resistance)) + 1]) %>%
  ggplot(aes(x = time, y = log2_cpm, shape = exposure, colour = dr_cat)) +
  stat_summary(position = position_dodge(0.5)) +
  facet_wrap(~gene_name)

analysis_data %>%
  group_by(gene_name) %>%
  mutate(colour_group = case_when(disease_resistance >= quantile(disease_resistance, 0.75) ~ 'High',
                                disease_resistance <= quantile(disease_resistance, 0.25) ~ 'Low',
                                TRUE ~ 'Medium')) %>%
  mutate(colour_group = factor(colour_group, levels = c('Low', 'Medium', 'High'), 
                               ordered = TRUE)) %>%
  ungroup %>%
  mutate(exposure = factor(exposure, levels = c('H', 'D'))) %>%
  ggplot(aes(x = time, y = log2_cpm, colour = exposure, shape = colour_group)) +
  stat_summary(position = position_dodge(0.5)) +
  scale_colour_manual(values = c('D' = 'red', 'H' = 'blue'),
                      labels = c('D' = 'Disease', 'H' = 'Healthy')) +
  labs(colour = 'Disease Exposure',
       shape = 'Disease Resistance') +
  facet_wrap(~gene_name) +
  guides(shape = guide_legend(title.position = 'top', title.hjust = 0.5),
         colour = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  theme_classic() +
  theme(#aspect.ratio = 1,
    legend.direction = "horizontal", 
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = 'plain', colour = 'black', size = 20),
    legend.text = element_text(face = 'plain', colour = 'black', size = 20),
    axis.text = element_text(face = 'plain', colour = 'black', size = 10),
    # axis.title = element_text(face = 'plain', colour = 'black', size = 16),
    axis.title = element_blank(),
    panel.border = element_rect(colour = 'black', fill = 'transparent'),
    strip.text = element_text(face = 'plain', colour = 'black', size = 16),
    strip.placement = "outside", 
    # strip.text = element_blank(),
    strip.background = element_blank())


disease_resistance_quants <- c(quantile(full_metadata$disease_resistance, 0.025, na.rm = TRUE),
                               median(full_metadata$disease_resistance, na.rm = TRUE),
                               quantile(full_metadata$disease_resistance, 0.975, na.rm = TRUE))

gene_models %>%
  reframe(emmeans::emmeans(model, ~time * exposure * disease_resistance, 
                           at = list(disease_resistance = disease_resistance_quants)) %>%
            as_tibble()) %>%
  
  mutate(disease_resistance = case_when(disease_resistance == min(disease_resistance) ~ 'Low',
                                        disease_resistance == max(disease_resistance) ~ 'High',
                                        TRUE ~ 'Medium')) %>%
  mutate(disease_resistance = factor(disease_resistance, levels = c('Low', 'Medium', 'High'), 
                               ordered = FALSE)) %>%
  ggplot(aes(x = time, y = emmean, ymin = lower.CL, ymax = upper.CL, 
             colour = exposure, shape = disease_resistance)) +
  geom_pointrange(position = position_dodge(0.5)) +
  scale_colour_manual(values = c('D' = 'red', 'H' = 'blue'),
                      labels = c('D' = 'Disease', 'H' = 'Healthy')) +
  labs(colour = 'Disease Exposure',
       shape = 'Disease Resistance',
       x = NULL,
       y = expression(~log[2]~(CPM))) +
  facet_wrap(~gene_name) +
  guides(shape = guide_legend(title.position = 'top', title.hjust = 0.5),
         colour = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  theme_classic() +
  theme(#aspect.ratio = 1,
    legend.direction = "horizontal", 
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = 'plain', colour = 'black', size = 18),
    legend.text = element_text(face = 'plain', colour = 'black', size = 18),
    axis.text = element_text(face = 'plain', colour = 'black', size = 16),
    axis.title = element_text(face = 'plain', colour = 'black', size = 18),
    panel.border = element_rect(colour = 'black', fill = 'transparent'),
    strip.text = element_text(face = 'italic', colour = 'black', size = 18),
    strip.placement = "outside", 
    strip.background = element_blank())
# ggsave('../Results/gene_expression.png', height = 7, width = 7)

gene_models$model[[2]] %>%
  summary


cpm = 5.6 + t7*0.31 + H*0.24 + dr*0.75 + t7*H*-1.89 + t7*dr*-2.62 + H*dr*-1.27 + t7*H*dr*6.9


gene_models %>%
  reframe(emmeans(model, ~time * exposure) %>%
            contrast(list("T3" = c(1, 0, -1, 0),
                          'T7' = c(0, 1, 0, -1))) %>%
            broom::tidy(conf.int = TRUE)) %>%
  rename(time = contrast) %>%
  ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(label = if_else(p.value < 0.05, '*', '')),
            y = Inf, vjust = 1, colour = 'black',
            position = position_dodge(0.5),
            size = 6) +
  labs(x = NULL,
       y = expression(~log[2]~(CPM))) +
  facet_wrap(~gene_name) +
  theme_classic() +
  theme(#aspect.ratio = 1,
    legend.direction = "horizontal", 
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = 'plain', colour = 'black', size = 18),
    legend.text = element_text(face = 'plain', colour = 'black', size = 18),
    axis.text = element_text(face = 'plain', colour = 'black', size = 16),
    axis.title = element_text(face = 'plain', colour = 'black', size = 18),
    panel.border = element_rect(colour = 'black', fill = 'transparent'),
    strip.text = element_text(face = 'italic', colour = 'black', size = 18),
    strip.placement = "outside", 
    strip.background = element_blank())

gene_models %>%
  reframe(emmeans(model, ~time * exposure * disease_resistance, 
                  at = list(disease_resistance = disease_resistance_quants)) %>%
            contrast(list("T3_low" = c(1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                          'T7_low' = c(0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0),
                          
                          "T3_medium" = c(0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0),
                          'T7_medium' = c(0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0),
                          
                          "T3_high" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0),
                          'T7_high' = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1))) %>%
            broom::tidy(conf.int = TRUE) %>%
            as_tibble() %>%
            select(-term) %>%
            separate(contrast, into = c('time', 'resistance'))) %>%
  mutate(resistance = factor(resistance, levels = c('low', 'medium', 'high'))) %>%
  # filter(p.value < 0.05)
  
  ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high, 
             colour = resistance, shape = resistance)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(aes(label = if_else(p.value < 0.05, '*', '')),
            y = Inf, vjust = 1, colour = 'black',
            position = position_dodge(0.5),
            size = 6) +
  labs(colour = 'Disease Resistance',
       shape = 'Disease Resistance',
       x = NULL,
       y = expression(~log[2]~(CPM))) +
  facet_wrap(~gene_name) +
  guides(shape = guide_legend(title.position = 'top', title.hjust = 0.5),
         colour = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  theme_classic() +
  theme(#aspect.ratio = 1,
    legend.direction = "horizontal", 
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = 'plain', colour = 'black', size = 18),
    legend.text = element_text(face = 'plain', colour = 'black', size = 18),
    axis.text = element_text(face = 'plain', colour = 'black', size = 16),
    axis.title = element_text(face = 'plain', colour = 'black', size = 18),
    panel.border = element_rect(colour = 'black', fill = 'transparent'),
    strip.text = element_text(face = 'italic', colour = 'black', size = 18),
    strip.placement = "outside", 
    strip.background = element_blank())

library(emmeans)
emmeans::emmeans(gene_models$model[[2]], ~time * exposure * disease_resistance, 
                 at = list(disease_resistance = disease_resistance_quants)) %>% 
  contrast(list("T3_low" = c(1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                'T7_low' = c(0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0),
                
                "T3_medium" = c(0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0),
                'T7_medium' = c(0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0),
                
                "T3_high" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0),
                'T7_high' = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1))) %>%
  broom::tidy(conf.int = TRUE) %>%
  as_tibble() %>%
  select(-term) %>%
  separate(contrast, into = c('time', 'resistance'))


#### Make Model Plots ####
disease_resistance_levels <- modelr::seq_range(full_metadata$disease_resistance, 
                                               n = 100)

make_contrast <- function(timepoint, resistance_block, n_blocks){
  #timepoint is either T3 or T7
  #resistance_block is the order number of the resistance value
  #n_blocks is the number of total resistance values at
  #emmens must be T3_D, T7_D, T3_H, T7_H repeated for each block
  if(timepoint == 'T3'){
    contrast_block <- c(1, 0, -1, 0)
  } else if(timepoint == 'T7'){
    contrast_block <- c(0, 1, 0, -1)
  }
  
  out <- c(rep(rep(0, 4), resistance_block - 1), contrast_block, rep(rep(0, 4), n_blocks - resistance_block))
  # names(out) <- str_c(timepoint, '_R', resistance_block)
  list(out) %>%
    set_names(str_c(timepoint, '_R', resistance_block))
}

all_contrasts <- tibble(disease_resistance = disease_resistance_levels) %>%
  mutate(resistance_level = row_number()) %>%
  expand_grid(time = c('T3', 'T7')) %>%
  rowwise %>%
  mutate(contrast = make_contrast(time, resistance_level, length(disease_resistance_levels))) %>%
  ungroup %>%
  mutate(resistance_level = str_c('R', resistance_level))


tst <- gene_models %>%
  reframe(emmeans(model, ~time * exposure * disease_resistance, 
                  at = list(disease_resistance = disease_resistance_levels)) %>%
            contrast(all_contrasts$contrast) %>%
            broom::tidy(conf.int = TRUE) %>%
            as_tibble() %>%
            select(-term) %>%
            separate(contrast, into = c('time', 'resistance_level'))) %>%
  left_join(all_contrasts, by = c('time', 'resistance_level')) %>%
  nest_by(gene_name, time, .key = 'model_data') %>%
  left_join(gene_models %>%
              ungroup %>%
              select(gene_name, data) %>%
              unnest(data) %>%
              nest(data = -c(gene_name, time)) %>%
              rowwise %>%
              mutate(data = list(data %>%
                                   select(-tank) %>%
                                   pivot_wider(names_from = exposure, 
                                               values_from = 'log2_cpm') %>%
                                   mutate(fc = D - H) %>%
                                   filter(!is.na(fc)))) %>%
              ungroup,
            by = c('gene_name', 'time'))



tst_plot <- tst %>%
  select(-data) %>%
  unnest(model_data) %>%
  mutate(resistance_level = str_remove(resistance_level, 'R') %>%
           as.numeric()) %>%
  # filter(p.value < 0.05)
  
  ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high, 
             group = interaction(time, resistance_level))) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge(0.5)) +
  facet_wrap(~gene_name) 


plot_data <- as_tibble(ggplot_build(tst_plot)$data[[2]]) %>%
  left_join(tibble(PANEL = factor(1:4),
                   gene_name = unique(tst$gene_name)),
            by = 'PANEL') %>%
  select(-PANEL) %>%
  select(x, y, ymin, ymax, group, gene_name) %>%
  mutate(time = c('T3', 'T7')[(group + 1) %% 2 + 1]) %>%
  arrange(gene_name, time, x) %>%
  group_by(gene_name, time) %>%
  mutate(disease_resistance = disease_resistance_levels) %>%
  ungroup

(1:10 + 1) %% 2 + 1

plot_data %>%
  ggplot(aes(x = x, y = disease_resistance, colour = time)) +
  geom_point() +
  facet_wrap(~gene_name)

point_data <- plot_data %>%
  group_by(gene_name, time) %>%
  summarise(model = list(lm(x ~ disease_resistance)),
            .groups = 'drop') %>%
  left_join(tst %>%
              select(-model_data),
            by = c('gene_name', 'time')) %>%
  rowwise(gene_name, time) %>%
  reframe(mutate(data, x = predict(model, newdata = data))) %>%
  mutate(colour_group = case_when(disease_resistance >= quantile(full_metadata$disease_resistance, 0.75, na.rm = TRUE) ~ 'A',
                                  disease_resistance >= quantile(full_metadata$disease_resistance, 0.5, na.rm = TRUE) ~ 'B',
                                  disease_resistance <= quantile(full_metadata$disease_resistance, 0.25, na.rm = TRUE) ~ 'D',
                                  disease_resistance < quantile(full_metadata$disease_resistance, 0.5, na.rm = TRUE) ~ 'C',
                                  TRUE ~ 'error')) %>%
  mutate(colour_group = case_when(colour_group == 'A' ~ '(0,1]',
                                  colour_group == 'B' ~ '(1,2]',
                                  colour_group == 'error' ~ '(2,3]',
                                  colour_group == 'C' ~ '(3,4]',
                                  colour_group == 'D' ~ '(4,5]'),
         colour_group = factor(colour_group, levels = c('(0,1]', '(1,2]', '(2,3]', '(3,4]', '(4,5]')))


wesanderson::wes_palette("Zissou1", 9, type = "continuous")
colour_options <- set_names(c(wesanderson::wes_palette("Zissou1", 9, type = "continuous")[c(1,5,9)]),
                            c('A', 'B', 'C'))

wesanderson::wes_palette("Zissou1", 9, type = "continuous")
colour_options <- set_names(c(wesanderson::wes_palette("Zissou1", 9, type = "continuous")[c(1,3)],
                              'gray50', 
                              wesanderson::wes_palette("Zissou1", 9, type = "continuous")[c(5,9)]),
                            c('A', 'B', 'error', 'C', 'D'))

plot_data %>%
  ggplot(aes(x = x, y = y, ymin = ymin, ymax = ymax, group = group %% 2)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_ribbon(alpha = 0.1) +
  geom_line() +
  geom_point(data = point_data, inherit.aes = FALSE,
             aes(x = x, y = fc, colour = colour_group)) +
  facet_wrap(~gene_name) +
  scale_colour_manual(values = unname(colour_options), #breaks = c(0, 1, 2), 
                      labels = c('Highly\nResistant','' , '', 'Highly\nSusceptible'),
                      drop = FALSE) +
  scale_x_continuous(breaks = c(1, 2), minor_breaks = NULL, labels = c('Day 3', 'Day 7')) +
  
  labs(x = NULL,
       y = expression(~log[2]~("Fold Change")),
       # y = 'log2(CPM D) - log2(CPM H)',
       # y = 'log2 Fold Change',
       colour = NULL) +
  guides(colour = guide_colorsteps(reverse = TRUE)) +
  theme_classic() +
  theme(legend.title = element_text(face = 'plain', colour = 'black', size = 18),
        legend.text = element_text(face = 'plain', colour = 'black', size = 18),
        axis.text = element_text(face = 'plain', colour = 'black', size = 16),
        axis.title = element_text(face = 'plain', colour = 'black', size = 18),
        panel.border = element_rect(colour = 'black', fill = 'transparent'),
        strip.text = element_text(face = 'italic', colour = 'black', size = 18),
        strip.placement = "outside", 
        strip.background = element_blank())
ggsave('../Results/gene_expression.png', height = 7, width = 7)




#### Model FC ####
gene_models2 <- gene_models %>%
  mutate(data = list(pivot_wider(data, 
                                 names_from = exposure, 
                                 values_from = log2_cpm) %>%
                       mutate(fc = D - H) %>%
                       filter(time == 'T7'))) %>%
  mutate(model = list(lm(fc ~ disease_resistance,
                           data = data)))

map(gene_models2$model, ~summary(.x)) %>%
  set_names(gene_models$gene_name) %>%
  map(~as_tibble(.x, rownames = 'Term')) %>%
  bind_rows(.id = 'Gene Name') %>%
  mutate(across(where(is.double), ~round(., digits = 3)),
         across(where(is.double), ~sprintf("%.3f", .))) %>%
  mutate(df = str_c('F(', NumDF, ', ', DenDF, ') = ', `F value`),
         .keep = 'unused') %>%
  select(-`Mean Sq`) %>%
  rename(SS = `Sum Sq`,
         p = `Pr(>F)`) %>%
  relocate(df, .before = 'p') %>%
  mutate(Term = str_replace_all(Term, ':', ' x '),
         Term = str_replace_all(Term, '_', ' '),
         Term = str_to_title(Term)) %>% View


tibble(disease_resistance = modelr::seq_range())
gene_models2$model[[2]] %>%
  bind_cols()
  

gene_models2 %>%
  mutate(broom::glance(model))

tst <- gene_models2 %>%
  mutate(pred_data = list(tibble(disease_resistance = modelr::seq_range(data$disease_resistance, 10))),
         pred_data = list(bind_cols(pred_data, predict(model, newdata = pred_data, se.fit = TRUE))))

tst %>% 
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_ribbon(data = . %>% unnest(pred_data),
              aes(x = disease_resistance, ymin = fit - se.fit, ymax = fit + se.fit),
              alpha = 0.5) +
  geom_line(data = . %>% unnest(pred_data),
            aes(x = disease_resistance, y = fit)) +
  # geom_point(data = . %>% unnest(data),
  #            aes(x = disease_resitance, y = fc)) +
  facet_wrap(~gene_name)



tst$data[[2]]
