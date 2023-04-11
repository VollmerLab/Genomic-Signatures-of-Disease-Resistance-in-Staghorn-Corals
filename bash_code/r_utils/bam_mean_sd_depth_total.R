
if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  the_dir <- args[1]
  out_file <- args[2]
  sequence_length <- as.integer(args[3])
} else {
  the_dir <- '/scratch/j.selwyn/variant_calling/flyeHQ_25-July-2022/genotyping'
  out_file <- '/scratch/j.selwyn/variant_calling/flyeHQ_25-July-2022/genotyping/depth.stats'
  sequence_length <- 140
}

suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tibble)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(forcats)))
suppressMessages(suppressWarnings(library(ggplot2)))

options(bitmapType="cairo")
out_dir <- str_remove(out_file, '/[a-z\\.]+$')

individual_depth_stats <- list.files(the_dir, pattern = 'depths$') %>%
  tibble(ind = .) %>%
  mutate(file = str_c(the_dir, ind, sep = '/'),
         ind = str_remove(ind, '.depths'),
         ind = str_remove_all(ind, 'HVFGNDSX3_s1_UDP[0-9]+_i7_8bp-UDP[0-9]+_i5_8bp2_7067-AT-0001_'),
         ind = str_remove_all(ind, '-fp-fqscrn-repr')) %>%
  rowwise(ind) %>%
  summarise(read_delim(file, delim = '\t', col_names = c('chrom', 'length', 'mapped', 'unmapped'),
                       col_types = 'cnnn'),
            .groups = 'drop') %>%
  mutate(depth = (mapped * sequence_length) / length) %>%
  filter(chrom != '*') %>%
  select(ind, chrom, depth) %T>%
  write_csv(str_c(out_dir, 'individual_contig_depths.csv', sep = '/')) 

#Output overall stats for angsd
individual_depth_stats %>%
  group_by(chrom) %>% 
  summarise(depth = sum(depth),
            .groups = 'drop') %>%
  summarise(mean_depth = mean(depth),
            sd_depth = sd(depth)) %>%
  pivot_longer(cols = everything()) %>%
  pull(value) %>%
  write_lines(out_file)

dir.create(str_c(out_dir, 'depth_histograms', sep = '/'))
depth_histograms <- individual_depth_stats %>%
  nest_by(ind) %>%
  mutate(mean_depth = mean(data$depth, na.rm = TRUE),
         sd_depth = sd(data$depth, na.rm = TRUE),
         se_depth = sd_depth / sqrt(nrow(data)),
         depth_lwr = mean_depth - 1.96 * se_depth,
         depth_upr = mean_depth + 1.96 * se_depth) %>%
  summarise(plots = list(ggplot(data, aes(x = depth)) +
                           geom_histogram(bins = 50) +
                           geom_vline(xintercept = mean_depth - depth_lwr, linetype = 'dashed') +
                           geom_vline(xintercept = mean_depth + depth_upr, linetype = 'dashed') +
                           geom_vline(xintercept = mean_depth) +
                           scale_x_continuous(trans = scales::pseudo_log_trans(base = 10), labels = scales::comma_format()) +
                           labs(title = ind,
                                x = 'log10(Contig Depth)',
                                y = 'Count') +
                           theme_classic()),
            .groups = 'drop') %>%
  mutate(out_file = str_c(out_dir, '/depth_histograms/', ind, '_depth_hist.png')) %>%
  rowwise %>%
  mutate(save = ggsave(out_file, plot = plots, width = 15, height = 15))


individual_depth_comparison <- individual_depth_stats %>%
  group_by(ind) %>%
  summarise(mean_depth = mean(depth),
            sd_depth = sd(depth),
            n = n(),
            se_depth = sd_depth / sqrt(n)) %>%
  mutate(depth_lwr = mean_depth - 1.96 * se_depth,
         depth_upr = mean_depth + 1.96 * se_depth,
         ind = fct_reorder(ind, mean_depth)) %>%
  ggplot(aes(y = ind, x = mean_depth, xmin = depth_lwr, xmax = depth_upr)) +
  geom_vline(data = . %>% summarise(se_depth = sd(mean_depth) / sqrt(n()),
                                    mean_depth = mean(mean_depth)), 
             aes(xintercept = mean_depth - 1.96 * se_depth), linetype = 'dashed') +
  geom_vline(data = . %>% summarise(se_depth = sd(mean_depth) / sqrt(n()),
                                    mean_depth = mean(mean_depth)), 
             aes(xintercept = mean_depth + 1.96 * se_depth), linetype = 'dashed') +
  geom_vline(data = . %>% summarise(mean_depth = mean(mean_depth)), 
             aes(xintercept = mean_depth)) +
  geom_pointrange() +
  labs(x = 'Mean Contig Depth',
       y = NULL) +
  theme_classic()
ggsave(str_c(out_dir, '/individual_depthComparison.png'), plot = individual_depth_comparison, width = 15, height = 15)
