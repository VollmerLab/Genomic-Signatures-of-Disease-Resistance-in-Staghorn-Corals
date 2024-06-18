if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  linkage_dir <- args[1]
  lfmm_results_file <- args[2]
} else {
  linkage_dir <- '/scratch/j.selwyn/science/Genomic-Signatures-of-Disease-Resistance-in-Staghorn-Corals/variant_calling/linkage'
  lfmm_results_file <- '/scratch/j.selwyn/science/Genomic-Signatures-of-Disease-Resistance-in-Staghorn-Corals/Results/lfmm_results.csv'
}

setwd(linkage_dir)

library(tidyverse)

#### Functions ####
read_in_linkage <- function(x, pos){
  filter(x, r2 >= 0.5)
}

clumpify_snps <- function(lfmm, linkage, min_r2, max_dist, seed_p){
  linkage <- filter(linkage, r2 >= min_r2, dist <= max_dist)
  
  store <- dplyr::select(lfmm, locus, fdr) %>%
    arrange(fdr) %>%
    mutate(clump_id = numeric(nrow(.)),
           clump_loci = vector('list', nrow(.))) %>%
    filter(fdr <= seed_p)
  
  i <- 0
  while(i != nrow(store) & i <= nrow(lfmm)){
    i <- i+1
    clumped_loci <- filter(linkage, locus1 == store$locus[i] | locus2 == store$locus[i]) %>%
      dplyr::select(locus1, locus2) %>%
      as.matrix() %>%
      as.character() %>% unique %>%
      str_subset(store$locus[i], negate = TRUE)
    
    store2 <- filter(store, !locus %in% clumped_loci)
    store2$clump_id[i] <- i
    store2$clump_loci[[i]] <- clumped_loci
    
    linkage <- filter(linkage, 
                      !locus1 %in% clumped_loci, !locus2 %in% clumped_loci, 
                      locus1 != store$locus[i], locus2 != store$locus[i])
    
    if(interactive()){
      message('Round ', i, ': clumped ', scales::comma(length(clumped_loci)), ' loci including ', 
              scales::comma(nrow(store) - nrow(store2)), ' potential seed loci. ', scales::comma(nrow(store2) - i), ' remaining')
    }
    store <- store2
  }
  
  
  dplyr::select(store, locus, clump_id, clump_loci) %>%
    unnest(clump_loci, keep_empty = TRUE) %>%
    pivot_longer(cols = c(locus, clump_loci),
                 names_to = 'is_index',
                 values_to = 'locus') %>%
    distinct %>%
    mutate(is_index = is_index == 'locus') %>%
    filter(!is.na(locus)) %>%
    mutate(clump_id = str_c('clump', clump_id, sep = '_'))
}
# locus_data <- lfmm_clumping$locus_data; ld_file <- lfmm_clumping$ld_file
clump_snps <- function(locus_data, ld_file, file_size, chunkSize){
  
  if(file_size < 500000000){
    linkage_disequilibrium <- read_delim(ld_file, delim = '\t', skip = 1,
                                         col_names = c('locus1', 'locus2', 'dist', 'r2_expG', 
                                                       'D', 'D_prime', 'r2', 'sample_size',
                                                       'maf1', 'maf2', 'hap00', 'hap01', 
                                                       'hap10', 'hap11', 'hap_maf1', 
                                                       'hap_maf2', 'chi2', 'loglike', 'nIter'),
                                         col_types = cols(.default = col_double(),
                                                          locus1 = col_character(),
                                                          locus2 = col_character())) %>%
      filter(r2 >= 0.5) %>%
      mutate(across(starts_with('locus'), ~str_replace(., ':', '-')))
  } else {
    linkage_disequilibrium <- read_delim_chunked(ld_file, delim = '\t', skip = 1,
                                                 col_names = c('locus1', 'locus2', 'dist', 'r2_expG', 
                                                               'D', 'D_prime', 'r2', 'sample_size',
                                                               'maf1', 'maf2', 'hap00', 'hap01', 
                                                               'hap10', 'hap11', 'hap_maf1', 
                                                               'hap_maf2', 'chi2', 'loglike', 'nIter'),
                                                 col_types = cols(.default = col_double(),
                                                                  locus1 = col_character(),
                                                                  locus2 = col_character()),
                                                 callback = DataFrameCallback$new(read_in_linkage),
                                                 chunk_size = chunkSize) %>%
      filter(r2 >= 0.5) %>%
      mutate(across(starts_with('locus'), ~str_replace(., ':', '-')))
  }
  
  gc()
  
  out <- clumpify_snps(lfmm = locus_data, linkage = linkage_disequilibrium, 
                min_r2 = 0.5, max_dist = 250000, seed_p = 1)
  gc()
  out
}

#### Get Linkage Results Files ####
lfmm_results <- read_csv(lfmm_results_file, 
                         show_col_types = FALSE) %>%
  mutate(locus = str_c(chromosome, position, sep = '-')) %>%
  nest(locus_data = -c(chromosome)) %>%
  mutate(ld_file = str_c(chromosome, '.ld'),
         exists = file.exists(ld_file)) %>%
  
  filter(exists)

# (n(n-1))/2

lfmm_clumping <- lfmm_results %>%
  # sample_n(10) %>%
  # slice(1) %>%
  # filter(ld_file == 'Acerv_scaffold_108.ld.gz') %>%
  rowwise(chromosome) %>%
  mutate(size = file.info(ld_file)$size,
         n_pairs = (nrow(locus_data) * (nrow(locus_data) - 1)) / 2) %>%
  # filter(size > 1000000000) %>% arrange(size)
  summarise(clump_snps(locus_data, ld_file, file_size = size, chunkSize = 5e8),
            .groups = 'drop')

write_csv(lfmm_clumping, 'lfmm_snp_clumps.csv.gz')