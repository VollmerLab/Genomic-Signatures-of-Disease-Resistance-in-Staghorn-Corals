#### Libraries ####
library(tidyverse)
library(magrittr)
library(adegenet)
library(poppr)
library(SNPRelate)

#### Settings - and order of filtering ####
confidence_cutoff <- 0.99 #require probability of genotype to be greater than this
#Identify Clone Groups 
#Remove Clones based on which has least missing loci
individual_missingness <- 0.3 #Maximum % of loci missing in an individual to keep individual
#Split into two datasets one to compare A. cerv and A. palm & one to compare within A. cerv
#Remove if entirely missing within one species/location
#Remove if monomorphic
locus_missingness <- 0.1 #Maximum % of individuals missing in a locus to keep locus
minor_allele_frequency <- 0.05
LD_cutoff <- 0.5 #Linkage disequilibrium cutoff 

#### Functions ####
process_angsd <- function(genotype_file, probability_file, sample_names, cutoff){
  filter_below <- function(genotype, probabilites, cutoff){
    to_swap_to_na <- probabilites < cutoff
    probabilites[to_swap_to_na] <- NA
    
    message('Swapped ', scales::percent(sum(to_swap_to_na, na.rm = TRUE) / prod(dim(probabilites)), accuracy = 1), 
            ' genotypes to NA due to having probability less than ', cutoff)
    
    genotype * round(probabilites, 0)
  }
  
  genotypes <- read_csv(genotype_file, 
                        col_types = cols(.default = col_integer(),
                                         contig = col_character(),
                                         ref = col_character(),
                                         alt = col_character())) %>%
    mutate(locus_number = row_number(),
           genotype = cbind(!!!syms(sample_names)),
           .keep = 'unused')
  
  
  
  if(str_detect(probability_file, '\\*')){
    genotype_probabilities <- list.files(path = dirname(probability_file), 
               pattern = str_remove(probability_file, str_c(dirname(probability_file), '/')),
               full.names = TRUE) %>%
      map(~read_csv(.x, 
                    col_types = cols(.default = col_double(),
                                     contig = col_character(),
                                     ref = col_character(),
                                     alt = col_character()))) %>%
      bind_rows() %>%
      mutate(locus_number = row_number(),
             genotype_probability = cbind(!!!syms(sample_names)),
             .keep = 'unused')
  } else {
    genotype_probabilities <- read_csv(probability_file, 
                                       col_types = cols(.default = col_double(),
                                                        contig = col_character(),
                                                        ref = col_character(),
                                                        alt = col_character())) %>%
      mutate(locus_number = row_number(),
             genotype_probability = cbind(!!!syms(sample_names)),
             .keep = 'unused')
  }
  
  
  full_join(genotypes, genotype_probabilities,
            by = c('contig', 'position', 'ref', 'alt', 'locus_number')) %>%
    
    #Filter out genotypes below confidence threshold
    mutate(genotype = filter_below(genotype = genotype, 
                                   probabilites = genotype_probability, 
                                   cutoff = cutoff), 
           .keep = 'unused') %>%
    
    #Remove Loci fixed after confidence cutoff
    mutate(alt_freq = rowSums(genotype, na.rm = TRUE) / (2 * (rowSums(!is.na(genotype)))),
           ref_freq = 1 - alt_freq, 
           .before = genotype) %>%
    filter(ref_freq != 1, alt_freq != 1) %>%
    arrange(locus_number) %>%
    select(-locus_number, -ref_freq, -alt_freq) %>%
    
    #Convert to genlight to more easily use in adegenet etc.
    mutate(snp_loc = str_c(ref, alt, sep = '/'),
           loc_name = str_c(contig, position, sep = '-'),
           .before = genotype) %$%
    
    new('genlight', gen = t(genotype), ploidy = 2,
        loc.all = snp_loc, loc.names = loc_name, 
        chromosome = contig, position = position)
}

flip_loci <- function(gl){
  #Use to flip loci that the minor allele was called as the major in acerv due to calls being done for both species
  tmp <- as_tibble(strata(gl)) %>%
    mutate(genotype = as.matrix(gl))
  
  just_acerv <- tmp %>%
    filter(data_origin == 'vollmer',
           species == 'Ac')
  
  
  
  n_alleles <- 2 * colSums(!is.na(just_acerv$genotype))
  # p <- (2 * colSums(just_acerv$genotype == 0, na.rm = TRUE) + colSums(just_acerv$genotype == 1, na.rm = TRUE)) / n_alleles
  q <- (2 * colSums(just_acerv$genotype == 2, na.rm = TRUE) + colSums(just_acerv$genotype == 1, na.rm = TRUE)) / n_alleles
  loci_flip <- names(q)[q > 0.5]
  
  
  to_flip <- tmp$genotype[,loci_flip]
  to_flip[to_flip == 2] <- 3
  to_flip[to_flip == 0] <- 2
  to_flip[to_flip == 3] <- 0
  
  tmp$genotype[,loci_flip] <- to_flip
  
  out <- new('genlight', gen = tmp$genotype, ploidy = 2,
             loc.all = gl@loc.all, loc.names = gl@loc.names, 
             chromosome = gl@chromosome, position = gl@position)
  strata(out) <- strata(gl)
  out
}

remove_monomorphic <- function(snp_clone){
  snp_mat <- as.matrix(snp_clone)
  #n_ind <- dim(snp_mat)[1]
  
  all_0 <- colSums(snp_mat == 0, na.rm = TRUE) / colSums(!is.na(snp_mat))
  # all_2 <- colSums(snp_clone == 2, na.rm = TRUE) / n_ind
  loci_to_remove <- unname(which(all_0 == 1 | all_0 == 0))
  message(scales::comma(length(loci_to_remove)), ' monomorphic loci removed')
  snp_clone[,-loci_to_remove]
}

remove_missing_pop <- function(snp_clone, pop){
  setPop(snp_clone) <- pop
  
  missing_loci_in_pops <- seppop(snp_clone) %>%
    map(~as.matrix(.x) %>%
          is.na(.) %>%
          colSums(.) %>%
          equals(nInd(.x)) %>%
          which()) %>%
    unlist %>%
    names %>%
    unique %>%
    str_remove(str_c(levels(pop(snp_clone)), collapse = '|')) %>%
    str_remove('^\\.') %>%
    unique
  
  message('Removed ', scales::comma(length(missing_loci_in_pops)), 
          ' loci which were all NAs in one of the populations')
  
  snp_clone[,!snp_clone@loc.names %in% missing_loci_in_pops]
}

genlight_to_gds <- function(gen_light, out_name){
  SNPRelate::snpgdsCreateGeno(gds.fn = out_name, 
                              genmat = as.matrix(gen_light), 
                              sample.id = gen_light@ind.names, 
                              snp.id = gen_light@loc.names, 
                              snp.chromosome = str_extract(gen_light@chromosome, '[0-9]+') %>% as.integer,  
                              snp.position = gen_light@position, 
                              snp.allele = gen_light@loc.all, 
                              snpfirstdim = FALSE)
  
  genofile <- SNPRelate::snpgdsOpen(out_name)
  SNPRelate::snpgdsSummary(genofile)
  SNPRelate::snpgdsClose(genofile)
  out_name
}

get_gds <- function(gds, column){
  z <- 0
  if(is.character(gds)){
    gds <- snpgdsOpen(gds)
    z <- 1
  }
  out <- gdsfmt::index.gdsn(gds, index = column) %>%
    gdsfmt::read.gdsn()
  
  if(z == 1){
    snpgdsClose(gds)
  }
  out
}

choose_snps <- function(genomic_data, samples_use = NULL, ...){
  genofile <- snpgdsOpen(genomic_data)
  if(is.null(samples_use)){
    samples_use <- get_gds(genofile, 'sample.id')
  }
  
  snps_used <- snpgdsSelectSNP(genofile, autosome.only = FALSE,
                               sample.id = samples_use, ...)
  snpgdsClose(genofile)
  snps_used
}

find_unlinked_loci <- function(gds_file, ...){
  genofile <- snpgdsOpen(gds_file)
  
  prune_snps <- snpgdsLDpruning(genofile, autosome.only = FALSE, ...)
  
  snpgdsClose(genofile)
  
  unlist(prune_snps) %>%
    unname
}

preprocess_loci <- function(genomic, samples_use, population, save_dir, name_base, linkage_filter){
  gc(verbose = FALSE)
  
  #Subset loci & remove monomorphic/entirely missing loci
  subsample_genomic <- genomic[genomic@ind.names %in% samples_use] %>%
    remove_monomorphic() %>%
    remove_missing_pop(population)
  subsample_genomic_gds <- genlight_to_gds(subsample_genomic, str_c(save_dir, '/initial_', name_base, '.gds'))
  
  #Filter loci on MAF/missingness threshold
  locus_keep <- choose_snps(subsample_genomic_gds, 
                            maf = minor_allele_frequency, 
                            missing.rate = locus_missingness)
  filtered_subsample <- subsample_genomic[,locus_keep]
  
  #Remove linked loci
  if(linkage_filter){
    unlinked_filtered_subsample_loci <- find_unlinked_loci(subsample_genomic_gds, ld.threshold = LD_cutoff, 
                                                           method = 'composite', start.pos = 'random', 
                                                           snp.id = locus_keep, slide.max.bp = 250000)
    unlinked_filtered_subsample <- filtered_subsample[,unlinked_filtered_subsample_loci]
    write_rds(unlinked_filtered_subsample, str_c(save_dir, '/unlinked_preprocessed_', name_base, '.rds'), compress = 'xz')
    filtered_subsample_gds <- genlight_to_gds(unlinked_filtered_subsample, str_c(save_dir, '/unlinked_preprocessed_', name_base, '.gds'))
    return(unlinked_filtered_subsample)
    
  } else {
    write_rds(filtered_subsample, str_c(save_dir, '/preprocessed_', name_base, '.rds'), compress = 'xz')
    filtered_subsample_gds <- genlight_to_gds(filtered_subsample, str_c(save_dir, '/preprocessed_', name_base, '.gds'))
    return(filtered_subsample)
  }
}

#### Read in data ####
sample_data <- read_csv('../intermediate_files/collected_sample_metadata.csv', 
                        show_col_types = FALSE)

resistance_data <- read_csv('../intermediate_files/disease_resistance.csv', show_col_types = FALSE)

#Use if skipping SNP calling step
# genomic_data <- process_angsd(genotype_file = '../Data/genotypes.csv.xz',
#                               probability_file = '../Data/genotype_probabilities_.*.csv.xz',
#                               sample_names = sample_data$ID, 
#                               cutoff = confidence_cutoff)

#Use if calling SNPs on your own 
genomic_data <- process_angsd(genotype_file = '../variant_calling/genotyping/genotypes.csv',
                              probability_file = '../variant_calling/genotyping/genotype_probabilities.csv',
                              sample_names = sample_data$ID,
                              cutoff = confidence_cutoff)
strata(genomic_data) <- sample_data
genomic_data <- flip_loci(genomic_data)

write_rds(genomic_data, '../intermediate_files/initial_full_genomic.rds')

#### Identify Clones ####
png('../Results/allSP_clone_removal.png', width = 7, height = 7, units = 'in', res = 125)
clone_data <- as.snpclone(genomic_data)
filter_stat_data <- filter_stats(clone_data, distance = 'bitwise.dist', plot = TRUE)
cutoff_clones <- cutoff_predictor(filter_stat_data$farthest$THRESHOLDS) 
abline(v = cutoff_clones, col = 'red') 
dev.off()

mlg.filter(clone_data, distance = 'bitwise.dist', algorithm = "farthest") <- cutoff_clones
# mlg.table(clone_data)

updated_metadata <- mutate(sample_data, 
                           clone_group = str_c('clone', 
                                               mll(clone_data, type = 'contracted'), 
                                               sep = '_'))

#### Initial Missing Data per Individual ####
# choose clone to keep least missing data initially
metadata_clones_removed <- as.matrix(genomic_data) %>%
  is.na %>%
  rowSums() %>%
  divide_by(nLoc(genomic_data)) %>%
  enframe(name = 'ID', value = 'pct_missing') %>%
  left_join(updated_metadata, ., by = 'ID') %>%
  left_join(resistance_data %>%
              filter(treatment == 'D') %>%
              select(-location, -treatment),
            by = 'gen_id') %T>%
  write_csv('../intermediate_files/clone_metadata.csv') %>%
  
  group_by(clone_group) %>%
  filter(pct_missing == min(pct_missing)) %>%
  ungroup %>%
  filter(pct_missing < individual_missingness) %>%
  filter((data_origin == 'vollmer' & !is.na(disease_resistance)) | 
           data_origin == 'baum')

write_csv(select(metadata_clones_removed, -pct_missing), '../intermediate_files/preprocessed_metadata.csv')

#### PreProcess Loci ####
preprocessed_loci <- expand_grid(pop_group = c('apalm_acerv', 'acerv'),
                                 only_unlinked = c(TRUE, FALSE)) %>%
  mutate(pop_form = case_when(pop_group == 'apalm_acerv' ~ list(~species),
                              TRUE ~ list(~location)),
         samples = case_when(pop_group == 'apalm_acerv' ~ list(filter(metadata_clones_removed, species != 'Apr') %>% pull(ID)),
                             pop_group == 'acerv' ~ list(filter(metadata_clones_removed, data_origin == 'vollmer') %>% pull(ID)),
                             pop_group == 'acerv_panama' ~ list(filter(metadata_clones_removed, data_origin == 'vollmer', location == 'Panama') %>% pull(ID)),
                             pop_group == 'acerv_florida' ~ list(filter(metadata_clones_removed, data_origin == 'vollmer', location == 'Florida') %>% pull(ID)),
                             TRUE ~ list(NA_character_))) %>%
  rowwise %>%
  mutate(subset_population = list(preprocess_loci(genomic = genomic_data,
                                                  samples_use = samples,
                                                  population = pop_form,
                                                  save_dir = '../intermediate_files',
                                                  name_base = pop_group,
                                                  linkage_filter = only_unlinked)))

preprocessed_loci$subset_population[[4]]
