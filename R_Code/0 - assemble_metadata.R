#### Script to Initially Assemble metadata for Project ####

library(tidyverse)
library(magrittr)
library(xml2)
library(adegenet)
library(poppr)
library(reticulate)

#### Functions ####
convert_to_decimal <- function(x){
  out <- str_split(x, 'W|N') %>%
    unlist(recursive = FALSE) %>%
    as.numeric() %>%
    magrittr::divide_by(c(1, 60)) %>%
    sum
  
  if_else(str_detect(x, 'W'), -1 * out, out)
}

srr_to_srs <- function(data){
  pysradb <- import('pysradb')
  db = pysradb$SRAweb()
  
  srr_ids <- filter(data, !is.na(SRR_id)) %>%
    pull(SRR_id)
  
  out <- db$srr_to_srs(srr_ids, detailed = FALSE, 
                       sample_attribute = FALSE, 
                       expand_sample_attributes = FALSE) %>%
    as_tibble() %>%
    select(run_accession, sample_accession) %>% 
    rename(SRS_id = sample_accession, SRR_id = run_accession)
  full_join(data, out, by = 'SRR_id')
}

#### Gather Collection Location Data ####
florida_locations <- read_csv('../Data/SP.pop.data.csv', show_col_types = FALSE) %>%
  select(Genotype, subregion, reefName, x.coord, y.coord) %>%
  rename(gen_id = Genotype,
         lat = x.coord,
         lon = y.coord) %>%
  mutate(reef = str_c(subregion, reefName, sep = '-'), .keep = 'unused')

panama_locations <- tribble(
  ~'reef', ~'lon_garmin', ~'lat_garmin', 
  'tetas', '82W6.068', '9N16.579',
  'sebastians', '82W7.631', '9N15.274',
  'holy shit', '82W6.929', '9N16.773',
  'CK14', '82W7.557', '9N15.238',
  'CK4', '82W7.625', '9N15.517'
) %>%
  rowwise %>% 
  mutate(lat = convert_to_decimal(lat_garmin),
         lon = convert_to_decimal(lon_garmin)) %>%
  ungroup %>%
  select(-ends_with('garmin')) %>%
  mutate(reef = str_replace_all(reef, c('tetas' = 'Tet', 'sebastians' = 'SR', 'holy shit' = 'HS')))

#### Process Baum Lab Locations ###
baum_locations <- read_xml('../Data/biosample_result.xml') %>%
  as_list() %>%
  as_tibble %>%
  unnest_wider(BioSampleSet) %>%
  rowwise %>%
  mutate(SRS_id = Ids[[3]][[1]]) %>%
  mutate(location = Attributes[[3]][[1]] %>% str_remove(' \\(.*\\)')) %>%
  select(SRS_id, location)


#### Read in Genetic file data & Join All data ####
all_individual_data <- read_delim('../Data/bam.list', delim = '\t', 
                                  col_names = 'file', show_col_types = FALSE) %>%
  mutate(ID = str_remove(file, dirname(file)) %>% str_remove('/') %>% str_remove('-fp.*bam$')) %>%
  select(-file) %>%
  mutate(SRR_id = str_extract(ID, 'SRR[0-9]+')) %>% 
  srr_to_srs() %>%
  left_join(baum_locations, by = 'SRS_id') %>%
  select(-SRR_id, -SRS_id) %>%
  mutate(gen_id = str_extract(ID, '[A-Za-z0-9]+$')) %>%
  left_join(florida_locations, by = 'gen_id') %>%
  mutate(data_origin = if_else(str_detect(ID, 'baum'), 'baum', 'vollmer'),
         species = str_extract(ID, 'A[cp]|ac|apr|apa') %>% str_to_sentence(),
         location = if_else(is.na(location), str_extract(ID, 'FL|PA'), location) %>% str_replace_all(c('^FL$' = 'Florida', '^PA$' = 'Panama')),
         reef = if_else(data_origin == 'vollmer' & location == 'Panama', str_extract(ID, 'CK4|CK14|HS|SR|Tet'), reef)) %>%
  select(ID, gen_id, data_origin, species, location, reef, lat, lon) %>%
  left_join(panama_locations, by = 'reef') %>%
  mutate(lat = coalesce(lat.x, lat.y),
         lon = coalesce(lon.x, lon.y),
         .keep = 'unused')

dir.create('../intermediate_files')
write_csv(all_individual_data, '../intermediate_files/collected_sample_metadata.csv')
