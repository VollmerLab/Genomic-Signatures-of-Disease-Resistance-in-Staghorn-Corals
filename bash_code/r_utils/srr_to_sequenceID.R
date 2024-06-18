

suppressMessages(suppressWarnings(library(tidyverse)))

rename_codes <- read_delim('Data/bam.list', 
                           delim = '\t', 
                           col_names = 'file',
                           show_col_types = FALSE) %>%
  mutate(ID = str_remove(file, dirname(file)) %>%
           str_remove('/') %>% 
           str_remove('-fp.*bam$')) %>%
  select(-file) %>%
  mutate(SRR_id = str_extract(ID, 'SRR[0-9]+')) %>%
  filter(SRR_id %in% read_lines('SRR.numbers_baum') | 
           is.na(SRR_id)) %>%
  left_join(read_csv('SRR.numbers', 
                     show_col_types = FALSE) %>%
              rename(ID = LibraryName,
                     SRR_id = Run),
            by = 'ID') %>%
  mutate(SRR_id = coalesce(SRR_id.x, SRR_id.y), 
         .keep = 'unused')

new_file_names <- rename_codes %>%
  rowwise(ID, SRR_id) %>%
  summarise(file = list.files(path = 'variant_calling/raw_reads', 
                            pattern = SRR_id, 
                            full.names = TRUE),
            .groups = 'drop') %>%
  mutate(new_name = str_replace(file, SRR_id, ID)) %>%
  select(file, new_name)


done <- with(new_file_names, file.rename(file, new_name))
message('All files renamed from SRR to metadata name')