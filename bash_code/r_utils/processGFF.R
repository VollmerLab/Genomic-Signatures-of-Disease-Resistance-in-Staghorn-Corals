if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  structural_annotation_file <- args[1]
  genome_file <- args[2]
  NCORES <- as.integer(args[3])
} else {
  structural_annotation_file <- '/scratch/j.selwyn/science/Genomic-Signatures-of-Disease-Resistance-in-Staghorn-Corals/genome/genomic.gff'
  genome_file <- '/scratch/j.selwyn/science/Genomic-Signatures-of-Disease-Resistance-in-Staghorn-Corals/genome/acerv_genome.fasta'
  NCORES <- 24
}

suppressMessages(suppressWarnings(library(IRanges)))
suppressMessages(suppressWarnings(library(Biostrings)))
suppressMessages(suppressWarnings(library(bioseq)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(multidplyr)))


#### Preprocess Inputs ####
outdir <- dirname(structural_annotation_file)
base_name <- str_remove(structural_annotation_file, str_c(outdir, '/')) %>%
  str_remove('.gff')


#### Read in File ####
last_line <- system(str_c("grep -n '^##FASTA$' ", structural_annotation_file), intern = TRUE) %>%
  str_extract('[0-9]+') %>%
  as.integer()

structural_annotation_raw <- read_delim(structural_annotation_file, delim = '\t',
                                        skip = 7, #n_max = last_line - 2,
                                        col_types = 'ccciidcic', na = c('.', ''),
                                        col_names = c('Seqid', 'Source', 'Type', 'Start', 'End', 
                                                      'Score', 'Strand', 'Phase', 'Attributes')) %>%
  filter(!str_detect(Seqid, '^#'))
#Parsing issues are fine. Just marking between scaffolds. Remove with filter

genome <- readBStringSet(genome_file)
names(genome) <- str_remove(names(genome), ' Acropora cervicornis isolate K2 Acerv_scaffold_[0-9]+, whole genome shotgun sequence')

#### Identify Gene IDs ####
gene_locations <- structural_annotation_raw %>%
  filter(Type == 'gene') %>%
  mutate(gene_id = str_extract(Attributes, 'P5673_[0-9]+')) %>%
  select(Seqid, gene_id, Start, End, Strand) %>%
  arrange(gene_id) %>%
  full_join(structural_annotation_raw %>%
              filter(str_detect(Type, '[mt]RNA')) %>%
              select(Seqid, Start, End, Strand, Type),
            by = c('Seqid', 'Start', 'End', 'Strand')) %>%
  rename(gene_start = Start,
         gene_end = End,
         gene_type = Type) 
write_csv(gene_locations, str_c(outdir, '/', base_name, '_gene_id_locations.csv'))

#### Put Maker Structural Elements within each Gene ####
maker_annotatons <- structural_annotation_raw %>%
  filter(!Type %in% c('contig', 'gene', 'mRNA', 'tRNA', 'region')) %>%
  mutate(gene_id = str_extract(Attributes, 'P5673_[0-9]+')) %>%
  select(-Source) %>%
  nest(maker_structural_elements = c(Type, Start, End, Score, Phase, Attributes)) %>%
  full_join(gene_locations,
            by = c('Seqid', 'gene_id', 'Strand')) %>%
  select(Seqid, gene_id, gene_start, gene_end, gene_type, Strand, 
         maker_structural_elements) %>%
  unnest(maker_structural_elements) %>%
  select(-Score) %>%
  rename(structure_type = Type,
         structure_start = Start,
         structure_end = End,
         structure_phase = Phase,
         structure_attributes = Attributes)
write_csv(maker_annotatons, str_c(outdir, '/', base_name, '_structural_annotation.csv'))


#### Identify Codons within each gene ####
cluster <- new_cluster(NCORES)
cluster_library(cluster, c('dplyr', 'bioseq', 'Biostrings'))
cluster_copy(cluster, c('genome'))

sort_strands <- function(data){
  strand <- unique(data$Strand)
  if(strand == '-'){
    out <- arrange(data, Seqid, -position)
  } else {
    out <- arrange(data, Seqid, position)
  }
  out
}

codons_positions <- maker_annotatons %>%
  filter(structure_type == 'CDS') %>%
  rowwise %>%
  partition(cluster) %>%
  mutate(sequence = subseq(genome[[Seqid]], structure_start, structure_end) %>% 
           as.character %>% dna) %>%
  collect %>%
  arrange(Seqid, structure_start) %>%
  
  mutate(cds_id = str_c(gene_id, structure_start, structure_end, sep = '_')) %>%
  select(Seqid, gene_id, cds_id, structure_start, structure_end, Strand, 
         structure_phase, sequence) %>%
  rowwise(everything(), -sequence) %>%
  partition(cluster) %>%
  summarise(tibble(position = seq(structure_start, structure_end, by = 1),
                   dna_base = seq_split_pattern(sequence, '') %>% unlist %>% dna)) %>%
  collect() %>%
  ungroup %>%
  group_split(Strand) %>%
  map(sort_strands) %>%
  bind_rows %>%
  group_by(gene_id) %>%
  mutate(codon = (row_number() - 1) %/% 3 + 1,
         codon_position = (row_number() - 1) %% 3 + 1) %>%
  
  group_by(Seqid, gene_id, codon) %>%
  partition(cluster) %>%
  mutate(codon_seq = seq_combine(dna_base, collapse = '')) %>%
  collect %>%
  ungroup %>%
  select(Seqid:cds_id, Strand, structure_phase, position, dna_base, 
         codon, codon_position, codon_seq)

# codons_positions <- tmp3
write_csv(codons_positions, str_c(outdir, '/', base_name, '_codon_positions.csv.gz'))
