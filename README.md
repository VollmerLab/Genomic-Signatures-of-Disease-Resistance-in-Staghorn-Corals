# Genomic-Signatures-of-Disease-Resistance-in-Endangered-Staghorn-Corals-
Code associated with manuscript analyzing the genomic basis of disease resistance in staghorn corals

Data needed to analyze this data from scratch:
  - the assembled genome (NCBI BioProj: PRJNA948411)
  - the *A. cervicornis* WGS reads to be used (NCBI BioProj: PRJNA950067)
  - Supplemental WGS for *A. palmata* and *A. prolifera* (NCBI BioProj: PRJNA473816)

Can skip SNP preprocessing and calling by starting with the already identified SNPs `Data/genotypes.csv.gz` & `Data/genotype_probabilities.csv.gz`. If so skip preprocessing and SNP calling steps.

## 1. Preprocess Sequences
  - all (bash scripts)[bash_code] are designed to run on a SLURM based HPC
  - scripts should be modified such that the `${scriptDir}` variable is assigned at the beginning of scripts to the `bash_code` directory.
  ```
    bash bash_code/preProcess.pipeline \
      variant_calling \
      genome/acerv_genome.fasta \
      PE \
      DNA \
      140 \
      variant_calling/raw_reads/*fastq.gz
  ```

## 2. Call SNPs
SNP calling requirements:
  - Minimum Depth = 0.33 * Number of Individuals
  - Maximum Depth = Mean Depth + 4 * Standard Deviation
  - Minimum Number of Individuals = 0.5 * Number of Individuals
  - Minimum Base Quality = 30
  - Minimum Minor Allele Frequency = 0.05
  - Minimum Mapping Quality = 30
  - SNP p-value = 1e-6
  ```
  bash bash_code/genotype.pipeline \
    variant_calling/genotyping \
    genome/acerv_genome.fasta \
    0.33 \
    4 \
    0.5 \
    30 \
    0.05 \
    30 \
    variant_calling/repair/*fastq.gz
  ```

## 3. Merge metadata together
Combine all metadata for samples into useable form using: `r_code/0 - assemble_metadata.R`

## 4. Estimate Disease Resistance
Model disease trials and estimate genotype level disease resistance using: `1 - calculate_disease_resistance.R`

## 5. Preprocess Called SNPs
Filter poorly represented SNPs and SNPs with low call confidence using: `2 - preprocess_genetics.R`

## 6. Identify Genetic Population Structure
Follow code in: `3 - genetic_structure-preLFMM.R` which include code to be used interactively on an HPC to us ANGSD with the post-filtering set of individuals and loci to identify population structure

## 7. Perform Latent-Factor Mixed Model analysis for GWAS
Run code: `4 - runLFMM.R`
