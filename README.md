# Genomic-Signatures-of-Disease-Resistance-in-Endangered-Staghorn-Corals-
Code associated with manuscript analyzing the genomic basis of disease resistance in staghorn corals

Data needed to analyze this data from scratch:
  - the assembled genome (NCBI BioProj: PRJNA948411)
  - the *A. cervicornis* WGS reads to be used (NCBI BioProj: PRJNA950067)
  - Supplemental WGS for *A. palmata* and *A. prolifera* (NCBI BioProj: PRJNA473816)

Can skip SNP preprocessing and calling by starting with the already identified SNPs `Data/genotypes.csv.gz` & `Data/genotype_probabilities.csv.gz`. If so skip preprocessing and SNP identifications steps.

## 1. Preprocess Sequences
  - all (bash scripts)[bash_code] are designed to run on a SLURM based HPC
  - scripts should be modified such that the `${scriptDir}` variable is assigned at the beginning of scripts to the `bash_code` directory.
  ```
    bash bash_pipeline/preProcess.pipeline \
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
bash /work/vollmer/software/jds_scripts/genotype.pipeline \
  /scratch/j.selwyn/variant_calling/k2_18-October-2022/genotyping \
  /scratch/j.selwyn/genome_assemblies/k2_18-October-2022/k2HybridAssembly.fasta \
  0.33 \
  4 \
  0.5 \
  30 \
  0.05 \
  30 \
  /scratch/j.selwyn/variant_calling/k2_18-October-2022/repair/*fastq.gz
```
