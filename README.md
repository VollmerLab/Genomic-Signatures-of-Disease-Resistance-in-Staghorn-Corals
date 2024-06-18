# Genomic-Signatures-of-Disease-Resistance-in-Endangered-Staghorn-Corals-
Code associated with manuscript analyzing the genomic basis of disease resistance in staghorn corals

Data needed to analyze this data from scratch:
  - the assembled genome & annotations (NCBI BioProj: PRJNA948411)
  - the *A. cervicornis* WGS reads to be used (NCBI BioProj: PRJNA950067)
  - Supplemental WGS for *A. palmata* and *A. prolifera* (NCBI BioProj: PRJNA473816)

There are two approaches to reproducing this analysis. The first involves preprocessing Illumina reads and identifying SNPs in the genome. The second option can skip ahead to already identified SNPs `Data/genotypes.csv.gz` & `Data/genotype_probabilities.csv.gz`. To start with already identified SNPs skip the initialization, preprocessing, and SNP calling steps. 

All code from step 3 onwards is designed to be run interactively inside RStudio (or your prefered IDE) except where noted within the code.

## Software Used
Here is the list of software installed and used in this complete analysis. Because of varying computer architectures and methods of installation you will need to look at each script and change the functions used to target where your specific installation of each of these software packages resides. These software are generally loaded as either a singularity container (e.g. `singularity run -B "/home:/home,/scratch:/scratch,/work:/work" ${RSTUDIO_IMAGE}` in script [`bash_code\runRscript.slurm`](bash_code\runRscript.slurm)), a conda environment (e.g. `source activate sratoolkit` in [initialize_dir.sh](initialize_dir.sh)), or a module (e.g. `module load parallel` in [initialize_dir.sh](initialize_dir.sh)). 

- [`initialize_dir.sh`](initialize_dir.sh)
  - [sratoolkit](https://github.com/ncbi/sra-tools)
  - [entrez-direct](https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/versions/21.6.20240308/README)
  - [gnu parallel](https://www.gnu.org/software/parallel/)
  - [R](https://cran.r-project.org/)
- [`bash_code/preProcess.pipeline`](bash_code/preProcess.pipeline)
- [`bash_code\assessSequencing.pipeline`](bash_code\assessSequencing.pipeline)
- [`bash_code\assessFastq.slurm`](bash_code\assessFastq.slurm)
  - [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  - [multiqc](https://multiqc.info/)
- [`bash_code\indexGenome.slurm`](bash_code\indexGenome.slurm)
  - [bwa](https://github.com/lh3/bwa)
  - [samtools](http://www.htslib.org/)
  - [gmap](http://research-pub.gene.com/gmap/)
- [`bash_code\assessMapping.slurm`](bash_code\assessMapping.slurm)
  - [bwa](https://github.com/lh3/bwa)
  - [samtools](http://www.htslib.org/)
  - [gmap](http://research-pub.gene.com/gmap/)
  - [R](https://cran.r-project.org/)
- [`bash_code\preProcessFastp.slurm`](bash_code\preProcessFastp.slurm)
  - [fastp](https://github.com/OpenGene/fastp)
- [`bash_code\preProcessFQScreen.slurm`](bash_code\preProcessFQScreen.slurm)
  - Database to filter against. Database construction described [here](https://stevenwingett.github.io/FastQ-Screen/). Customized to include *Symbiodinium* genomes (see [Table S3](Manuscript/Vollmer et al 2023 Supplement.pdf) for accession numbers)
    - Change `CONFFILE` variable in the [`bash_code\preProcessFQScreen.slurm`](bash_code\preProcessFQScreen.slurm) script to target your configuration file for fastq screen
  - [gnu parallel](https://www.gnu.org/software/parallel/)
  - [bowtie](https://bowtie-bio.sourceforge.net/index.shtml)
  - [fastq_screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
- [`bash_code\preProcessRepair.slurm`](bash_code\preProcessRepair.slurm)
  - [bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
  - [pigz](https://zlib.net/pigz/)
- [`bash_code\genotype.pipeline`](bash_code\genotype.pipeline)
- [`bash_code\calcDepthStats.slurm`](bash_code\calcDepthStats.slurm)
  - [samtools](http://www.htslib.org/)
  - [R](https://cran.r-project.org/)
- [`bash_code\callSNPs.slurm`](bash_code\callSNPs.slurm)
  - [angsd](https://www.popgen.dk/angsd/index.php/ANGSD)
- [`bash_code\callSNPsPost.slurm`](bash_code\callSNPsPost.slurm)
  - [pcangsd](http://www.popgen.dk/software/index.php/PCAngsd)
  - [R](https://cran.r-project.org/)
- [`bash_code\runRscript.slurm`](bash_code\runRscript.slurm)
  - [R](https://cran.r-project.org/)
- [`bash_code\runNGSLD_array.slurm`](bash_code\runNGSLD_array.slurm)
  - [ngsLD](https://github.com/fgvieira/ngsLD)

## R packages Used
- [adegenet](https://cran.r-project.org/web/packages/adegenet/index.html)
- [bioseq](https://cran.r-project.org/web/packages/bioseq/index.html)
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
- [broom](https://cran.r-project.org/web/packages/broom/index.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
- [dartR](https://cran.r-project.org/web/packages/dartR/index.html)
- [forcats](https://cran.r-project.org/web/packages/forcats/index.html)
- [ggdist](https://cran.r-project.org/web/packages/ggdist/index.html)
- [gghalves](https://cran.r-project.org/web/packages/gghalves/index.html)
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html)
- [gratia](https://cran.r-project.org/web/packages/gratia/index.html)
- [hierfstat](https://cran.r-project.org/web/packages/hierfstat/index.html)
- [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html)
- [janitor](https://cran.r-project.org/web/packages/janitor/index.html)
- [LEA](https://bioconductor.org/packages/release/bioc/html/LEA.html)
- [lme4](https://cran.r-project.org/web/packages/lme4/index.html)
- [lmerTest](https://cran.r-project.org/web/packages/lmerTest/index.html)
- [magrittr](https://cran.r-project.org/web/packages/magrittr/index.html)
- [mgcv](https://cran.r-project.org/web/packages/mgcv/index.html)
- [multidplyr](https://cran.r-project.org/web/packages/multidplyr/index.html)
- [patchwork](https://cran.r-project.org/web/packages/patchwork/index.html)
- [poppr](https://cran.r-project.org/web/packages/poppr/index.html)
- [progressr](https://cran.r-project.org/web/packages/progressr/index.html)
- [readr](https://cran.r-project.org/web/packages/readr/index.html)
- [reticulate](https://cran.r-project.org/web/packages/reticulate/index.html)
- [SNPRelate](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)
- [stringr](https://cran.r-project.org/web/packages/stringr/index.html)
- [tibble](https://cran.r-project.org/web/packages/tibble/index.html)
- [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)
- [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
- [xml2](https://cran.r-project.org/web/packages/xml2/index.html)

## 0. Initialize Directory
- all [bash scripts](bash_code) are designed to run on a SLURM based HPC
- scripts should be modified such that the `${scriptDir}` variable is assigned at the beginning of scripts to the `bash_code` directory.
- Each bash script will need to be modified depending on where/how the software is installed. 
- Clone the github repository and then run `sbatch initialize_dir.slurm` from within that directory to initialize it with all the NCBI sequences & genome. 
  - This will download ~1 TB of data from NCBI including:
    - *Acropora cervicornis* genome 
    - *Acropora cervicornis* genome annotations 
    - *Acropora cervicornis* short-read data

## 1. Preprocess Sequences
  ```
    bash bash_code/preProcess.pipeline \
      $(pwd)/variant_calling \
      $(pwd)/genome/acerv_genome.fasta \
      PE \
      DNA \
      140 \
      $(pwd)/variant_calling/raw_reads/*fastq.gz
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
    $(pwd)/variant_calling/genotyping \
    $(pwd)/genome/acerv_genome.fasta \
    0.33 \
    4 \
    0.5 \
    30 \
    0.05 \
    30 \
    $(pwd)/variant_calling/repair/*fastq.gz
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
