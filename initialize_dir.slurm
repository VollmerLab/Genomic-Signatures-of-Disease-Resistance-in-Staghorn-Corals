#!/bin/bash

#SBATCH --job-name=initDIR
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

## Download Genome
module load anaconda3; source activate ncbi_datasets
mkdir genome; cd genome
datasets download genome accession GCA_032359415.1 \
  --include genome,gff3,rna,cds,protein,seq-report \
  && unzip ncbi_dataset.zip \
  && rm -rf ncbi_dataset.zip \
  && mv ncbi_dataset/data/GCA_032359415.1/GCA_032359415.1_NEU_Acer_K2_genomic.fna ./acerv_genome.fasta \
  && mv ncbi_dataset/data/GCA_032359415.1/* ./ \
  && rm -rf ncbi_dataset \
  && cd ..

#Process annotations into useful formats for later
bash bash_code/runRscript.slurm bash_code/r_utils/processGFF.R genome/genomic.gff genome/acerv_genome.fasta ${SLURM_CPUS_PER_TASK}

## Get SRR Numbers and Names
module load anaconda3; source activate sratoolkit
esearch -db sra -query PRJNA950067 |\
  efetch -format runinfo |\
  cut -d "," -f 1,12 >\
  SRR.numbers

esearch -db sra -query PRJNA473816 |\
  efetch -format runinfo |\
  grep "2018-10-22\|2021-07-31" |\
  grep -v "SRR7235995" |\
  grep -v "SRR7235983" |\
  cut -d "," -f 1 >\
  SRR.numbers_baum

vollmer_reads=( $(cat SRR.numbers | cut -d "," -f 1 | tail -n +2) )
baum_reads=( $(cat SRR.numbers_baum) )
echo "${#vollmer_reads[@]} Vollmer lab Samples"
echo "${#baum_reads[@]} Baum lab Samples"

SRR_reads=( "${vollmer_reads[@]}" "${baum_reads[@]}" )
echo "${#SRR_reads[@]} Total Samples"

#Download Reads
mkdir variant_calling/raw_reads
module load parallel
parallel --jobs ${SLURM_CPUS_PER_TASK} \
  --tmpdir raw_reads \
  "fastq-dump -O variant_calling/raw_reads --split-files --origfmt --gzip {}" ::: ${SRR_reads[@]}
conda deactivate

#Rename from NCBI name to names matching metadata
bash bash_code/runRscript.slurm bash_code/r_utils/srr_to_sequenceID.R
#Extend to include the renaming portion

#Rename forward/reverse for preprocessing scripts
rename _2 .r2 raw_reads/*fastq.gz
rename _1 .r1 raw_reads/*fastq.gz
