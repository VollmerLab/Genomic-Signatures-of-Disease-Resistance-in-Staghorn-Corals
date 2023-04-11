#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=fastQC
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

# Script that takes in a list of fastq files and Runs FASTQC & MultiQC

outdir=${1};shift 1 #Either DNA or RNA (changes how mapping is done due to splice sites)
fastq_files="$@"

mkdir -p ${outdir}
cd ${outdir} #move to outdir for if temp files are made

#Run FastQC & MultiQC on fastq files
module load anaconda3; source activate fastqc
mkdir -p ${outdir}/FastQC_fastq

fastqc \
  -t ${SLURM_CPUS_PER_TASK} \
  -o ${outdir}/FastQC_fastq \
  ${fastq_files}

conda deactivate; module load anaconda3; source activate multiqc
multiqc ${outdir}/FastQC_fastq -n ${outdir}/fastq_report
conda deactivate
