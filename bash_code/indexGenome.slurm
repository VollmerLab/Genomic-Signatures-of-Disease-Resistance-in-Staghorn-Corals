#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=indexGenome
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

# Script that takes in a genome and indexes it as either DNA (with bwa index) or RNA (with gmap_build)

outdir=${1}
ref_genome=${2}
nucleotide=${3}


if [ ${nucleotide} != "DNA" ] && [ ${nucleotide} != "RNA" ]
then
  echo "Nucleotide (\$3) must be one of DNA or RNA"
  exit 0
else
  echo "Running with Nucleotide ${nucleotide}"
fi


mkdir -p ${outdir}
cd ${outdir} #move to outdir for if temp files are made
ln -s ${ref_genome} ${outdir}/reference.fa

module load bwa; module load samtools/1.10

if [ ${nucleotide} == "DNA" ]
then
  bwa index ${outdir}/reference.fa

else
  #Build Genome database
  module load singularity
  CONTAINER=/work/vollmer/software/gmap.sif

  singularity exec --bind /work,/scratch,/tmp ${CONTAINER} \
    gmap_build \
      --dir=${outdir}/genomeDB \
      -d=genomeDB \
      -t ${SLURM_CPUS_PER_TASK} \
      ${outdir}/reference.fa

fi

grep -o '^>.*' ${outdir}/reference.fa | sed -e 's/>//g' > ${outdir}/reference.contigs
samtools faidx ${outdir}/reference.fa
