#!/bin/bash

#SBATCH --job-name=LD
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load singularity
CONTAINER=/work/vollmer/software/ngsld.sif

#outdir=/scratch/j.selwyn/linkage
#ARRAY_ID=0

outdir=${1}
ARRAY_ID=${SLURM_ARRAY_TASK_ID}

cd ${outdir}

#get all position & genotype files
all_position_files=( $(ls ${outdir}/*_posFile.pos) )
all_genotype_files=( $(ls ${outdir}/*_genotypes.tsv.gz) )

#Both are sorted alphabetically so pick pair and run
posFile=${all_position_files[${ARRAY_ID}]}
genotypes=${all_genotype_files[${ARRAY_ID}]}

#Get chromosome
chromosome=$(echo ${genotypes##*/} | sed "s/_genotypes.tsv.gz//g")
echo ${chromosome}

#Get necessary info for NGSLD
nInd=$(zcat ${genotypes} | head -n 1 | tr $'\t' '\n' | wc -l)
nLoc=$(wc -l < ${posFile})

singularity exec --bind /work,/scratch,/tmp ${CONTAINER} \
  ngsLD \
    --n_threads ${SLURM_CPUS_PER_TASK} \
    --geno ${genotypes} \
    --pos ${posFile} \
    --n_ind ${nInd} \
    --n_sites ${nLoc} \
    --max_kb_dist 0 \
    --max_snp_dist 0 \
    --extend_out \
    --ignore_miss_data \
    --outH ${outdir}/${chromosome}.ld

#gzip ${outdir}/${chromosome}.ld
