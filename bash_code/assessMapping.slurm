#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=mapGenome
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

# Script that maps a list of FASTA files to a genome as either DNA (with bwa mem) or RNA (with gsnap)

outdir=${1}
se_or_pe=${2} #Either SE or PE
nucleotide=${3};shift 3 #Either DNA or RNA (changes how mapping is done due to splice sites)
fastq_files="$@"

cd ${outdir} #move to outdir if any temp file writing is needed

#Checking Conditions
if [ ${se_or_pe} != "SE" ] && [ ${se_or_pe} != "PE" ]
then
  echo "se_or_pe (\$3) must be one of SE or PE"
  exit 0
else
  echo "Running with ${se_or_pe} data"
fi

if [ ${nucleotide} != "DNA" ] && [ ${nucleotide} != "RNA" ]
then
  echo "Nucleotide (\$4) must be one of DNA or RNA"
  exit 0
else
  echo "Running with Nucleotide ${nucleotide}"
fi

# Sort out post-array & choose sample to run in array
tmp=( ${fastq_files} )
indir=${tmp%/*}
sample_names=( $(ls ${fastq_files} | sed -e 's!.*/!!' | sed -e 's/[\._][rR][12]\.*\..*//' | sort -u) )
extensions=( $(ls ${fastq_files} | sed -e 's!.*/!!' | sed -e 's/^.*\([\._][rR][12]\)/\1/g' | sort -u) )
number_of_individuals=${#sample_names[@]}
number_of_extensions=${#extensions[@]}

sample_name=${sample_names[${SLURM_ARRAY_TASK_ID}]}

echo "Mapping ${sample_name} to Genome"

#Map Reads to Reference Genome
module load bwa; module load samtools/1.10

if [ ${se_or_pe} == "SE" ]
then
  in_file=${indir}/${sample_name}${extensions[0]}
else
  in_file="${indir}/${sample_name}${extensions[0]} ${indir}/${sample_name}${extensions[1]}"
fi

#Used for testing. If the indexed bam file already esists then don't remake it
if [ ! -f ${outdir}/${sample_name}_sorted.bam.bai ]
  then
    if [ ${nucleotide} == "DNA" ]
    then

      bwa mem \
        -t ${SLURM_CPUS_PER_TASK} \
        ${outdir}/reference.fa \
        ${in_file} > ${outdir}/${sample_name}.sam
    else
      module load singularity
      CONTAINER=/work/vollmer/software/gmap.sif

      singularity exec --bind /work,/scratch,/tmp ${CONTAINER} \
        gsnap \
          --gunzip \
          -D ${outdir}/genomeDB  \
          -d genomeDB \
          -t ${SLURM_CPUS_PER_TASK} \
          ${forceSE} \
          --format=sam \
          -B 4 \
          -N 1 \
          -M 2 \
          -n 10 \
          --output-file=${outdir}/${sample_name}.sam \
          ${in_file}

    fi

    #Convert to BAM and Index
    MEMORY_PER_TASK=$(echo "scale=2 ; 0.8 * ${SLURM_MEM_PER_NODE} / ${SLURM_CPUS_PER_TASK}" | bc)
    MEMORY_PER_TASK=${MEMORY_PER_TASK%.*}

    samtools view -Sb -@ ${SLURM_CPUS_PER_TASK} ${outdir}/${sample_name}.sam > ${outdir}/${sample_name}.bam
    echo "Finished SAM -> BAM!"

    samtools sort -@ ${SLURM_CPUS_PER_TASK} -m ${MEMORY_PER_TASK}M -O BAM ${outdir}/${sample_name}.bam > ${outdir}/${sample_name}_sorted.bam
    echo "Finished BAM sorting!"

    samtools index -b -@ ${SLURM_CPUS_PER_TASK} ${outdir}/${sample_name}_sorted.bam
    echo "Finished BAM indexing!"
fi

#Post-process mapping and clean-up files
module load R
samtools flagstat -@ ${SLURM_CPUS_PER_TASK} -O tsv ${outdir}/${sample_name}_sorted.bam | \
  Rscript -e 'write.table(t(rev(read.table("stdin", sep = "\t")))[c(1,3), ], sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)' > ${outdir}/${sample_name}_mapstats2.tmp1

printf %s\\n sample ${sample_name} > ${outdir}/${sample_name}.tmp2
paste ${outdir}/${sample_name}.tmp2 ${outdir}/${sample_name}_mapstats2.tmp1 > ${outdir}/${sample_name}_mapstats.tsv
rm ${outdir}/${sample_name}*tmp*
rm ${outdir}/${sample_name}.sam
rm ${outdir}/${sample_name}.bam
echo "Finished Everything!"
