#!/bin/bash

#SBATCH --job-name=angsd
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

outdir=${1}
bam_list=${2}
contig_list=${3}
array_size=${4}

MINDP=${5} # Minimum depth filter - e.g. 0.33 x number of individuals - code as result of  0.33 * n_individuals (as integer)
MAXDP=${6} # Maximum depth filter - e.g = mean depth + 4 s.d. - code as 4 (e.g. the number of standard deviations to multiply by)
MININD=${7} # Minimum individual filter - e.g. 50% of individuals - code as result of  0.5 * n_individuals (as integer)
MINQ=${8} # Minimum quality filter - e.g. 20
MINMAF=${9} # Minimum minor allele frequency filter - e.g. 0.05
MINMAPQ=${10:-20} # Minimum mapping quality (alignment score) filter, default value is 20

n_ind=$(sed -n '$=' ${bam_list})

mkdir -p ${outdir}
cd ${outdir}

module load anaconda3; source activate angsd
#get reference genome
ref_genome=${outdir}/../reference.fa

#Scale maximum depth by mean/sd of read depth
readarray -t depth_stats < ${outdir}/../depth.stats
MAXDP=$(echo "scale=0 ; (${depth_stats[0]} + ${MAXDP} * ${depth_stats[1]})/1" | bc)

#Announce SNP filter settings
echo "Minimum total sequencing depth for a locus to be included: ${MINDP}"
echo "Maximum total sequencing depth for a locus to be included: ${MAXDP}"
echo "Minimum number of individuals locus must be found in to be included: ${MININD}"
echo "Minimum number base quality to include locus: ${MINQ}"
echo "Minimum minor allele frequency to include locus: ${MINMAF}"
echo "Minimum alignment score to include locus: ${MINMAPQ}"

#Find which contigs to use
readarray -t all_contigs < ${contig_list}

split_jobs () {
  #https://www.geeksforgeeks.org/split-the-number-into-n-parts-such-that-difference-between-the-smallest-and-the-largest-part-is-minimum/
  X=$1
  N=$2

  out=""

  if [ $((X % N)) -eq 0 ]
    then
      for ((i=1; i<=N; i++))
      do
        out="${out} $((X / N))"
      done

    else
      zp=$((N - (X % N)))
      pp=$((X / N))

      for ((i=1; i<=N; i++))
      do
        if [ ${i} -gt ${zp} ]
        then
          out="${out} $((pp + 1))"
        else
          out="${out} ${pp}"
        fi
      done
    fi
  echo "${out}"
}

contig_bins=( $(split_jobs ${#all_contigs[@]} ${array_size}) )

ARRAY_ID=${SLURM_ARRAY_TASK_ID}
starting_index=$(echo "0 ${contig_bins[@]:0:$((ARRAY_ID))} 0" | sed -e 's/  */+/g' | bc)
contig_names=${all_contigs[@]:${starting_index}:${contig_bins[${ARRAY_ID}]}}

echo ${contig_names}

#Run ANGSD
for contig_use in ${contig_names}; do
  echo "Starting ${contig_use}"
  angsd \
    -r ${contig_use}: \
    -nThreads ${SLURM_CPUS_PER_TASK} \
    -out ${outdir}/${contig_use}_genolike \
    -bam ${bam_list} \
    -ref ${ref_genome} \
    -remove_bads 1 \
    -uniqueOnly 1 \
    -only_proper_pairs 1 \
    -C 50 \
    -minMapQ ${MINMAPQ} \
    -minQ ${MINQ} \
    -doCounts 1 \
    -setMinDepth ${MINDP} \
    -setMaxDepth ${MAXDP} \
    -setMinDepthInd 5 \
    -minInd ${MININD} \
    -GL 2 \
    -doGlf 2 \
    -doMajorMinor 1 \
    -skipTriallelic 1 \
    -doMaf 3 \
    -minMaf ${MINMAF} \
    -SNP_pval 1e-6 \
    -doPost 1 \
    -doGeno 19 \
    -geno_minDepth 5 \
    -geno_minMM 0.1
done
