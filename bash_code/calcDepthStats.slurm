#!/bin/bash

#SBATCH --job-name=depthStats
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G

outdir=${1}
bam_list=${2}
MINQ=${3}
MINMAPQ=${4}
speed=${5:-fast} #accurate or fast

scriptDir=/work/vollmer/software/jds_scripts

#Scale maximum depth by mean/sd of read depth
mkdir -p ${outdir}
cd ${outdir}

module load samtools
module load R

#Checking Conditions
if [ ${speed} != "fast" ] && [ ${speed} != "accurate" ]
then
  echo "speed (\$5) must be one of fast or accurate. Instead we got ${speed}"
  exit 0
else
  echo "Running with ${speed} speed setting"
fi

if [ ${speed} == "accurate" ]
then
  #Currently this takes practically infinite time and doesn't calculate what you want it to. It calculates the mean and sd of the mean depths per individual you probably want the mean and sd of the total depths of each locus
  samtools depth -f ${bam_list} -d 0 -q ${MINQ} -Q ${MINMAPQ} | cut -f2- | cut -f2- > ${outdir}/bam.depth

  mean_depth=$(cat ${outdir}/bam.depth | Rscript -e 'cat(mean(rowMeans(read.table("stdin", sep = "\t"), na.rm = TRUE), na.rm = TRUE))') &
  waitID=$!
  sd_depth=$(cat ${outdir}/bam.depth | Rscript -e 'cat(sd(rowMeans(read.table("stdin", sep = "\t"), na.rm = TRUE), na.rm = TRUE))')
  wait ${waitID}

  echo ${mean_depth} | cat > ${outdir}/depth.stats
  echo ${sd_depth} | cat >> ${outdir}/depth.stats
else
  #https://bioinformatics.stackexchange.com/questions/13722/calculating-average-coverage-for-bam-files-sequence-data
  readarray -t all_bam < ${bam_list}

  nInd=${#all_bam[@]}

  for ((i=0; i<nInd; i++))
  do
      echo ${all_bam[${i}]}

      this_bam=${all_bam[${i}]}
      depth_name=$(echo ${this_bam##*/} | sed -e 's/_sorted\.bam/\.depths/g')

      samtools idxstats ${all_bam[$i]} > ${outdir}/${depth_name}
  done

  Rscript ${scriptDir}/r_utils/bam_mean_sd_depth_total.R \
    ${outdir} \
    ${outdir}/depth.stats \
    140
  #rm ${outdir}/*depths
fi
