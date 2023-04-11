#!/bin/bash

#SBATCH --job-name=repair
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

outdir=${1}
se_or_pe=${2}; shift 2 #Either SE or PE
fastq_files="$@"

mkdir -p ${outdir}/orphans

#Checking Conditions
if [ ${se_or_pe} != "SE" ] && [ ${se_or_pe} != "PE" ]
then
  echo "se_or_pe (\$3) must be one of SE or PE Instead we got ${se_or_pe}"
  exit 0
else
  echo "Running with ${se_or_pe} data"
fi

# Sort out post-array & choose sample to run in array
tmp=( ${fastq_files} )
indir=${tmp%/*}
sample_names=( $(ls ${fastq_files} | sed -e 's!.*/!!' | sed -e 's/[\._][rR][12]\.*\..*//' | sort -u) )
extensions=( $(ls ${fastq_files} | sed -e 's!.*/!!' | sed -e 's/^.*\([\._][rR][12]\)/\1/g' | sort -u) )
number_of_individuals=${#sample_names[@]}
number_of_extensions=${#extensions[@]}

sample_name=${sample_names[${SLURM_ARRAY_TASK_ID}]}

echo "Running Repairs on ${sample_name}"


if [ ${se_or_pe} == "PE" ]
then
  if [ -f "${outdir}/${sample_name}-repr${extensions[0]}" ] && [ -f "${outdir}/${sample_name}-repr${extensions[1]}" ]; then
    echo "Output files for ${sample_name} already created"
    exit 0
  fi

  module load anaconda3; source activate bbmap
	# repair the paired end files
  repair.sh \
    in1=${indir}/${sample_name}${extensions[0]} \
    in2=${indir}/${sample_name}${extensions[1]} \
    out1="${outdir}/${sample_name}-reprTMP${extensions[0]}" \
    out2="${outdir}/${sample_name}-reprTMP${extensions[1]}" \
    outs="${outdir}/orphans/${sample_name}-repr-orph${extensions[0]}" \
    repair

    #if zipped then rezip
    if [ $(echo ${extensions[0]} | grep -c '\.gz$') == 1 ]
    then
      conda deactivate; module load anaconda3; source activate pigz
      pigz -p ${SLURM_CPUS_PER_TASK} -c -d "${outdir}/${sample_name}-reprTMP${extensions[0]}" | \
      sed -e '/^@.*\/[12]#/{s!#! #!g}' | \
      pigz -p ${SLURM_CPUS_PER_TASK} -c > "${outdir}/${sample_name}-repr${extensions[0]}" &

      waitID=$!

      pigz -p ${SLURM_CPUS_PER_TASK} -c -d "${outdir}/${sample_name}-reprTMP${extensions[1]}" | \
      sed -e '/^@.*\/[12]#/{s!#! #!g}' | \
      pigz -p ${SLURM_CPUS_PER_TASK} -c > "${outdir}/${sample_name}-repr${extensions[1]}"
      wait ${waitID}
    else
      cat "${outdir}/${sample_name}-reprTMP${extensions[0]}" | \
      sed -e '/^@.*\/[12]#/{s!#! #!g}' >"${outdir}/${sample_name}-repr${extensions[0]}" &

      waitID=$!

      cat "${outdir}/${sample_name}-reprTMP${extensions[1]}" | \
      sed -e '/^@.*\/[12]#/{s!#! #!g}' >"${outdir}/${sample_name}-repr${extensions[1]}"
      wait ${waitID}
    fi
    
else

  if [ -f "${outdir}/${sample_name}-repr${extensions[0]}" ]; then
    echo "Output files for ${sample_name} already created"
    exit 0
  fi

  module load anaconda3; source activate pigz
  #Fix issue if there is no space between "/[12]" and "#" seperating which read it is and fltrscrn

  pigz -p ${SLURM_CPUS_PER_TASK} -c -d ${indir}/${sample_name}${extensions[0]} | \
  sed -e '/^@.*\/[12]#/{s!#! #!g}' | \
  pigz > "${outdir}/${sample_name}-repr${extensions[0]}"
fi

rm ${outdir}/${sample_name}*TMP*
