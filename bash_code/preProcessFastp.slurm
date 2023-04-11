#!/bin/bash

#SBATCH --job-name=fp
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load anaconda3; source activate fastp

outdir=${1}
se_or_pe=${2} #Either SE or PE
minLength=${3};shift 3
fastq_files="$@"


mkdir -p ${outdir}/failed
cd ${outdir} #move to outdir if any temp file writing is needed

#Checking Conditions
if [ ${se_or_pe} != "SE" ] && [ ${se_or_pe} != "PE" ]
then
  echo "se_or_pe (\$2) must be one of SE or PE. Instead we got ${se_or_pe}"
  exit 0
else
  echo "Running with ${se_or_pe} data"
fi

#
# Sort out post-array & choose sample to run in array
tmp=( ${fastq_files} )
indir=${tmp%/*}
sample_names=( $(ls ${fastq_files} | sed -e 's!.*/!!' | sed -e 's/[\._][rR][12]\.*\..*//' | sort -u) )
extensions=( $(ls ${fastq_files} | sed -e 's!.*/!!' | sed -e 's/^.*\([\._][rR][12]\)/\1/g' | sort -u) )
number_of_individuals=${#sample_names[@]}
number_of_extensions=${#extensions[@]}

sample_name=${sample_names[${SLURM_ARRAY_TASK_ID}]}

echo "Running FastP on ${sample_name}"

#Sort out single vs paired end
if [ ${se_or_pe} == "SE" ]
then
  in_file="--in1 ${indir}/${sample_name}${extensions[0]}"
  out_file="--out1 ${outdir}/${sample_name}-fp${extensions[0]}"

  if [ -f "${outdir}/${sample_name}-fp${extensions[0]}" ]; then
    echo "Output files for ${sample_name} already created"
    exit 0
  fi
else
  in_file="--in1 ${indir}/${sample_name}${extensions[0]} --in2 ${indir}/${sample_name}${extensions[1]}"
  out_file="--out1 ${outdir}/${sample_name}-fp${extensions[0]} --out2 ${outdir}/${sample_name}-fp${extensions[1]}"

  if [ -f "${outdir}/${sample_name}-fp${extensions[0]}" ] && [ -f "${outdir}/${sample_name}-fp${extensions[1]}" ]; then
    echo "Output files for ${sample_name} already created"
    exit 0
  fi
fi

#Run FastP
fastp \
  ${in_file} \
  ${out_file} \
  --failed_out "${outdir}/failed/${sample_name}-fail.fq.gz" \
  -w ${SLURM_CPUS_PER_TASK} \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 40 \
  --length_required ${minLength} \
  --low_complexity_filter \
  --complexity_threshold 30 \
  --detect_adapter_for_pe \
  --cut_tail \
  --cut_tail_window_size 1 \
  --cut_tail_mean_quality 20 \
  --cut_front \
  --cut_front_window_size 1 \
  --cut_front_mean_quality 20 \
  --cut_right \
  --cut_right_window_size 10 \
  --cut_right_mean_quality 20 \
  --trim_poly_g \
  --poly_g_min_len 10 \
  --trim_poly_x \
  --correction \
  --unqualified_percent_limit 40 \
  -h ${outdir}/${sample_name}.fastp.html \
  -j ${outdir}/${sample_name}.fastp.json \
  --report_title "fastp"
