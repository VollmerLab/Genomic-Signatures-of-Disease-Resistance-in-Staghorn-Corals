#!/bin/bash

#SBATCH --job-name=fqscrn
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=250G

module load parallel
module load bowtie
module load anaconda3; source activate fastq_screen

outdir=${1}; shift 1
fastq_files="$@"

TMPDIR=${outdir}/fastqscrn_${SLURM_JOB_USER}_${SLURM_JOBID}_${SLURM_ARRAY_TASK_ID}
CONFFILE=/work/vollmer/Databases/fastq_screenDatabases/fastq_screen_2022-August-24.conf
SUBSET=0
FILTER=00000000000000000000000000
ALIGNER=bowtie2

# Sort out post-array & choose sample to run in array
tmp=( ${fastq_files} )
indir=${tmp%/*}
full_file=( $(ls ${fastq_files} | sed -e 's!.*/!!') )

in_file=${full_file[${SLURM_ARRAY_TASK_ID}]}

#WHat should output look like
#out_file=$(echo ${in_file} | sed -e 's/\(\.f.*\.gz\)/.tagged_filter\1/g')
#rename_out=$(echo ${out_file} | sed -e 's/\.tagged_filter//' | sed -e 's/\(\.[rR][12]\.f.*\.gz\)/-fqscrn\1/g')

out_file=$(echo ${in_file} | sed -e 's/\(\.f.*$\)/.tagged_filter\1/g')
rename_out=$(echo ${out_file} | sed -e 's/\.tagged_filter//' | sed -e 's/\(\.[rR][12]\.f.*$\)/-fqscrn\1/g')

echo "Running FQScreen on ${in_file}"

if [ -f "${outdir}/${rename_out}" ]; then
	echo "Output files for ${in_file} already created"
	exit 0
fi

#Run FQScreen
#Added while loop because occassionally file wouldn't exist after run.
#This is possibly a result of changing the script file mid-run but its unclear. So just loop until output file exists
while [ ! -f ${outdir}/${rename_out} ]
do
	mkdir -p ${TMPDIR}
	fastq_screen \
		--aligner ${ALIGNER} \
		--conf ${CONFFILE} \
		--threads ${SLURM_CPUS_PER_TASK} \
		--tag \
		--force \
		--filter ${FILTER} \
		--subset ${SUBSET} \
		--outdir ${TMPDIR} \
		${indir}/${in_file}

	mv ${TMPDIR}/* ${outdir}/
	rm -rf ${TMPDIR}
	mv ${outdir}/${out_file} ${outdir}/${rename_out}
done
