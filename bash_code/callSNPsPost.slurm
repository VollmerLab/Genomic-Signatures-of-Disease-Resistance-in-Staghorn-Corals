#!/bin/bash

#SBATCH --job-name=mergeAngsd
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

outdir=${1}
scriptDir=/work/vollmer/software/jds_scripts

#Get all types of output from angsd
output_types=$(ls ${outdir}/contig_snps/*gz | sed -e 's/.*\(\..*\.gz\)/\1/g' | sort -u)
out_name=$(ls ${outdir}/contig_snps/*gz | sed -e 's/.*contig_.*_//g' | sed -e 's/\..*//g' | sort -u)

for output in ${output_types}; do
  echo ${output}
  cat ${outdir}/contig_snps/*${output} > ${outdir}/contig_snps/${out_name}_TMP${output}

  #remove rows that match first row excluding first row
  header_line=$(zcat ${outdir}/contig_snps/${out_name}_TMP${output} | head -n 1)
  gzip -c -d ${outdir}/contig_snps/${out_name}_TMP${output} | sed "s/${header_line}//g" | sed '/^$/d' | sed "1 i ${header_line}" | \
    gzip -c > ${outdir}/${out_name}${output}
  rm ${outdir}/contig_snps/${out_name}_TMP${output}
done

rm -rf ${outdir}/contig_snps

module load anaconda3; source activate pcangsd
cd ${outdir}
pcangsd \
  -b ${outdir}/genolike.beagle.gz \
  -t ${SLURM_CPUS_PER_TASK} \
  -o output.pcangsd \
  --tree \
  --maf_save \
  --pi_save \
  --admix \
  --admix_auto 10
conda deactivate

module load R/4.2.0
Rscript ${scriptDir}/r_utils/post_genotyping.R ${outdir}
