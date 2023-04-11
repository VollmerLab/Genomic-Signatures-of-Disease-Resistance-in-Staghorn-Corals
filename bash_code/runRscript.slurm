#!/bin/bash

#SBATCH --job-name=Rscript
#SBATCH --partition=short
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

to_run=$1
echo "${to_run}"

#module load R/4.2.0

#Rscript "${to_run}" "${@:2}"

module load singularity/3.5.3

#RSTUDIO_IMAGE="/shared/container_repository/rstudio/rocker-tidyverse-4.2.1.sif"
RSTUDIO_IMAGE="/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif"
singularity run -B "/home:/home,/scratch:/scratch,/work:/work" ${RSTUDIO_IMAGE} Rscript "${to_run}" "${@:2}"

#singularity shell --bind /work,/scratch,/tmp ${RSTUDIO_IMAGE}
