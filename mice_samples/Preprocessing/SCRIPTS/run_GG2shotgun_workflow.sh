#!/bin/bash

module load snakemake/7.17.1
#module load snakemake
source ~/.bashrc
source ~/.bash_profile



snakemake  -s ./workflow/Snakefile_gg2_shotgun  --cluster \
    " sbatch --account=project_2008347 \
    --partition small --mem 150G --time=3-00:00:00 \
    --ntasks=1 --cpus-per-task=32 \
     -e LOGS/shotgun_err_%A_%a.txt" -j 400 \
     --stats ./comp_jobs_stats.json --keep-going --rerun-incomplete

