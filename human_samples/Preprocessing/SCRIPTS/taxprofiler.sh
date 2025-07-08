#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --partition=small
#SBATCH --account=project_200XXXX
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=80G

export SINGULARITY_TMPDIR=$PWD
export SINGULARITY_CACHEDIR=$PWD
unset XDG_RUNTIME_DIR

# Activate  Nextflow on Puhti
module load nextflow

#run the nextflow
#for testing (in profile params use test,singularity flag )

nextflow run nf-core/taxprofiler -r 1.1.5 -resume \
   -profile singularity --max_cpus 4 \
   --input ./config/samplesheet.csv \
   --databases ./config/database.csv \
   --outdir ./RESULTS  \
   --perform_shortread_qc --save_preprocessed_reads\
   --perform_shortread_complexityfilter  --shortread_complexityfilter_tool bbduk \
   --perform_shortread_hostremoval \
   --hostremoval_reference /scratch/project_200XXXX/nextflow_metagenomics/DB/T2T-CHM13v2.0.zip \
   --shortread_hostremoval_index /scratch/project_200XXXX/nextflow_metagenomics/DB/human_CHM13/ \
   --save_hostremoval_unmapped \
   --save_analysis_ready_fastqs \
   --run_profile_standardisation \
   --run_motus --run_metaphlan --run_kraken2

