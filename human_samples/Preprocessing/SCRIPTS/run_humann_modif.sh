#!/bin/bash
# SBATCH --job-name=humann_modif 
# SBATCH --account= project_200XXXX
# SBATCH --time=10:00:00
# SBATCH --mem-per-cpu=80G 
# SBATCH --partition=small
# SBATCH -e LOGS/human_err_%A_%a.txt -o LOGS/human_out_%A_%a.txt

source ~/.bashrc
source ~/.bash_profile
module load humann/3.8

#create LOGS folder to make structured output
#mkdir -p ./LOGS

bash ./SCRIPTS/humann_regroup_funct.sh 
bash ./SCRIPTS/humann_renorm_rename.sh
