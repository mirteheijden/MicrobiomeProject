#!/bin/bash

#SBATCH --job-name=humann_regroup_funct.sh
#SBATCH --time=10:00:00          #adjust for the time needed
#SBATCH --partition=small        #for multinode, large partition should be use
#SBATCH --nodes=1                #max nodes allowed 26
#SBATCH --account=project_200XXXX
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=80G

 
#This is the script to regroup the Uniref 90 into various function database Due to small sample size, the process will be done in 
# whole table if samples more than 1000, consider doing the step for each sample individually and merge after due to memmory requirement issue.

echo "start"

# Define the parent dir
parent="/scratch/project_200XXXX/human_feces_batch1"
module load humann/3.8
input="$parent/RESULTS/humann3/merged/genefamilies.txt"
    
# uniref90-ko
regroup_KO() { 
	output="$parent/RESULTS/humann3/final/genefamilies_Uniref90_KO.txt" 
    	DB="$parent/DATABASE_DIR/utility_mapping/map_ko_uniref90.txt.gz" 
	humann_regroup_table --input "$input" --output "$output" -c $DB
}

# uniref90-EggNOG
regroup_eggnog() { 
	output="$parent/RESULTS/humann3/final/genefamilies_Uniref90_EggNOG.txt" 
    	DB="$parent/DATABASE_DIR/utility_mapping/map_eggnog_uniref90.txt.gz" 
	humann_regroup_table --input "$input" --output "$output" -c $DB
}
# uniref90-GO
regroup_GO() { 
	output="$parent/RESULTS/humann3/final/genefamilies_Uniref90_GO.txt" 
    	DB="$parent/DATABASE_DIR/utility_mapping/map_go_uniref90.txt.gz" 	humann_regroup_table --input "$input" --output "$output" -c $DB
}

# uniref90-Pfam
regroup_pfam() { 
	output="$parent/RESULTS/humann3/final/genefamilies_Uniref90_Pfam.txt" 
    	DB="$parent/DATABASE_DIR/utility_mapping/map_pfam_uniref90.txt.gz" 
	humann_regroup_table --input "$input" --output "$output" -c $DB
}
# uniref90-EC
regroup_ec() { 
	output="$parent/RESULTS/humann3/final/genefamilies_Uniref90_EC.txt" 
    	DB="$parent/DATABASE_DIR/utility_mapping/map_level4ec_uniref90.txt.gz" 
	humann_regroup_table --input "$input" --output "$output" -c $DB
}
# uniref90-MetaCyc
regroup_metacyc() { 
	output="$parent/RESULTS/humann3/final/genefamilies_Uniref90_MetaCyc.txt" 
    	DB="$parent/DATABASE_DIR/utility_mapping/map_level4ec_uniref90.txt.gz" 
	humann_regroup_table --input "$input" --output "$output" -c $DB -g uniref90_rxn
}


regroup_KO 
echo "done regroup_KO" 
regroup_eggnog 
echo "done regroup_eggnog" 
regroup_GO 
echo "done regroup_GO" 
regroup_pfam 
echo "done regroup_pfam" 
regroup_ec 
echo "done regroup_EC" 
regroup_metacyc 
echo "done regroup_metacyc"
