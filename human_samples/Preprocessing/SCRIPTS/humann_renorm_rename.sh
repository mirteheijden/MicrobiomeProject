#!/bin/bash
#SBATCH --job-name=humann_renorm_rename.sh
#SBATCH --time=10:00:00          
#SBATCH --partition=small        
#SBATCH --nodes=1                
#SBATCH --account=project_200XXXX
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=80G

# This is the script to renorm, rename, and split the function output

echo "start"
# Define the parent dir
parent="/scratch/project_200XXXX/human_feces_batch1"

module load humann/3.8

# Input directory containing the original files
input_dir="$parent/RESULTS/humann3/final/"

# Output directory to store the renamed files
output_dir_split="$parent/RESULTS/humann3/final/splitstrat"

# Find relevant input files and iterate over them
find "$input_dir" -type f -name "genefamilies_Uniref90_*.txt" | while read -r input_file; do
    # Extract the filename without extension
    filename=$(basename "$input_file" .txt)
    
    # Output filenames
    output_renorm="$input_dir/Renorm_${filename}.txt"

    # Run humann_renorm_table command
    humann_renorm_table --input "$input_file" --output "$output_renorm" --units relab --update-snames
    
    #Run humann_split_stratified_table command on renorm output
    humann_split_stratified_table --input "$output_renorm" --output "$output_dir_split"


    # Optionally, check if the command was successful
    if [ $? -eq 0 ]; then
        echo "File $input_file processed successfully."
    else
        echo "Error processing file $input_file."
    fi
done


#rename the renorm table and split
# uniref90-ko
rename_KO() {
    input="$parent/RESULTS/humann3/final/Renorm_genefamilies_Uniref90_KO.txt"
    output="$parent/RESULTS/humann3/final/RenormRename_genefamilies_Uniref90_KO.txt"
    DB="$parent/DATABASE_DIR/utility_mapping/map_ko_name.txt.gz"
    humann_rename_table --input "$input" --output "$output" -c $DB
}

# uniref90-EggNOG
rename_eggnog() {
    input="$parent/RESULTS/humann3/final/Renorm_genefamilies_Uniref90_EggNOG.txt"
    output="$parent/RESULTS/humann3/final/RenormRename_genefamilies_Uniref90_EggNOG.txt"
    DB="$parent/DATABASE_DIR/utility_mapping/map_eggnog_name.txt.gz"
    humann_rename_table --input "$input" --output "$output" -c $DB
}

# uniref90-GO
rename_GO() {
	input="$parent/RESULTS/humann3/final/Renorm_genefamilies_Uniref90_GO.txt"
    output="$parent/RESULTS/humann3/final/RenormRename_genefamilies_Uniref90_GO.txt"
    DB="$parent/DATABASE_DIR/utility_mapping/map_go_name.txt.gz"
    humann_rename_table --input "$input" --output "$output" -c $DB
}

# uniref90-Pfam
rename_pfam() {
	input="$parent/RESULTS/humann3/final/Renorm_genefamilies_Uniref90_Pfam.txt"
    output="$parent/RESULTS/humann3/final/RenormRename_genefamilies_Uniref90_Pfam.txt"
    DB="$parent/DATABASE_DIR/utility_mapping/map_pfam_name.txt.gz"
    humann_rename_table --input "$input" --output "$output" -c $DB
}

# uniref90-EC
rename_ec() {
    input="$parent/RESULTS/humann3/final/Renorm_genefamilies_Uniref90_EC.txt"
    output="$parent/RESULTS/humann3/final/RenormRename_genefamilies_Uniref90_EC.txt"
    DB="$parent/DATABASE_DIR/utility_mapping/map_ec_name.txt.gz"
    humann_rename_table --input "$input" --output "$output" -c $DB
}


# uniref90-MetaCyc
rename_metacyc() {
	input="$parent/RESULTS/humann3/final/Renorm_genefamilies_Uniref90_MetaCyc.txt"
    output="$parent/RESULTS/humann3/final/RenormRename_genefamilies_Uniref90_MetaCyc.txt"
    DB="$parent/DATABASE_DIR/utility_mapping/map_level4ec_uniref90.txt.gz"
    humann_rename_table --input "$input" --output "$output" -c $DB --names metacyc-rxn
}

rename_KO
echo "done rename_KO"
rename_eggnog
echo "done rename_eggnog"
rename_GO
echo "done rename_GO"
rename_pfam
echo "done rename_pfam"
rename_ec
echo "done rename_EC"
rename_metacyc
echo "done rename_metacyc" 


#split rename table
# Find relevant input files and iterate over them
find "$input_dir" -type f -name "RenormRename_genefamilies_Uniref90_*.txt" | while read -r input_file; do
    # Extract the filename without extension
    filename=$(basename "$input_file" .txt)
    
    #Run humann_split_stratified_table command on renorm output
    humann_split_stratified_table --input "$input_file" --output "$output_dir_split"


    # Optionally, check if the command was successful
    if [ $? -eq 0 ]; then
        echo "File $input_file processed successfully."
    else
        echo "Error processing file $input_file."
    fi
done

