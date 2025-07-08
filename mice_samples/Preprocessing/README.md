# mice_feces_IV2
Adapted from @erawijantari

Puhti - pre-processing

Set working directory: config, DATA, DB, LOGS, RESULTS, SCRIPTS, workflow, DATABASE_DIR

download databases:

Mice

download FASTA file from here: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964188535.1/ to local dir. 
Go to NCBI_dataset --> data --> GCA_964188535.1. Copy file .fna file to DB map, change file.fna to mice.FASTA and copy to Puhti DB using winSCP
$ --hostremoval_reference /scratch/project_200XXXX/mice_feces_IV2/DB/mice.fasta
Let use --save_hostremoval_unmapped just to be safe to store the unmapped host reads (typically this is the final processed reads) ???
Human:

download FASTA file from here: https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_009914755.4/
If you downloaded it from local dir copy it to puhti using scp
$ unzip T2T-CHM13v2.0.zip
$ --hostremoval_reference /scratch/project_200XXX/nextflow_metagenomics/DB/human_CHM13
Let use --save_hostremoval_unmapped just to be safe to store the unmapped host reads (typically this is the final processed reads) ???

MetaPhlan4:

Follow set-up here: https://docs.csc.fi/apps/metaphlan/
Greengenes2: adapted from @TuomasBorman Setting and Running GG2 plugin - FOR SHOTGUN

For shotgun reads, we need newest QIIME2 tools with GG2 plugin and Woltka

Go to directory where you want to save the config for installation (e.g in project or home, not scratch)

Download and save the environment config file needed for installation wget https://raw.githubusercontent.com/qiime2/distributions/dev/2023.9/shotgun/released/qiime2-shotgun-ubuntu-latest- conda.yml

Install QIIME2:

$ module load tykky
$ mkdir qiime2-shotgun-2023.09
$ conda-containerize new --mamba --prefix qiime2-shotgun-2023.09/ qiime2-shotgun-ubuntu-latest-conda.yml
Create a file for plugin installation

create a text file named post_install_plugins_shotgun.txtcontaining the following info: (you can use nano, vim, or other favorite text editor)

$ pip install q2-greengenes2 pip install woltka conda install -c bioconda bowtie2 pip install https://github.com/knights-lab/SHOGUN/archive/master.zip pip install https://github.com/qiime2/q2-shogun/archive/master.zip qiime dev refresh-cache

Install plugins

$ conda-containerize update qiime2-shotgun-2023.09/ --post-install post_install_plugins_shotgun.txt
Add the software path so that the software is executable

$ export PATH="qiime2-shotgun-2023.09/bin:$PATH" (preferably using full path, here is the example of current directory)

Download the WoL2 database: http://ftp.microbio.me/pub/wol2/ https://github.com/qiyunzhu/woltka/blob/master/doc/wol.md#the-wol-database

Let's put this in the DB directory, remember to change the project number in the script to match yours

$ cd ./DB
$ chmod +x download_wol2.sh
$ sbatch ../SCRIPTS/download_wol2.sh


NEXTFLOW

Cretae the database sheet: "config/database.csv"

Create in Excel, save as .CSV, move to Putty. Database is only real database (MetaPhlan4, mOTU), not the fasta.file documents (mice, human_CHM13)
Crease the samplesheet: "config/samplesheet.csv"

Run part one of taxprofiler_part1.sh script for cleaning data + removing mice contamination (host contamination):

$ sbatch SCRIPTS/taxprofiler_part1.sh
When done, create a second samplesheet: "config/samplesheet2.csv" and run the second part of taxprofiler for removing human contamination (during sample processing):

$ sbatch SCRIPTS/taxprofiler_part2.sh
Be sure to update the pipeline regularly:

$ nextflow pull nf-core/taxprofiler
MetaPhlan4 results "metaphlan_db_meta4_combined_reports.txt" can be found RESULTS/Metaphlan. Use WinSCP to transfer to harddrive or to Rstudio. If they are not merged yet, merge the .metaphlan_profile.txt files yourself:

$ module load metaphlan

$ cd /scratch/project_200XXXX/mice_feces_IV2/RESULTS/metaphlan

$ merge_metaphlan_tables.py db_meta4/*metaphlan_profile.txt > metaphlan_db_meta4_combined_reports.txt

Check the removal of the mice and human reads by:

Calculating the reads in the Raw data (reads_count_preprocess), after mice removal () and after human removal ():
$ cd ./SCRIPTS
Create the scripts using nano ( $ nano reads_count_preprocess, $ nano reads_count_mice_removal_postprocess, nano reads_count_mice_and_human_removal_postprocess)
Run the scripts, one by one:
$ chmod +x reads_count_preprocess
$ ./reads_count_preprocess
$ chmod +x reads_count_mice_removal_postprocess
$ ./reads_count_mice_removal_postprocess
$ chmod +x reads_count_mice_and_human_removal_postprocess
$ ./reads_count_mice_and_human_removal_postprocess
Results can be found in RESULTS/reads_count_preprocess.tsv, RESULTS/reads_count_mice_removal_postprocess and RESULTS/reads_count_mice_and_human_removal_postprocess
Move .tsv files to PC using WinSCP

SNAKEMAKE
Snakemake (Greengenes2) Adapted from @erawijantari and @TuomasBorman
Create the "Snakefile_gg2_shotgun" under the "workflow" directory. Adjust the path to the qiime config, project number, and files locations according to your need. We will use the pre- processed reads from nextflow that has been stored in RESULTS/analysis_ready_fastqs

Create the "run_GG2shotgun_workflow.sh" bash script in SCRIPTS to execute snakemake Please see SCRIPTS/run_GG2shotgun_workflow.sh, modify the project number

Exectue run:

$ chmod +x ./SCRIPTS/run_GG2shotgun_workflow.sh
$ ./SCRIPTS/run_GG2shotgun_workflow.sh
Results "counts.qza" and "taxonomy.qza" can be found at RESULTS/gg2. Use WinSCP to transfer to harddrive or to Rstudio

HUMANN3 (functional profiling) 
adapted from @erawijantari, @TuomasBorman and @KatariinaParnanen
Create the "Snakefile_humann" using nano under the "workflow" directory. Remember to adjust the path, project allocation, and file locations.

Create a bash script "SCRIPTS/run_humann_workflow.sh" to execute snakemake.

Execute run:

$ chmod +x ./SCRIPTS/run_humann_workflow.sh
$ ./SCRIPTS/run_humann_workflow.sh
IMPORTANT NOTES: always remember to check the tools version especially if running in batch. For running samples more than 100, it is more convenience to use the group function of snakemake.

When done, continue with merging genefamilies, pathabundance and pathcoverage files:

$ module load humann/3.8
$ cd /scratch/project_2008347/mice_feces_IV2/RESULTS/humann3
$ mkdir merged
$ humann_join_tables --input /scratch/project_2008347/mice_feces_IV2/RESULTS/humann3/raw --file_name genefamilies.tsv --output /scratch/project_2008347/mice_feces_IV2/RESULTS/humann3/merged/genefamilies.txt
$ humann_join_tables --input /scratch/project_2008347/mice_feces_IV2/RESULTS/humann3/raw --file_name pathabundance.tsv --output /scratch/project_2008347/mice_feces_IV2/RESULTS/humann3/merged/pathabundance.txt
$ humann_join_tables --input /scratch/project_2008347/mice_feces_IV2/RESULTS/humann3/raw --file_name pathcoverage.tsv --output /scratch/project_2008347/mice_feces_IV2/RESULTS/humann3/merged/pathcoverage.txt
Look if results are presnet in RESULTS/humann3/merged

$ cd /scratch/project_2008347/mice_feces_IV2/RESULTS/humann3/merged

Create humann_regroup_funct.sh and humann_renorm_rename.sh script and put in the SCRIPTS directory (see biobakery for info, https://github.com/biobakery/humann?tab=readme-ov-file#guides-to-humann-utility-scripts) 

$ humann_databases --download utility_mapping full /scratch/project_2008347/mice_feces_IV2/DATABASE_DIR
Run humann_regroup_funct.sh script
$ bash ./SCRIPTS/humann_regroup_funct.sh

Run humann_renorm_rename.sh script
$ bash ./SCRIPTS/humann_renorm_rename.sh

Output files can be found in RESULTS/humann3/final/splitstrat

