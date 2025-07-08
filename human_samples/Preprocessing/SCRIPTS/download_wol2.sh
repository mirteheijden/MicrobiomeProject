#!/bin/bash
#SBATCH --job-name=download_data           # Name of the job visible in the queue.
#SBATCH --account=project_200XXXX       # Choose the billing project. Has to be defined!
#SBATCH --partition=small          # Job queues: test, interactive, small, large, longrun, hugemem, hugemem_longrun
#SBATCH --time=24:00:00           # Maximum duration of the job. Max: depends of the partition. 
#SBATCH --mem=1G                  # How much RAM is reserved for job per node.
#SBATCH --ntasks=1                # Number of tasks. Max: depends on partition.
#SBATCH --cpus-per-task=1         # How many processors work on one task. Max: Number of CPUs per node.#!/bin/sh

# See
# http://ftp.microbio.me/pub/wol2/
# http://ftp.microbio.me/greengenes_release/
# https://github.com/qiyunzhu/woltka/tree/master/woltka/tests/data/function

wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/taxonomy/00README

# Taxonomy
#wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/taxonomy/taxid.map
#wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/taxonomy/nodes.dmp
#wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/taxonomy/names.dmp
# Use Greengenes2 database for taxonomic annotation
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2/taxonomy/gg2 http://ftp.microbio.me/greengenes_release/2022.10/2022.10.taxonomy.asv.nwk.qza
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2/taxonomy/gg2 http://ftp.microbio.me/greengenes_release/2022.10/2022.10.phylogeny.asv.nwk.qza
# Load the prebuilt bowtie2 database instead of building it from genomes
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/databases/bowtie2/00README
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/databases/bowtie2/WoLr2.1.bt2l
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/databases/bowtie2/WoLr2.2.bt2l
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/databases/bowtie2/WoLr2.3.bt2l
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/databases/bowtie2/WoLr2.4.bt2l
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/databases/bowtie2/WoLr2.rev.1.bt2l
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/databases/bowtie2/WoLr2.rev.2.bt2l
# Function
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/function/uniref/orf-to-uniref.map.xz
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/function/uniref/uniref_name.txt.xz
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/function/go/uniref/process.map.xz
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2 http://ftp.microbio.me/pub/wol2/function/go/go_name.txt
# The function database is missing gene coordinate file
wget -np -nH --cut-dirs=2 -nc -O wol2/function/coords.txt.xz https://raw.githubusercontent.com/qiyunzhu/woltka/master/woltka/tests/data/function/coords.txt.xz
# For the amplicon data mapping, load sequences that are in Greengenes2 tree 
wget -r -np -nH --cut-dirs=2 -nc --directory-prefix=wol2/taxonomy/gg2 http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza
