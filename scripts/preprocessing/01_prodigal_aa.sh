#!/bin/bash

#######################
#SBATCH --job-name=prod
#SBATCH --account=emm3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --error=data/logs/prodigal_%J.err
#SBATCH --output=data/logs/prodigal_%J.out
#######################

# Programs
module load prodigal

prodigal -a data/aa_db/BBMO.v2_95id_clust.faa -i data/raw/gc/BBMO.v2_95id_clust.fasta

# We change the name since the one by default from prodigal can put us in trouble 
sed 's/_[0-9] # .*//g' data/aa_db/BBMO.v2_95id_clust.faa >> data/aa_db/BBMO.v2_95id_clust_cleaned.faa
