#!/bin/bash

#######################
#SBATCH --job-name=parsing
#SBATCH --account=emm3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --error=data/logs/blast_%J.err
#SBATCH --output=data/logs/blast_%J.out

#######################

# We take all the cogs specified in the gene database and search for them 
# we only keep the ones presenting a length X 

COG_CHECK="/mnt/lustre/scratch/aauladell/biogem_7years_bigoutputs/cog_annotation/BBMO.v2_95id.faa.cog.annot.description.txt"

grep -f <(cut -f14 results/summary_charact/database_genes.txt | sort -u | grep 'COG') ${COG_CHECK} | \
    awk '{if (($6 - $5) > 83) print $0}' \
    > data/annotation/cog_annot/putative_function.txt 

