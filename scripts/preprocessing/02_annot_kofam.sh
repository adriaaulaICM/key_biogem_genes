#!/bin/bash

#######################
#SBATCH --job-name=hmm
#SBATCH --account=emm3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --error=data/logs/hmm_%J.err
#SBATCH --output=data/logs/hmm_%J.out
#######################

#### This script takes a HMM from KOFAM list, 
#### searches them in the KOFAM DB and performs a search in the 
#### BBMO gene catalog 

hmm_db='data/databases/hmm_gene_selection'
out_kofam="data/annotation/kofam_search"

# Programs
module load hmmer

# remove initial selection (in case i change it)
rm ${hmm_db}/*
# Retrieve HMM
for kegg in $(cut -f14 results/summary_charact/database_genes.txt | sort -u);
    do cp ~/scratch/databases/kofam_profiles/profiles_bact/${kegg}.hmm ${hmm_db}/. ;
    done


hmm_prefix=$(ls -1 ${hmm_db} | sed 's/.hmm//g')

# Do the search with a threshold specified in the KO_list file 
for hmm in ${hmm_prefix};
    do hmmsearch -T $( awk "/$hmm/{print $1}" ~/scratch/databases/kofam_profiles/ko_list | awk '{print $2}') \
        --cpu 40 \
        --tblout ${out_kofam}/hmm_tables/${hmm}_results.tsv \
        ${hmm_db}/${hmm}.hmm  \
        data/whole_catalog_bbmo_aa/BBMO.v2_95id_clust_cleaned.faa ;
    done

# Save results 
awk '/BL/{print $1,$3}' ${out_kofam}/hmm_tables/K*  > ${out_kofam}/putative_function.txt

