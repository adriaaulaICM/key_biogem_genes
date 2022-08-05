#!/bin/bash

#######################
#SBATCH --job-name=pr
#SBATCH --account=emm3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --error=data/logs/pr_%J.err
#SBATCH --output=data/logs/pr_%J.out
#######################

# By using MicRhoDE as a Proteorhodopsin DB, a search is performed in the database 
# with diamond 

module load diamond

PRdb='data/databases/PR_db'
resdir='data/annotation/PR_search'
scratch='/mnt/lustre/scratch/aauladell'

diamond makedb --in ${PRdb}/PR_aa_only.faa -d ${PRdb}/PR_diamond_db.dmnd

diamond blastp -d data/whole_catalog_bbmo_aa/wholeaa_diamond_db.dmnd \
               -q ${PRdb}/PR_aa_only.faa \
               --id 70 \
               --query-cover 80 \
               --subject-cover 80 \
               --threads 48 \
               -o ${resdir}/pr_search_diamond.tsv

cut -f2 ${resdir}/pr_search_diamond.tsv | sort -u | awk '{print $0,"PR"}' > ${resdir}/putative_pr.txt

cut -f1,3,7,11,13,15,16 --output-delimiter=$'\t' ${PRdb}/MicRhoDE_051214.txt  | 
    grep -f <(cut -f1 ${resdir}/pr_search_diamond.tsv | sort -u | sed 's/_1$//g') > ${resdir}/tax_info.tsv

# Putting the names to the cool tax doc 
head=$(cut -f1,3,7,11,13,15,16 --output-delimiter=$'\t'  ${PRdb}/MicRhoDE_051214.txt | head -n1)
(echo ${head} && cat ${resdir}/tax_info.tsv ) > ${resdir}/tax_info2.tsv && mv ${resdir}/tax_info2.tsv ${resdir}/tax_info.tsv
 
