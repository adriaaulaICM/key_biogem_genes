#!/bin/bash

#######################
#SBATCH --job-name=agg
#SBATCH --account=emm3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --error=data/logs/agg_%J.err
#SBATCH --output=data/logs/agg_%J.out
#######################

module load R 

# Fasta sequences
whole_gc_dna='data/whole_catalog_bbmo_nuc/whole_catalog_bbmo_nuc'
whole_gc_aa='/mnt/lustre/bio/shared/blgc/gc/BBMO.v2_95id_min250.faa'


# kofam, cog and PR  putative genes annotation
kofam='data/annotation/kofam_search/putative_function.txt'
cog='data/annotation/cog_annot/putative_function.txt'
pr='data/annotation/PR_search/putative_pr.txt'

#out path 
out_path='data/annotation/putative_genes'

# putting everything together
temp="${out_path}/sel_genes.txt"

awk '{print $1}' ${kofam} > ${temp}
awk '{print $1}' ${cog} >> ${temp}
awk '{print $1}' ${pr} >> ${temp}

# Selecting the genes from the whole catalog 
seqkit grep -n -j 24 -f ${temp} ${whole_gc_dna}  | seqkit seq -m 250 > ${out_path}/final_min250.fna

# i keep the names
seqkit seq -n ${out_path}/final_min250.fna > ${out_path}/names_seq_min250.txt 

seqkit grep -n -j 24 -f ${out_path}/names_seq_min250.txt ${whole_gc_aa} > ${out_path}/final_min250.faa

# EXTRA: cleaning COG sequences
Rscript scripts/helper_subset_cog_sequences.sh ${out_path}/final_min250.fna  \
    ${out_path}/final_min250.faa \
    $cog \
    ${out_path}/final_min250_cogg_corr



