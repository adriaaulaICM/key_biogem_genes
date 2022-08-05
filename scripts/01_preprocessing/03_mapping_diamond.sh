#!/bin/bash

#######################
#SBATCH --job-name=famli2
#SBATCH --account=emm3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=150G
#SBATCH --array=16,17,20,21,33,65,66,67,70,76,77,78,81,83,84%5
#SBATCH --output=data/logs/famli2/famli_mapping_%A_%a.out
#SBATCH --error=data/logs/famli2/famli_mapping_%A_%a.err
#######################

module load diamond
module load famli

#beforehand, diamonddb build
#diamond makedb --in data/annotation/putative_genes/final_min250_cogg_corr.faa -d data/annotation/putative_genes/genes_db.dmnd

CPU=20
# Files
READS_DIR='data/raw/clean.reads'
SAMPLE=$(ls ${READS_DIR}/*fastq.gz | awk -F'/' '{print $NF}' | awk -F'_' '{print $1}' | sort -u | awk "NR == ${SLURM_ARRAY_TASK_ID}")
DIAMONDDB='data/annotation/putative_genes/genes_db.dmnd'

# Check 
echo $SAMPLE 

# Output redirection
OUT_DIR='/mnt/lustre/scratch/aauladell/biogem_7years_bigoutputs/mapping_diamond/'
FAM_DIR='/mnt/lustre/scratch/aauladell/biogem_7years_bigoutputs/mapping_parse_famli/'

#cat ${READS_DIR}/${SAMPLE}_R1.clean.fastq.gz ${READS_DIR}/${SAMPLE}_R2.clean.fastq.gz > ~/scratch/samswhole/${SAMPLE}.fastq.gz

#diamond blastx --query ~/scratch/samswhole/${SAMPLE}.fastq.gz \
#               --out ${OUT_DIR}/${SAMPLE}.aln \
#               --min-score 20 \
#               --threads ${CPU} --db ${DIAMONDDB} \
#               --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
#               --query-cover 90 \
#               --id 95 --top 5 --block-size 9 \
#               --query-gencode 11 --unal 0

famli filter --input ${OUT_DIR}/${SAMPLE}.aln \
             --output-aln ${FAM_DIR}/${SAMPLE}.aln \
             --threads ${CPU} 

