#!/bin/bash

#######################
#SBATCH --job-name=tax
#SBATCH --account=emm3
#SBATCH --ntasks=1
#SBATCH --mem=250G
#SBATCH --cpus-per-task=48
#SBATCH --constraint=supermicro
#SBATCH --error=data/logs/tax_%J.err
#SBATCH --output=data/logs/tax_%J.out
#SBATCH --exclude=c[11-22]

#######################

dbdir='/mnt/lustre/scratch/aauladell/databases/gtdb-r89_54k/'
unidir='/mnt/lustre/scratch/aauladell/databases/unirefdb/'
mmseqs='scripts/mmseqs/mmseqs/bin/mmseqs'

# Only once 
#${mmseqs} createdb ${dbdir}/gtdb-r89_54k.rep_seq.fasta.gz ${dbdir}/gtdb_54kDB
#${mmseqs} createtaxdb ${dbdir}/gtdb_54kDB ${dbdir}/tmp  --ncbi-tax-dump ${dbdir}/taxonomy/ --tax-mapping-file ${dbdir}/gtdb_taxidmapping 

#we also create the db from the aa genes 
#./${mmseqs} createdb data/annotation/putative_genes/final_min250.faa data/annotation/putative_genes/aa_mmseq_db

# the taxonomy
./${mmseqs} taxonomy --search-type 1 --threads 48 --lca-ranks "species;genus;family;order;class;phylum;superkingdom" \
        data/annotation/putative_genes/aa_mmseq_db  \
        ${dbdir}/gtdb_54kDB data/04_taxgenes/taxresdb ${dbdir}/tmp/
#
## extract the results 
./${mmseqs} createtsv data/annotation/putative_genes/aa_mmseq_db data/04_taxgenes/taxresdb data/04_taxgenes/taxonomyResult.tsv

# we will also generate a taxonomy from unirefDB
./${mmseqs} taxonomy --search-type 1 --threads 48 --lca-ranks "species;genus;family;order;class;phylum;superkingdom" \
        data/annotation/putative_genes/aa_mmseq_db  \
        ${unidir}/unirefDB data/04_taxgenes/unirefdb ${unidir}/tmp/

# extract the results 
./${mmseqs} createtsv data/annotation/putative_genes/aa_mmseq_db data/04_taxgenes/unirefdb data/04_taxgenes/unirefResult.tsv
