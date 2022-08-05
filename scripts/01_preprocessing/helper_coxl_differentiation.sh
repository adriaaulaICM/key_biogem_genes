#!/bin/bash

# this script takes all the coxL and differentiates both of them as CODHI or CODHII

seqkit grep -f (grep 'K03520' data/annotation/kofam_search/putative_function.txt | awk {'print $1'}) data/annotation/putative_genes/final_min250.faa > /tmp/codh.fasta

seqkit grep -s -r -p  'AY.CSFR' /tmp/codh.fasta | seqkit seq -n  | awk  '{print $0, "CODHI"}' OFS='\t' > data/annotation/kofam_search/codh1_class.tsv
