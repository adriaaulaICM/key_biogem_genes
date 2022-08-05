library(tidyverse)

source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')


sag.groups <- readxl::read_xlsx('data/04_taxgenes/emi14896-sup-0002-tables1.xlsx',
                                skip = 2) %>% 
  select(Genome, Subclade, Genomospecies) %>% 
  rename( sag = Genome)

sag.tax.df <- read_tsv('data/04_taxgenes/pelagibacter_sags_tax.aln',
                       col_names =  c("qseqid", "sseqid", "pident", "length",
                                      "mismatch", "gapopen", "qstart", "qend",
                                      "sstart", "send", "evalue", "bitscore",
                                      "qlen", "slen")) %>% 
  rename( genename = qseqid)  %>% 
  mutate(cov = length / qlen) %>% 
  separate(col = sseqid, into = c('geneid', 'sag'), sep = '_') %>% 
  filter(cov >= 0.8) %>% 
  distinct()

sag.tax.df %>% 
  left_join(sag.groups, by = 'sag') %>% 
  select(genename, sag, pident, cov, Subclade, Genomospecies)  %>% 
  write_tsv('data/04_taxgenes/taxonomy_genevars_sags_haro.tsv')
