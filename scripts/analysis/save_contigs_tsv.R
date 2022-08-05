library(tidyverse)
ab.alr <- read_rds('data/04_table_gen/abtblraw_all.rds')
contigs.df <- ab.alr %>%
  select(genename) %>% 
  mutate(contigs = str_extract(genename, 'BL[0-9]*\\.k127_[0-9]*')) %>%
  separate(col = contigs,
           into = c('sam_origin', 'contig_id'),
           sep = '\\.')

contigs.df %>%
  select(sam_origin, contig_id) %>%
  distinct() %>%
  write_tsv('data/annotation/putative_genes/contigs_w_putative_genes.tsv')


contigs.df %>% 
  select(sam_origin, contig_id) %>%
  distinct() %>%
  group_by(sam_origin) %>% 
  nest() %>% 
  rowwise() %>% 
  pwalk(~.x %>% write_tsv(x = .x %>% select(data),
                   path = 'data/annotation/putative_genes/',
                   ample, '.tsv'))
  


contigs.df %>% 
  select(sam_origin, contig_id) %>%
  head(n = 1500) %>% 
  group_by(sam_origin) %>% 
  nest() %>% 
  pwalk(~write_tsv(list(...)$data, paste0(list(...)$sam_origin, ".tsv")))

