library(tidyverse)


tax.contig <- readRDS('data/04_table_gen/taxonomy_contigwise.rds')
tax <- readRDS('data/04_table_gen/taxonomy_subset.rds')


na.phylum <- tax.contig %>% 
  as_tibble() %>% 
  filter(is.na(phylum)) %>%
  # filter(is.na(last_rank)) %>% 
  pull(genename) 

ranks.from.sgene <- tax %>% 
  filter(genename %in% na.phylum) 


taxonomy.mostcomplete <- tax.contig %>% 
  filter(!genename %in% na.phylum) %>% 
  bind_rows(ranks.from.sgene)

write_rds(x = taxonomy.mostcomplete, 
          file = 'data/04_table_gen/taxonomy_merged_contig_sgene_gtdb.rds')
