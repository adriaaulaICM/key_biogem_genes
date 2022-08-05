library(tidyverse)
library(DT)


ab.raw.filt <- read_rds('data/04_table_gen/abtblraw_filt0rem_8sams.rds') 

annot <- read_rds('data/intermediate_files/annotation.rds') %>% 
  select(genename, cycle, gene_name)

vars <- ab.raw.filt %>% 
  left_join(annot, by = 'genename') %>% 
  group_by(gene_name, sample) %>% 
  mutate(relab = count / sum(count))  %>% 
  group_by(gene_name, genename) %>% 
  summarise(m.relab = mean(relab)) %>% 
  group_by(gene_name) %>% 
  arrange(gene_name, -m.relab) %>% 
  top_n(n = 5, wt = m.relab) %>% 
  pull(genename)

tax <- readRDS('data/04_table_gen/taxonomy_subset.rds') %>% 
  filter(genename %in% vars)


tax.uniref <- read_tsv('data/04_taxgenes/unirefResult.tsv',
                col_names = c('genename', 'genomeid', 'last_rank',
                              'rank_name', 'taxonomy')) %>% 
  filter(genename %in% vars) %>% 
  separate(col = taxonomy,
           into = c("species","genus","family",
                    "order","class","phylum","superkingdom"), sep = ";") %>% 
  mutate( label = case_when(is.na(phylum) ~ rank_name,
                            !is.na(phylum) ~ str_c(phylum, family,sep = ", "),
                            TRUE ~ rank_name)) %>% 
  distinct(genename, .keep_all = T)


interactive.tab <- datatable(tax)
saveWidget(interactive.tab, 'results/tables/taxonomy_gene_gtdb.html')

interactive.tab <- datatable(tax.uniref)
saveWidget(interactive.tab, 'results/tables/taxonomy_gene_uniref.html')
