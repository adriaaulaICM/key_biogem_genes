library(tidyverse)

source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')

# we want the total relative abundance for each gene variant inside each gene
# this is helpful to quantify how abundant is one of the variants
# que cosas eh? 

# Data import -------------------------------------------------------------

ab.raw.filt <- read_rds('data/04_table_gen/abtblraw_all.rds')
annot <- read_rds('data/intermediate_files/annotation.rds')
  
rank.df <- ab.raw.filt %>% 
  group_by(genename) %>% 
  summarize(count = sum(count)) %>% 
  left_join(annot %>% select(genename, gene_name),
            by = 'genename') %>% 
  group_by(gene_name) %>% 
  mutate(total.gene = sum(count),
         relab.total = count / total.gene) %>% 
  arrange(gene_name, -count)  %>% 
  mutate(cumsum = cumsum(relab.total))
  
saveRDS(rank.df, 'data/intermediate_files/total_relab_genev.rds')
