library(rhmmer) 
library(tidyverse)

# This script checks the results obtained for an specific KEGG (the one from dmdA) and 
# compares the matches against the curated HMM from José Gonzalez et al. 2019 

tables <- list.files('~/icm/tunnel/hmm_tables',full.names = T)
dmdajose <- rhmmer::read_tblout('~/icm/tunnel/resdmda.txt')
res <- tables %>% 
  map_df(~rhmmer::read_tblout(.x))

dmda <- res %>% filter(query_name == 'K17486') %>% 
  mutate( jose = ifelse( domain_name %in% dmdajose$domain_name, 'yes', 'no'))


ggplot(dmda, aes(sequence_score)) + 
  geom_histogram(bins = 60, aes(fill = jose)) + 
  geom_vline(xintercept =546.2 )


bind_rows(dmda, dmdajose, .id = 'casa') %>% 
ggplot(aes(sequence_score)) + 
  geom_histogram(bins = 60, aes(fill = jose)) + 
  geom_vline(xintercept =546.2 ) + 
  facet_wrap(~casa, ncol = 1)

res %>% filter(query_name == 'K17486') %>% 
  View()
  filter(sequence_score > 527.20 *1.1)
