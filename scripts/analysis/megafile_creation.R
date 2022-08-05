# we want a file with all the information together to generate good and/or fast
# vizualizations of parts of the dataset.    

library(tidyverse)

ab.raw.filt <- read_rds('data/04_table_gen/abtblalr_filt0rem_8sams.rds') 
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')
annot <- read_rds('data/intermediate_files/annotation.rds')
  
tax <- readRDS('data/04_table_gen/taxonomy_contigwise.rds')  %>% 
  filter(genename %in% ab.raw.filt$genename) 

lomb <- readRDS('data/intermediate_files/genes_seasonal.rds') 


megafile <- ab.raw.filt %>% 
  left_join(annot, by = 'genename') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name'))  %>% 
  left_join(tax, by = 'genename') %>% 
  left_join(lomb, by = 'genename')

small <- megafile %>% 
  select(genename, sample, count, gene_name, day_of_year, season,
         species, genus, family, order, class, peak)

write_rds(small, path = 'data/intermediate_files/megafile_8samples.rds')
