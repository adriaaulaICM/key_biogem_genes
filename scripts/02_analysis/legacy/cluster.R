library(tidyverse)

ab.alr.filt <- read_rds('data/04_table_gen/abtblalr_filt0rem_25sams.rds') 

envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')

lomb <- readRDS('data/intermediate_files/genes_seasonal.rds')

ab.alr <- ab.alr.filt %>% 
  filter(genename %in% lomb$genename)


themat <- ab.alr %>% 
  pivot_wider(names_from = sample, values_from = count) %>% 
  column_to_rownames( var = 'genename') %>% 
  as.matrix()

sim <- themat / sqrt(rowSums(themat * themat))
sim <- sim %*% t(sim)

D_sim <- as.dist(1 - sim)

clusts <- hclust(D_sim)

clusters <- kmeans(D_sim, centers = 12)
