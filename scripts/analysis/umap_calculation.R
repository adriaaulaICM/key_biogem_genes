library(tidyverse)
library(umap)
library(parallelDist)

set.seed(42)

ab.raw.filt <- read_rds('data/04_table_gen/abtblalr_filt0rem_8sams.rds') 
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')
annot <- read_rds('data/intermediate_files/annotation.rds')
  
tax <- readRDS('data/04_table_gen/taxonomy_subset.rds')  %>% 
  filter(genename %in% ab.raw.filt$genename) 

lomb <- readRDS('data/intermediate_files/genes_seasonal.rds') 


# Selection and subsetting ------------------------------------------------

genes.wtax <- tax %>% 
  filter(genename %in% lomb$genename) %>% 
  filter(family %in% family.order) %>%
  pull(genename)

ab.alr.long <- ab.raw.filt  %>% 
  filter(genename %in% genes.wtax) 


# Calculation for all the data --------------------------------------------
gene.mat <- ab.raw.filt %>% 
  filter(genename %in% lomb$genename) %>% 
  pivot_wider(names_from = sample, values_from = count, id_cols =  genename) %>% 
  column_to_rownames(var = 'genename') %>% 
  as.matrix()

custom.config <-  umap.defaults
custom.config$n_neighbors <- 20
custom.config$min_dist <- 0.5
gene.umap <- umap(gene.mat)


coord <- gene.umap$layout %>% 
  as_tibble(rownames = 'genename') %>% 
  left_join(annot, by = 'genename')

write_rds(coord, path = 'data/intermediate_files/umap_coord.rds')

# Subsetting for genes with taxa ------------------------------------------

gene.mat.fam <- ab.alr.long %>% 
  pivot_wider(names_from = sample, values_from = count, id_cols =  genename) %>% 
  column_to_rownames(var = 'genename') %>% 
  as.matrix()

gene.umap <- umap(gene.mat.fam)


coord.fam <- gene.umap$layout %>% 
  as_tibble(rownames = 'genename') %>% 
  left_join(annot, by = 'genename')

write_rds(coord.fam, path = 'data/intermediate_files/umap_coord_fam.rds')
