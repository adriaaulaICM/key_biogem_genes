library(tidyverse)
library(patchwork)
library(umap)
library(speedyseq)
library(parallelDist)

set.seed(42)

source('../202001_ts_arizona_chap1/src/util/params-graphs.R')

# Data import -------------------------------------------------------------


ab.raw.filt <- read_rds('data/04_table_gen/abtblalr_filt0rem_8sams.rds') 
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')
annot <- read_rds('data/intermediate_files/annotation.rds')

bl.phy <- read_rds('../202001_ts_arizona_chap1/data/cleaned/blphyloseq.rds')

# UMAP genes --------------------------------------------------------------

gene.mat <- ab.raw.filt %>% 
  pivot_wider(names_from = sample,
              values_from = count,
              id_cols =  genename,
              values_fill = 0) %>% 
  column_to_rownames(var = 'genename') %>% 
  as.matrix() %>% 
  t()

custom.config <-  umap.defaults
custom.config$n_neighbors <- 20
custom.config$min_dist <- 0.5
gene.umap <- umap(gene.mat)


# UMAP 16S  ---------------------------------------------------------------
otu.mat <- subset_samples(bl.phy, year >= 2009) %>%
  otu_table() %>% 
  as.matrix()

otu.umap <- umap(otu.mat)


gene.coord <- gene.umap$layout %>% 
  as_tibble(rownames = 'Sample_name') %>% 
  left_join(envdata, by = 'Sample_name')


otu.coord <- otu.umap$layout %>% 
  as_tibble(rownames = 'Sample_name') %>% 
  left_join(envdata, by = 'Sample_name')


# Plot them together ------------------------------------------------------

umap_plot <- function(dataset){
  dataset %>% 
  ggplot(aes(V1,V2)) + 
  geom_point(aes(color = season),
             size = 3, alpha =0.9) + 
  theme_bw() + 
  sea.colScale
  
}

gene.plot <- gene.coord %>% 
  umap_plot() + 
  ggtitle('Functional dataset')
  
otu.plot <- otu.coord %>% 
  umap_plot() + 
  ggtitle('16S dataset')

composite.plot <- gene.plot + 
  otu.plot  + 
  plot_layout(guides = 'collect') 
  

ggsave('results/figures/ordination_variants/16S_vs_genes.pdf',
       width = 9,
       height = 6)
