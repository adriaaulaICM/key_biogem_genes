library(tidyverse)
library(cowplot)
library(patchwork)
library(ggtext)
library(ggforce)
source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')

set.seed(42)


# Data import -------------------------------------------------------------
fig.path <- 'results/figures/ordination_variants'
dir.create(fig.path)


lomb <- readRDS('data/intermediate_files/genes_seasonal.rds') 

ab.raw <- read_rds('data/04_table_gen/abtblalr_filt0rem_8sams.rds')  
ab.raw.filt <- ab.raw  %>% 
  filter(genename %in% lomb$genename)
  
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')
annot <- read_rds('data/intermediate_files/annotation.rds')
  
tax <- readRDS('data/04_table_gen/taxonomy_contigwise.rds')  %>% 
  filter(genename %in% ab.raw.filt$genename) 

sel.season <- read_rds('data/intermediate_files/season_maxima.rds')

relab.total <- readRDS('data/intermediate_files/total_relab_genev.rds') %>% 
  ungroup() %>% 
  select(genename, relab.total)

coord.dat <- read_rds('data/intermediate_files/umap_coord_fam.rds') %>% 
  left_join(sel.season, by = 'genename')  %>% 
  left_join(relab.total, by = 'genename') 



