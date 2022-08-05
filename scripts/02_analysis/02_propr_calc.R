library(tidyverse)
library(propr)

source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')

annot <- read_rds('data/intermediate_files/annotation.rds')
tax <- readRDS('data/04_table_gen/taxonomy_contigwise.rds')

ab.raw.filt <- read_rds('data/04_table_gen/abtblraw_filt0rem_8sams.rds')

ab.scg.filt <- read_rds('data/04_table_gen/abtblraw_allonlyscg.rds')  %>% 
  group_by(genename) %>% 
  filter( n() >= 8) %>% 
  ungroup()


order.sel <- tax %>%
  group_by(order) %>%
  summarize( count = n()) %>%
  filter(count >= 50, !is.na(order)) %>% 
  pull(order)

# we will make a test with SAR86 
genename.sel <- tax %>% 
  filter(family == 'Puniceispirillaceae') %>% 
  pull(genename)


genename.list <- tax %>% 
  select(genename, order) %>% 
  filter(order %in% order.sel) %>% 
  split(.$order) %>% 
  map(~pull(.x, genename))


obtain_proportionality <- function(sel){
  
  dataset <-  ab.raw.filt %>% 
    filter(genename %in% sel) %>% 
    group_by(genename) %>%
    pivot_wider(names_from = genename,
                values_from = count, values_fill = list(count = 0) )  %>% 
    column_to_rownames(var = 'sample')
  
  propr.res <- propr(counts = dataset,
                     metric = 'rho',
                     symmetrize = T,
                     alpha = NA,
                     ivar = "clr")
  
  
  updat.pr <- updateCutoffs(propr.res,
                            cutoff = seq(0, 1, .1))
  
  
  min.cutoff <- updat.pr@fdr %>%
    as_tibble() %>%
    filter(FDR <= 0.05) %>%
    pull(cutoff) %>%
    min()
  
  res.propr <- getResults(updat.pr[">", min.cutoff]) %>% 
    as_tibble()
  
  return(res.propr)
  
}

results.whole <- genename.list %>% 
  map(~obtain_proportionality(.x))

saveRDS(results.whole, 'data/intermediate_files/rho_betweenorders.rds')
