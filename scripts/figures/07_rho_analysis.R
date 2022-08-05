library(tidyverse)
library(patchwork)
library(ggtext)
library(ggforce)
source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')


rho <- readRDS('data/intermediate_files/rho_betweenorders.rds')


rho.df <- rho %>% 
  bind_rows(.id = 'order') %>% 
  filter(propr >= .7)

