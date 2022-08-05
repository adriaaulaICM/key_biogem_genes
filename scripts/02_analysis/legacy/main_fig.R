library(tidyverse)
library(knitr)
library(gridExtra)
library(patchwork)
source('scripts/analysis/sourcefiles/backbone_functions.R')


mainA <- read_rds('data/intermediate_files/plotA.rds')
mainAB <- read_rds('data/intermediate_files/plotB.rds')
mainB <- read_rds('data/intermediate_files/variant_callingplot.rds')
labels.plot <- read_rds('data/intermediate_files/variant_labelplot.rds')
tableC <- read_rds('data/intermediate_files/finaltable.rds')
# polarD <- read_rds('data/intermediate_files/maxtrendsplots.rds')
barD <- read_rds('data/intermediate_files/bartrendsplots.rds')


mylayout <- "
AAACCCDD
BBBCCCFF
BBBCCCFF
BBBEEEFF
"

split_sentence <- function(char){
  strsplit(char, "(?<=.{24})", perl = TRUE)[[1]] %>% str_c(collapse = ' \n ')
}


tableC <- tableC %>%
  rowwise() %>% 
  mutate( `basic function` = split_sentence(`basic function`)) %>% 
  mutate(season_maximums  = str_replace_all(season_maximums, pattern = '//', '\n')) %>% 
  mutate(class_pred  = str_replace_all(class_pred, pattern = '//', '\n')) %>% 
  mutate(fam_pred  = str_replace_all(fam_pred, pattern = '//', '\n'))

compose_plots <- function(gene){
  
  p1 <- mainAB[[gene]] 
  p2 <- mainA[[gene]] + ylab('Log10 ratio (gene/SCG)') + theme_bw(base_size = 15)
  p3 <- mainB[[gene]]
  p4 <- labels.plot[[gene]]
  p5 <- barD[[gene]] 
  p6 <- tableGrob(tableC %>%
                    filter(gene_name == gene) %>%
                    mutate(relab_seasonal = str_c(relab_seasonal, " %")) %>%
                    t())
  
  
  p1 + p2 + p3 + p4 + p5 + p6 + 
    plot_layout(design = mylayout)
  
}


save_mult_plot(plots = names(mainB) %>% map(~compose_plots(.x)),
              vec_files = names(mainB),
              common = "results/figures/main/whole_",
              endfile = '.pdf',width = 18, height = 9)
  
