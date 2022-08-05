library(tidyverse)
library(patchwork)

source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')

env.raw <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017_old.tsv')

labels <- readxl::read_xlsx('data/metadata-raw/labels_pretty.xlsx',
                            sheet = 2,
                            col_names = c('Variable', 'Measure'))  %>% 
  mutate(all = str_c(Variable, "~", Measure ))


par.gen <-  c("Day_length", "Temperature", "Salinity", "Secchi",
                "NO2", "NO3", "PO4", "Si", 'Chla_3um',"Chla_total", "BP_FC1.55",
                "Bacteria_joint", "HNA", "LNA", "Synechococcus", "Prochlorococcus_FC",
                "Peuk1", "Peuk2", "PNF_Micro", "PNF2_5um_Micro", "PNF_5um_Micro",
                "Cryptomonas", "Micromonas", "HNF_Micro")

params <- c("Day_length", "Temperature", "Salinity", "Secchi",
    "NO2", "NO3", "PO4", "Si", 'Chla_3um', "Chla_total", 'BP_FC1.55',
    "Bacteria_joint", "Synechococcus", "Prochlorococcus_FC",
    "PNF_Micro", "PNF2_5um_Micro", "PNF_5um_Micro",
    "Cryptomonas", "Micromonas", "HNF_Micro")

env.whole <- env.raw %>% 
  select(year, day_of_year, season, month, one_of(par.gen)) %>% 
  gather(key = 'param', value = 'val',
         -month, -year, -day_of_year, -season,
         factor_key = T) %>%
  # We will only plot from 2009 to 2015
  filter( year %in% c(2009:2015)) %>%
  mutate( parameter = factor(param, levels = par.gen,
                             labels = labels$all)) 

num.days.mnt <- c(0,31,28,31,30,31,30,31,31,30,31,30)
cumnum <- cumsum(num.days.mnt)

env <- env.whole %>% 
  filter(param %in% params[1:8])


bio <- env.whole %>% 
  filter(param %in% c(params[c(9:15)], 'Peuk1', 'Peuk2'))


plotting_trend <- function(df){
  
  df %>% 
    ggplot( aes(day_of_year, val)) + 
    geom_point(alpha = 0.7, aes(color = year))  +
    stat_smooth(aes(group = parameter),
                method = "gam",
                color = 'gray26',
                formula = y ~ s(x, k =12, bs = 'cc'),
                se = T, size = 1.2,
                alpha = 0.8) + 
    facet_wrap(~parameter, ncol = 3,
               labeller = label_parsed,
               scales = 'free_y') +
    ylab(NULL) +
    xlab('Day of the year (month label)') +
    scale_x_continuous(breaks = cumnum, labels = str_to_title(month.order)) +
    lab.flip + 
    lil.strip + 
    leg.bottom + 
    guides(color = guide_colourbar(title = 'Year',
                                   title.vjust = 1, 
                                   barwidth = unit(20, 'lines'),
                                   barheight = unit(.5, 'lines'))) + 
    theme( strip.text = element_text(size = 8)) 
  
}


envplot <- env  %>% plotting_trend()

ggsave('results/figures/metadata_env.pdf',
       plot = envplot, 
       width = 8, 
       height = 9)

options(scipen=-3)


bioplot <- bio %>% plotting_trend()

ggsave('results/figures/metadata_bio.pdf',
       plot = bioplot, 
       width = 8, 
       height = 9)

composite <- envplot + bioplot +
  plot_layout(ncol = 2, guides = 'collect') +
  plot_annotation(tag_levels = 'A') & 
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.position = 'bottom')


ggsave('results/figures/metadata_all.pdf',
       plot = composite, 
       width = 14, 
       height = 8)
