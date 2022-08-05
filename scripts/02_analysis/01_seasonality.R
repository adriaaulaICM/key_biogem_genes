library(tidyverse)
library(lomb)
library(lubridate)
source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')

#For each gene variant and for each gene aggregating all the variants 
# we calculate the seasonality based in the lomb scargle stat

# Data import -------------------------------------------------------------
annot <- read_rds('data/intermediate_files/annotation.rds')
ab.alr.filt <- read_rds('data/04_table_gen/abtblalr_filt0rem_8sams.rds')  %>% 
  left_join(annot, by = 'genename') %>% 
  filter(gene_name != 'dddP')

envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')

ab.df0 <- ab.alr.filt  %>%
  correct_samnames() %>%
  left_join(envdata, by = c('sample' = 'Sample_name')) %>%
  mutate(decimaldat = decimal_date(Date))


lomb0 <- ab.df0 %>%
  split(.$genename) %>%
  map(~randlsp( x =.x$count,
                times = .x$decimaldat,
                from = 0.5, to = 3,
                type = 'period',
                plot = F))

resgene0 <- tibble( genename = names(lomb0),
                   pval = map_dbl(lomb0, ~.x$p.value),
                   peak = map_dbl(lomb0, ~.x$peak),
                   interval.range  = map(lomb0, ~.x$peak.at) )  %>%
  mutate( interval = map_dbl(interval.range, ~mean(.x))) %>%
  mutate(qval = fdrtool::fdrtool(pval, statistic="pvalue")$qval)



write_rds(resgene0, 'data/intermediate_files/genes_seasonality_test.rds')

# Saving only the seasonal ones 
resgene0 %>% 
  filter(peak >= 8, qval <= 0.05, interval <= 1.00106515 ) %>% 
  saveRDS(object = ., file = 'data/intermediate_files/genes_seasonal.rds')

# Aggregated gene seasonality ---------------------------------------------
# The clr is already applied in 00_dataimport
ab.agg <- read_rds('data/04_table_gen/alragg.tbl.rds') 

agg.res <- ab.agg %>% 
  pivot_longer(names_to = 'sample',
               values_to = 'logratio',
               cols = -gene_name) %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) %>% 
  dplyr::arrange(decimal_date) %>% 
  split(.$gene_name) %>% 
  map(~randlsp( x =.x$logratio,
                times = .x$decimal_date,
                from = 0.5, to = 3,
                type = 'period',
                plot = F))

resdf <- tibble( gene_name = names(agg.res ),
                 pval = map_dbl(agg.res , ~.x$p.value),
                 peak = map_dbl(agg.res , ~.x$peak),
                 interval.range  = map(agg.res , ~.x$peak.at),
                 interval.start = map_dbl(agg.res, ~.x$peak.at[1]),
                 interval.end = map_dbl(agg.res, ~.x$peak.at[2]))  %>% 
  mutate( interval = map_dbl(interval.range, ~mean(.x)))
  
plot_periodogram <- function(sel, lombobj){
  
  map(lombobj[names(lombobj) %in% sel], ~tibble( scanned = .x$scanned,
                             power = .x$power)) %>%
    bind_rows(.id = 'sel') %>%
    ggplot(aes(scanned, power)) +
    geom_point() + 
    geom_hline(yintercept = 10, linetype =2 , color = 'grey2') + 
    geom_hline(yintercept = 8, linetype =2 , color = 'dodgerblue') + 
    geom_vline(xintercept = 1, linetype =2 , color = 'grey2') + 
    geom_line(aes(group = sel)) +
    facet_wrap(~sel, scales = 'free_y') + 
    scale_x_continuous(breaks = c(0,1,2,3,4,5), limits = c(0,5)) 
  
}

periodoplot <- plot_periodogram(gene.order, agg.res) + 
  lil.strip

ggsave(filename = 'results/figures/gene_agg_seasonality/periodogram_genes.pdf',
       plot = periodoplot, 
       width = 10, height = 8)

resdf %>% 
  select(gene_name, peak,  interval, interval.start, interval.end, pval) %>% 
  arrange(pval) %>% 
  write_tsv('results/stats/seasonality_agg.tsv')
