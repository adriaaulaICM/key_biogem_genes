library(tidyverse)
source('scripts/analysis/sourcefiles/backbone_params-graphs.R')
source('scripts/analysis/sourcefiles/backbone_functions.R')
source('scripts/analysis/sourcefiles/params_biogem.R')

# This scripts prints the temporal trend of the genes aggregating all the variants 
# it also plots on top the richness trend. 

theme_set(theme_bw(base_size = 19))

# Data import -------------------------------------------------------------
ab.agg <- read_rds('data/04_table_gen/alragg.tbl.rds') 
descriptions <- readxl::read_xlsx('results/summary_charact/database_genes.xlsx')
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')

outscg <- descriptions %>% filter(str_detect(cycle, 'scg')) %>% pull(annotation)
ab.agg <- ab.agg %>% filter(!annotation %in% outscg)

# Merging everything together 
megafile.agg <-  ab.agg %>% 
  gather(key = 'sample', value = 'logratio', -annotation) %>% 
  left_join(descriptions, by = 'annotation') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name'))

plot_genes <- function(df){
  ggplot(data = df, aes( x = month, y = logratio)) +
    geom_jitter( aes(fill = gene_name),
                 shape = 21, alpha = 0.6, size = 3.2) + 
    stat_smooth(aes(x = month,
                    color = gene_name,
                    group = gene_name), # continuous x-axis
                method = "lm",
                formula = y ~ poly(x, 3)) + 
    ggtitle( unique(df$facets_genes)) + 
    scale_color_manual(values = palette.all) +
    scale_fill_manual(values = palette.all) +
    guides(fill = FALSE, color = FALSE) + 
    lil.strip + 
    scale_x_month + 
    theme(axis.text.y = element_text(size = 17),
          plot.title = element_text(size = 22, face = 'italic'),
          strip.text = element_text(size = 21,face = 'bold',hjust = 0),
          plot.background = element_blank()) + 
    ylab('Log10 ratio (gene / single copy genes)')

  }


megafile.agg %>% 
  filter(annotation == 'CODHI') %>% 
  ggplot(aes( x = month, y = logratio)) +
  geom_jitter( aes(fill = gene_name),
               shape = 21, alpha = 0.6, size = 3.2) + 
  geom_boxplot(aes(fill = gene_name), outlier.colour = 'transparent') + 
  scale_color_manual(values = palette.all) +
  scale_fill_manual(values = palette.all) +
  guides(fill = FALSE, color = FALSE) + 
  lil.strip + 
  scale_x_month + 
  theme(axis.text.y = element_text(size = 17),
        plot.title = element_text(size = 22, face = 'italic'),
        strip.text = element_text(size = 21,face = 'bold',hjust = 0),
        plot.background = element_blank()) + 
  ylab('Log10 ratio (gene / single copy genes)')


trends <- megafile.agg %>%
  split(.$gene_name) %>% 
  map(~plot_genes(.))

write_rds(trends, 'data/intermediate_files/plotA.rds')

filenames <- str_c('results/figures/agg_annotation/',
                   gsub("[^0-9A-Za-z///']","" ,
                        names(trends) ,ignore.case = TRUE),
                   '_distribution_genes.pdf')

# Rich ness  --------------------------------------------------------------

ab.raw <- read_rds('data/04_table_gen/abtblraw_all.rds')
annotation <- read_rds('data/annotation/putative_genes/annotation.rds')

sub <- ab.raw %>% 
  left_join(annotation %>% select(genename, annotation), by = 'genename') %>% 
  left_join(descriptions, by = 'annotation') 
  
rich.raw <- sub %>% 
  group_by(gene_name, sample) %>% 
  summarize( richness = n()) %>% 
  left_join(envdata, by = c("sample" = "Sample_name")) 

richtrends <- rich.raw %>% filter(!is.na(gene_name)) %>%
  pull(gene_name) %>% unique() %>% 
  set_names( map(., ~ggplot( rich.raw %>% filter(gene_name == .x), aes(month, richness)) +
                   geom_jitter( shape = 21, color = 'grey', size =2) +
                   stat_smooth(aes(x = month,
                                   group = gene_name), # continuous x-axis
                               method = "lm",
                               formula = y ~ poly(x, 3)) + 
                   ggtitle(.x) + 
                   theme(axis.line.x = element_blank(),
                         axis.text.x = element_blank(),
                         axis.title.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         panel.background=element_blank(),
                         panel.border=element_blank(),
                         panel.grid.major=element_blank(),
                         panel.grid.minor=element_blank(),
                         plot.background=element_blank()) + 
                   theme(axis.line.y = element_line())) , .)

library(patchwork)

write_rds(richtrends, 'data/intermediate_files/plotB.rds')

main <- seq(1,length(richtrends)) %>% 
  map(~ richtrends[[.x]] +  trends[[.x]] +
        plot_layout(ncol = 1, heights = c(0.2, 0.8)))

names <- rich.raw$gene_name %>% unique() %>% .[-30]  

main <- set_names(main, names)
saveRDS(main, 'data/intermediate_files/plotAB.rds')
map2(filenames, main, ggsave, width = 6 , height = 8 )
