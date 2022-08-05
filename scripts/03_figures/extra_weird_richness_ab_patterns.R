library(tidyverse)
library(mgcv) 
library(patchwork)
source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')

path.fig <- 'results/figures/gene_agg_seasonality/'
dir.create(path = path.fig)

ab.agg <- read_rds('data/04_table_gen/alragg.tbl.rds') 
descriptions <- readxl::read_xlsx('results/gene_info/database_genes.xlsx')
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv') 
seasonality <- read_tsv('results/stats/seasonality_agg.tsv') %>% 
  mutate( seasonal = ifelse(peak >= 8 & pval <= 0.05 & interval.start <= 1,
                            TRUE, FALSE))

sea.genes <- seasonality %>% filter(seasonal) %>% pull(gene_name)

# Merging everything together 
megafile.agg.all <-  ab.agg %>% 
  pivot_longer( names_to = 'sample', values_to = 'logratio', -gene_name) %>% 
  left_join(descriptions, by = 'gene_name') %>% 
  left_join(seasonality, 
            by = c('gene_name')) %>% 
  left_join(envdata, by = c('sample' = 'Sample_name'))  %>% 
  mutate(ratio = 2 ^ logratio)

megafile.agg <- megafile.agg.all %>% 
  filter(gene_name == 'psbA')   

# General pattern ---------------------------------------------------------

megafile.agg %>% 
  ggplot(aes(month, ratio)) +
  geom_point() + 
  geom_smooth() 


psbA_comparison <- megafile.agg %>% 
  ggplot(aes(factor(season, levels = season.order), ratio)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
  xlab('Season') + 
  ggtitle(label = 'psbA seasonality pattern', 
          subtitle = str_c('ANOVA doesnt raise any significant comparison.',
          'At the month level, only one group is significant,',
          'with the significance based in the outlier', sep = '\n') ) + 
  
megafile.agg.all %>% 
  filter(gene_name == 'pufM') %>% 
  ggplot(aes(factor(season, levels = season.order), ratio)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
  xlab('Season') +
  plot_layout(ncol = 1) +
  ggtitle(label = 'pufM seasonality') + 

megafile.agg.all %>% 
  filter(gene_name == 'narB') %>% 
  ggplot(aes(factor(season, levels = season.order), ratio)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
  xlab('Season') +
  plot_layout(ncol = 1) +
  ggtitle(label = 'narB seasonality') 

ggsave(filename = 'results/figures/general/psbA_comparison.pdf',
       psbA_comparison, 
       width = 7, height = 10)

megafile.agg %>% 
  ggplot(aes(factor(month), ratio)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(alpha = 0.8)  +
  scale_x_month

# Is there an statistical difference? 

aov(ratio ~ season, data = megafile.agg) %>% 
  TukeyHSD()



aov(ratio ~ as.factor(month), data = megafile.agg %>% filter(ratio < 0.75)) %>% 
  TukeyHSD() %>% 
  .$`as.factor(month)` %>%
  as_tibble(rownames = 'month') %>% 
  filter(`p adj` <= 0.05)


# and for narB ? 
aov(ratio ~ as.factor(month), data =  megafile.agg.all %>%  filter(gene_name == 'narB', ratio < 0.75)  )  %>%
  TukeyHSD() %>%
  .$`as.factor(month)` %>%
  as_tibble(rownames = 'month') %>%
  filter(`p adj` <= 0.05)

# Comparison richness, ratio and chlA pattern ------------------------------
rich.df <- read_tsv('data/intermediate_files/richness_gene_trends.tsv')

comparison <- megafile.agg %>% 
  select(gene_name, Date,  Chla_3um, Synechococcus) %>% 
  left_join(rich.df %>% select(gene_name, Date, richness),
            by = c('gene_name', 'Date') )  %>% 
  mutate(bloom = ifelse(Synechococcus >= 40000, TRUE, FALSE)) %>% 
  pivot_longer(names_to = 'par',
               values_to = 'val',
               cols = -c(gene_name, Date, bloom)) 


compar.synecho <- comparison %>% 
  ggplot( aes(Date, val)) + 
  geom_line() + 
  geom_point(aes(color =bloom)) +
  facet_wrap(~par, ncol = 1, scale = 'free_y') + 
  ggtitle('The peaks in richness and ratio are coming from Synechococcus')

ggsave(filename = 'results/figures/general/psbA_pattern.pdf',
       plot = compar.synecho,
       width = 9, height = 9)


# The same with dmdA? -----------------------------------------------------


megafile.agg.all %>% 
  filter(gene_name == 'dmdA') %>% 
  select(gene_name, Date, ratio, Chla_3um, Chla_total ) %>% 
  left_join(rich.df %>% select(gene_name, Date, richness),
            by = c('gene_name', 'Date') )  %>% 
  mutate(bloom = ifelse(Chla_total >= 0.75, TRUE, FALSE)) %>%
  filter(Chla_total <= 1) %>% 
  pivot_longer(names_to = 'par',
               values_to = 'val',
               cols = -c(gene_name, Date, bloom))  %>% 
  ggplot( aes(Date, val)) + 
  geom_line() + 
  geom_point(aes(color =bloom)) +
  facet_wrap(~par, ncol = 1, scale = 'free_y') 
