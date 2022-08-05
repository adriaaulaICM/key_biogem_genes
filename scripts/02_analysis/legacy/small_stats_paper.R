library(tidyverse)

# This script will calculate various statistics and necessities I need to write down
# the damn paper
# eventually some cool figures will be putted elsewhere

source('scripts/analysis/sourcefiles/data_import.R')
source('scripts/analysis/sourcefiles/params_biogem.R')
source('scripts/analysis/sourcefiles/backbone_params-graphs.R')

  
table <- read_rds('data/intermediate_files/finaltable.rds')

# Number of annotated genes  ----------------------------------------------
annot <- annot %>% 
  distinct(genename, .keep_all = T) 
  
annot %>% pull(genename) %>% length()

# without hao and nifH and the SCG 

annot %>% 
  filter(!str_detect(cycle, 'scg')) %>% 
  print() %>% 
  filter(!gene_name %in% c('hao', 'nifH')) %>%
  pull(genename) %>% length()


# Num genes present with read counts  -------------------------------------
mega.ab <- ab.raw %>% 
  left_join(annot, by = 'genename')  %>% 
  filter(!gene_name %in% c('hao', 'nifH', 'chiA')) %>%
  filter(!str_detect(cycle, 'scg')) 


mega.ab %>% 
  pull(genename) %>% unique %>% length()
  

# Richness  trend  --------------------------------------------------------

# How the richness behaves. 

richscaled <- mega.ab %>% 
  group_by(gene_name,sample) %>% 
  summarise(richness = n()) %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) %>% 
  mutate( rich.scaled = scale(richness)) 

richscaled %>% 
  group_by(gene_name, month) %>% 
  summarize(meantrend = mean(rich.scaled))


ggplot( richscaled, aes(month, rich.scaled, fill = gene_name)) + 
  geom_boxplot(show.legend = F) + 
  facet_wrap(~gene_name) + 
  scale_fill_manual(values = palette.all) + 
  lil.strip + 
  scale_x_month + 
  ylab('Scaled richness')

ggsave('results/figures/richness_trend.pdf', width = 15, height = 8)

# We see that narB, nasA and UreC present a contrary tendency. 
# is this due a lack of variants? LeTs checK it OUt 
table  %>% 
  filter(gene_name %in% c('narB', 'nasA', 'ureC', 'psbA')) 

# they are good enough


# Aggregated seasonality  -------------------------------------------------

# Check and compare the seasonality differences between the main groups behaving 
seaagg.genes <- table  %>% filter(agg_seasonality == 'yes') %>% pull(gene_name)

megafile.agg.sea <-  ab.agg %>% 
  gather(key = 'sample', value = 'logratio', -annotation) %>% 
  left_join(descriptions, by = 'annotation') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) %>% 
  filter(gene_name %in% seaagg.genes)  %>% 
  filter(!str_detect(cycle, 'scg')) %>% 
  mutate(logratio = ifelse(logratio < 0, 0, logratio))
 
ggplot(data = megafile.agg.sea, aes( x = month, y = logratio)) +
  geom_jitter( aes(fill = gene_name),
               shape = 21, alpha = 0.6, size = 3.2) + 
  stat_smooth(aes(x = month,
                  color = gene_name,
                  group = gene_name), # continuous x-axis
              method = "lm",
              formula = y ~ poly(x, 3)) + 
  facet_wrap(~cycle, scales = 'free_y') + 
  scale_color_manual(values = palette.all) +
  scale_fill_manual(values = palette.all) +
  # guides(fill = FALSE, color = FALSE) + 
  lil.strip + 
  scale_x_month + 
  theme(axis.text.y = element_text(size = 17),
        plot.title = element_text(size = 22, face = 'italic'),
        strip.text = element_text(size = 21,face = 'bold',hjust = 0),
        plot.background = element_blank()) + 
  ylab('Log10 ratio (gene / single copy genes)')
