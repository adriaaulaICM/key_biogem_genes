library(tidyverse)
library(patchwork)
library(ggtext)
library(ggforce)
source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')

# Barplot with each variant and each seasonal maxima. 

set.seed(42)

# Data import -------------------------------------------------------------

lomb <- readRDS('data/intermediate_files/genes_seasonal.rds') 

ab.raw.filt <- read_rds('data/04_table_gen/abtblalr_filt0rem_8sams.rds')  %>%
  filter(genename %in% lomb$genename)
relab.total <- readRDS('data/intermediate_files/total_relab_genev.rds') %>% 
  ungroup() %>% 
  select(genename, relab.total)
  
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')
annot <- read_rds('data/intermediate_files/annotation_key_genes.rds')
  
tax <- readRDS('data/04_table_gen/taxonomy_contigwise.rds')  %>% 
  filter(genename %in% lomb$genename) 

sel.season <- read_rds('data/intermediate_files/season_maxima.rds')


pattern.dist.sea <-  annot %>% 
  filter(genename %in% lomb$genename) %>% 
  select(genename, cycle, gene_name) %>% 
  left_join(sel.season, by = 'genename') %>% 
  left_join(relab.total, by = 'genename') %>% 
  left_join(tax, by = 'genename') %>% 
  mutate(family = case_when( family %in% family.selection ~ family,
                             last_rank %in% c("no rank", 'superkingdom')  ~ 'unclassified // no assignation',
                             TRUE ~ 'Other'))  %>% 
  mutate(family = case_when( family == 'Other' &
                               class == 'Alphaproteobacteria' ~ 'Other Alpha',
                             family == 'Other' &
                               class == 'Gammaproteobacteria' ~ 'Other Gamma',
                             TRUE ~ family)) %>%
  mutate(family = factor(family,
                         levels = family.order,
                         labels = family.labels), 
         gene_name = factor(gene_name,
                             levels = gene.order,
                             labels = gene.italics),
         season = factor(season, 
                         levels = season.trans)) %>% 
  filter(!is.na(season)) %>% 
  mutate(season = as.character(season),
         season = case_when(
             season == 'Winter-Spring' ~ 'Winter',
             season == 'Spring-Summer' ~ 'Spring',
             season == 'Summer-Autumn' ~ 'Summer',
             season == 'Autumn-Winter' ~ 'Autumn',
             TRUE ~ season))


pattern.plot <- pattern.dist.sea %>% 
  ggplot(aes(y = fct_rev(factor(season,
                                levels = str_to_title(season.order))),
             x = relab.total)) + 
  geom_bar(stat = 'identity', aes(fill = family)) +
  facet_wrap(~gene_name,
             scales = 'fixed') + 
  lil.strip +  
  ylab('Season') + 
  xlab('Total relative abundance (count variant / Total gene count)') + 
  scale_fill_manual(values = family.colors.euk, name = 'Family') + 
  scale_x_continuous(label = scales::percent) + 
  coord_cartesian(xlim = c(0, 0.55)) + 
  theme(strip.text.x = element_markdown(), panel.grid.major.y = element_blank(),
        legend.text = element_markdown(), 
        legend.position = c(0.8, 0.09),
        legend.direction = 'horizontal') + 
  guides(fill = guide_legend(nrow = 5, title.position = 'top'))


ggsave(filename = 'results/figures/ordination_variants/gene_distribution_seasons.pdf',
       plot = pattern.plot,
       width = 14, 
       height = 10)


taxver.plot <- pattern.dist.sea %>% 
  mutate(cycle = ifelse(cycle == 'phosporous', 'phosphorous', 'other')) %>% 
  ggplot(aes(y = fct_rev(factor(season,
                                levels = str_to_title(season.order))),
             x = relab.total)) + 
  geom_bar(stat = 'identity', aes(fill = gene_name)) +
  facet_wrap(cycle~family,
             scales = 'free') + 
  lil.strip +  
  ylab('Season') + 
  xlab('Total relative abundance (count variant / Total gene count)') + 
  scale_x_continuous(label = scales::percent) + 
  scale_fill_manual(values = palette.it, name = 'Gene name') + 
  theme(strip.text.x = element_markdown(), panel.grid.major.y = element_blank(),
        legend.text = element_markdown(), 
        legend.position = c(0.8, 0.09),
        legend.direction = 'horizontal') 

ggsave(filename = 'results/figures/ordination_variants/gene_tax_distribution_seasons.pdf',
       plot = taxver.plot, 
       width = 15, 
       height = 10)
