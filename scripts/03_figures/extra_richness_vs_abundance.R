library(tidyverse)
library(ggtext)
library(patchwork)

source('scripts/utils/params_biogem.R')
source('scripts/utils/backbone_params-graphs.R')

ab.alr <- read_rds('data/04_table_gen/abtblraw_all.rds') 
ab.agg <- read_rds('data/04_table_gen/alragg.tbl.rds') %>% 
  pivot_longer(names_to = 'sample', values_to = 'logratio', cols = -gene_name)
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')
descriptions <- readxl::read_xlsx('results/gene_info/database_genes.xlsx')
annot <- read_rds('data/intermediate_files/annotation.rds')
  
tax <- read_rds('data/04_table_gen/taxonomy_formatted.rds')


richness <- ab.alr %>% 
  left_join(annot, by = 'genename') %>% 
  group_by(gene_name, sample) %>% 
  summarize(richness = n()) %>% 
  ungroup() %>% 
  print()


comparison.df <- richness %>% 
  left_join(ab.agg, by = c('gene_name', 'sample')) %>% 
  left_join(envdata, by = c('sample' = 'Sample_name'))  %>% 
  mutate(gene_name = factor(gene_name, levels = gene.order, 
                            labels = ifelse(gene.order %in% gene.italics.sel,
                                            str_c('*', gene.order, '*'),
                                            gene.order))) 

stat.df <- comparison.df %>% 
  group_by(gene_name) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(richness ~ logratio, data = .x))) %>% 
  mutate( tidy.df = map(model, ~broom::tidy(.x)) ) %>% 
  unnest(tidy.df) %>% 
  filter(term == 'logratio') %>% 
  mutate(estimate = round(estimate, digits = 1)) %>% 
  select(gene_name, estimate,std.error, p.value)

rich.vs.ab.plot <- comparison.df  %>% 
  left_join(stat.df, by = 'gene_name') %>% 
  ggplot(aes(richness, logratio)) + 
  geom_point(aes(fill = gene_name), shape = 21, show.legend = F) + 
  geom_smooth(method = 'lm', 
              aes(group = gene_name,
                  color = p.value <= 0.01)) + 
  facet_wrap(~str_c(gene_name, ', beta: ', estimate), scales = 'free')  +
  scale_fill_manual(values = palette.it) +  
  lil.strip + 
  ylab('Log2 ratio (gene read count / geomean(single copy genes))') + 
  xlab('Richness') + 
  theme(strip.text.x = ggtext::element_markdown()) + 
  theme(legend.position = c(0.9, 0.1))

ggsave(filename = 'results/figures/general/rich_vs_ratio.pdf',
       plot = rich.vs.ab.plot, 
       width = 11,
       height = 9)
