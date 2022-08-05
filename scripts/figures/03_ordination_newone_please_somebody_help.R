library(tidyverse)
library(cowplot)
library(patchwork)
library(ggtext)
library(ggforce)
source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')

set.seed(42)


# Data import -------------------------------------------------------------
fig.path <- 'results/figures/ordination_variants'
dir.create(fig.path)


lomb <- readRDS('data/intermediate_files/genes_seasonal.rds') 

ab.raw <- read_rds('data/04_table_gen/abtblalr_filt0rem_8sams.rds')  
ab.raw.filt <- ab.raw  %>% 
  filter(genename %in% lomb$genename)
  
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')
annot <- read_rds('data/intermediate_files/annotation.rds')
  
tax <- readRDS('data/04_table_gen/taxonomy_merged_contig_sgene_gtdb.rds')  %>% 
  filter(genename %in% ab.raw.filt$genename) 

sel.season <- read_rds('data/intermediate_files/season_maxima.rds')
sel.months <- read_rds('data/intermediate_files/month_maxima.rds')

relab.total <- readRDS('data/intermediate_files/total_relab_genev.rds') %>% 
  ungroup() %>% 
  select(genename, relab.total)

coord.dat <- read_rds('data/intermediate_files/umap_coord.rds') %>% 
  left_join(sel.season, by = 'genename')  %>% 
  left_join(relab.total, by = 'genename') 


# Plot things -------------------------------------------------------------


# Main ordination ---------------------------------------------------------


season.max <- ptol_pal()(4)
names(season.max) <- str_to_title(season.order)

season.max <- c(season.max, 'transition' = 'transparent')

# labels 

label.pos <- data.frame( season = str_to_title(season.order),
                         V1 = c(1.8, 0, 2, -2),
                         V2 = c(-5, 2.5, 8, 0))




distribution <- coord.dat %>% 
  mutate(season = as.character(season)) %>% 
  mutate(season = ifelse(str_detect(season, '-'), 'transition', season)) %>% 
  filter(season != 'transition') %>% 
  ggplot(aes( V1, V2, fill = season)) + 
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level..),
                  show.legend = F) +
  geom_point(data = . %>% filter(season == 'transition'),
             size = 3, shape = 21, alpha = 0.7,
              show.legend = F)  + 
  geom_point(data = . %>% filter(season != 'transition'),
             size = 3, shape = 21, alpha = 0.7,
              show.legend = F)  + 
  geom_label( data = label.pos,
              aes(V1,V2,
                  label = season,
                  fill = season),
              color = 'white',
              inherit.aes = F,
              show.legend = F) +
  theme_void() +
  theme(plot.title = element_text(face = 'bold', size = 14)) + 
  leg.bottom + 
  ggtitle('Annual patterns of seasonal gene variants') + 
  scale_fill_manual(values = season.max)


# Month -------------------------------------------------------------------

month.distr <- coord.dat %>% 
  left_join(sel.months, by = 'genename') %>% 
  mutate(month = factor(month,
                        levels = 1:12,
                        labels = month.order %>% str_to_title())) %>% 
  filter(!is.na(month)) %>% 
  ggplot(aes( V1, V2, fill = month)) + 
  geom_point(size = 3, shape = 21, alpha = 1, show.legend = T)  + 
  theme(plot.title = element_text(face = 'bold', size = 14)) + 
  theme_void() + 
  scale_fill_ptol() + 
  ggtitle('Annual patterns of seasonal gene variants (monthly colored)') +
  leg.bottom 




# Phosphorous -------------------------------------------------------------


coord.dat %>% 
  mutate(season = as.character(season)) %>% 
  mutate(season = ifelse(str_detect(season, '-'), 'transition', season)) %>% 
  filter(cycle == 'phosporous') %>% 
  ggplot(aes( V1, V2, fill = season)) + 
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level..),
                  show.legend = F) +
  geom_point(aes(size = relab.total),
             shape = 21, alpha = 0.7,
              show.legend = F)  + 
  facet_wrap(~gene_name) + 
  scale_fill_manual(values = season.max)

# Table like plot ---------------------------------------------------------

nvar.sea.fam <- lomb %>% 
  left_join(tax %>% select(genename, species, genus, family), by = 'genename') %>% 
  mutate(family = ifelse(is.na(family) | !(family %in% family.order) , 'Other', family)) %>% 
  mutate(family = factor(family,
                         levels = rev(c(family.order, 'Other',  'Total')),
                         labels = rev(c(family.labels, 'Other', 'Total')))) %>% 
  left_join(sel.season, by = 'genename')  %>% 
  mutate(season = as.character(season),
         season = case_when(
             season == 'Winter-Spring' ~ 'Winter',
             season == 'Spring-Summer' ~ 'Spring',
             season == 'Summer-Autumn' ~ 'Summer',
             season == 'Autumn-Winter' ~ 'Autumn',
             TRUE ~ season)) %>% 
  filter(!is.na(season)) %>% 
  group_by(season,family) %>% 
  summarize(count = n()) 

total  <- nvar.sea.fam %>% 
  group_by(season) %>% 
  summarize(count = sum(count)) %>% 
  mutate(family = factor('Total'))

nvar.sea.fam.total <-  bind_rows(nvar.sea.fam, total) %>% 
  group_by(family) %>% 
  mutate(max.col = ifelse(count == max(count), season, 'notmax')) %>% 
  mutate(season = factor(season, levels = str_to_title(season.order))) 


nvar.plot <- nvar.sea.fam.total %>% 
  ggplot(aes(x = season, y = family)) + 
  geom_label(aes(label = count, color = max.col),
             size = 4, label.size = NA) +
  guides(color = FALSE) + 
  scale_color_manual(values = c(season.max, 'notmax' = 'grey5')) + 
  theme_minimal_hgrid() + 
  theme(axis.text.y = element_markdown(size = 12),
        axis.text.x = element_markdown(size = 12)) +
  ggtitle('Count of seasonal variants') + 
  ylab(NULL) + 
  xlab('Season')




# Barplot -----------------------------------------------------------------

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
  xlab('Total relative abundance (read count variant / total count gene)') + 
  scale_fill_manual(values = family.colors.euk, name = 'Family') + 
  scale_x_continuous(label = scales::percent) + 
  coord_cartesian(xlim = c(0, 0.55)) + 
  theme(strip.text.x = element_markdown(), panel.grid.major.y = element_blank(),
        legend.text = element_markdown(), 
        # legend.position = c(0.8, 0.09),
        legend.position = 'bottom',
        legend.direction = 'horizontal') 
  # guides(fill = guide_legend(nrow = 5, title.position = 'top'))



# Composition --------------------------------------------------------------


compo1 <- distribution +  nvar.plot + 
  plot_layout(widths = c(0.70, 0.30)) 

compo2 <- compo1 / pattern.plot + 
  plot_annotation(tag_levels = 'A') &
  trans.back


ggsave(filename = 'results/figures/ordination_variants/ordination_composition_gene.pdf',
       plot = compo2,
       width = 13, 
       height = 14)


