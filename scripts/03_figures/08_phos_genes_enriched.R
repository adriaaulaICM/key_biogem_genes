library(tidyverse)
library(patchwork)
library(ggtext)
source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')

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
  
tax <- readRDS('data/04_table_gen/taxonomy_merged_contig_sgene_gtdb.rds')  %>% 
  filter(genename %in% lomb$genename) %>% 
  #TODO change the taxonomy origin for gene-based tax
  mutate(order = ifelse(order == 'Synechococcales', 'PCC-6307', order))

sel.season <- read_rds('data/intermediate_files/season_maxima.rds')

fam.data16S <- read_rds('~/projects/202001_ts_arizona_chap1/data/analysis/family_agg_relab.rds')
fam.data.total16S <- read_rds('~/projects/202001_ts_arizona_chap1/data/analysis/family_agg_total.rds')

names(palette.phos) <- str_c('*', names(palette.phos), '*')

# Calculating groups with phosphorous favored during summer ---------------

data.relab.by.gene <- tax %>% 
  select(genename, genus, family, order) %>% 
  left_join(annot %>% select(genename, gene_name, cycle), by = 'genename') %>% 
  left_join(sel.season, by = 'genename') %>% 
  left_join(relab.total, by = 'genename') %>% 
  mutate(season = as.character(season),
         season = case_when(
             season == 'Winter-Spring' ~ 'Winter',
             season == 'Spring-Summer' ~ 'Spring',
             season == 'Summer-Autumn' ~ 'Summer',
             season == 'Autumn-Winter' ~ 'Autumn',
             TRUE ~ season))  %>% 
  mutate(season = factor(season, levels = str_to_title(season.order)))

families.enriched <-  data.relab.by.gene %>% 
  filter(cycle == 'phosporous') %>% 
  group_by(family, season) %>% 
  summarise(total.relab.season = sum(relab.total)) %>% 
  filter(total.relab.season == max(total.relab.season),
         str_detect(season, 'Summer')) 

write_rds(families.enriched,
          path = 'data/intermediate_files/families_wphos_enriched.rds')


# Setting up the families to be plotted -----------------------------------
selected.over100 <- fam.data.total16S %>% 
  select(Abundance, family) %>% 
  group_by(family) %>% 
  summarize(total = sum(Abundance)) %>% 
  filter(total >= 100) %>% 
  pull(family)

selected.pgenes.enrich <- families.enriched$family

fams.more10genes <- data.relab.by.gene %>% 
  group_by(family) %>% 
  filter(n() >= 10)

sel.phos <- intersect(selected.over100,
                      selected.pgenes.enrich) %>% 
  intersect(., fams.more10genes$family)

# Regarding the seasonality of the group as a whole and not for the phosphorous genes,
# what is the season with the most abundance?
  
families.other <- data.relab.by.gene %>% 
  filter(cycle != 'phosporous') %>% 
  group_by(family, season) %>% 
  summarise(total.relab.season = sum(relab.total)) %>% 
  filter(total.relab.season == max(total.relab.season)) 


families.other %>% 
  select(family, season) %>% 
  left_join(families.enriched, by = c('family', 'season')) %>% 
  filter(!is.na(total.relab.season))



# Viz ---------------------------------------------------------------------

selection.df <- data.relab.by.gene %>% 
  filter(family %in% sel.phos) %>% 
         # str_detect(family, 'uc_', negate = T)) %>% 
  filter(!is.na(season)) %>% 
  mutate(gene_name = factor(gene_name, levels = gene.order, 
                            labels = ifelse(gene.order %in% gene.italics.sel,
                                            str_c('*', gene.order, '*'),
                                            gene.order))) 


bycycle.patterns <- selection.df %>% 
  ggplot(aes(season, relab.total)) + 
  geom_bar(stat = 'identity', aes(fill = cycle), show.legend = T) +
  facet_wrap(~str_c(order, ', ', family), scales = 'free_y') +  
  lil.strip + 
  lab.flip
  

bygene.patterns <- selection.df %>% 
  filter(cycle == 'phosporous') %>% 
  ggplot(aes(season, relab.total)) + 
  geom_bar(stat = 'identity', aes(fill = gene_name), show.legend = T) +
  facet_wrap(~str_c(order, ', ', family), scales = 'free_y') +  
  scale_fill_manual(values = palette.phos, name = 'Gene name') + 
  scale_y_continuous(labels = scales::percent) + 
  lil.strip + 
  lab.flip +
  leg.bottom + 
  theme(legend.text = element_markdown()) + 
  ylab('Total gene relative abundance') + 
  xlab('Seasonal maxima')


bygene.patterns.other <- selection.df %>% 
  filter(cycle != 'phosporous') %>% 
  ggplot(aes(season, relab.total)) + 
  geom_bar(stat = 'identity', aes(fill = gene_name), show.legend = T) +
  facet_wrap(~str_c(order, ', ', family), scales = 'free_y') +  
  scale_fill_manual(values = palette.phos, name = 'Gene name') + 
  scale_y_continuous(labels = scales::percent) + 
  lil.strip + 
  lab.flip +
  leg.bottom + 
  theme(legend.text = element_markdown()) + 
  ylab('Total gene relative abundance') + 
  xlab('Seasonal maxima')


# Composition with 16S  ---------------------------------------------------

plot.16S <- fam.data16S %>% 
  mutate(family = factor(family, levels = selection.df$family %>% unique())) %>% 
  filter(!is.na(family)) %>% 
  filter(year %in% c(2009:2015)) %>% 
  ggplot(aes(day_of_year,Abundance)) + 
  geom_jitter( alpha = 0.4) + 
  stat_smooth(aes(x = day_of_year,
                  group = OTU),
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = F, size = 2,
              show.legend = F, alpha = 0.7) + 
  facet_wrap(~str_c(order, ', ', family), scales = 'free_y') +  
  scale_dayyear_shrt + 
  scale_y_continuous(labels = scales::percent) + 
  ylab('Relative abundance') + 
  xlab('Month') + 
  lil.strip 


composite.plot <- plot.16S + 
  bygene.patterns + 
  plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = 'A')

ggsave(filename = 'results/figures/phosphorous_fams/phosphorous_compar_16S_all.pdf',
       plot = composite.plot, 
       width = 12,
       height = 12)

# Composition specific groups ---------------------------------------------

famsel <-  c('Pelagibacteraceae', 'D2472')

fam.sel.patt <- fam.data16S %>% 
  mutate(family = factor(family, levels = selection.df$family %>% unique())) %>% 
  filter(!is.na(family)) %>% 
  filter(year %in% c(2009:2015)) %>% 
  filter(family %in% famsel) %>% 
  ggplot(aes(day_of_year,Abundance)) + 
  geom_jitter( alpha = 0.4) + 
  stat_smooth(aes(x = day_of_year,
                  group = OTU),
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = F, size = 2,
              show.legend = F, alpha = 0.7) + 
  facet_wrap(~str_c(order, ', ', family), scales = 'free_y') +  
  scale_dayyear_shrt + 
  scale_y_continuous(labels = scales::percent) + 
  xlab('Month') + 
  ylab('Relative abundance') + 
  lil.strip + 
  theme(panel.grid.minor = element_blank(),
        strip.text.x =element_text(size = 14,
                                   margin = margin(.05, 0, .1, 0, "cm"))
        ) 


gene.sel.patt <- selection.df %>% 
  filter(cycle == 'phosporous',
         family %in% famsel) %>% 
  ggplot(aes(season, relab.total)) + 
  geom_bar(stat = 'identity', aes(fill = gene_name), show.legend = T) +
  facet_wrap(~str_c(order, ', ', family), scales = 'free_y') +  
  scale_fill_manual(values = palette.phos, name = 'Gene name') + 
  scale_y_continuous(labels = scales::percent) + 
  lil.strip  +
  theme(legend.text = element_markdown()) + 
  theme(panel.grid.minor = element_blank(),
        strip.text.x =element_text(size = 14,
                                   margin = margin(.05, 0, .1, 0, "cm"))
        ) + 
  ylab('Total gene relative abundance') + 
  xlab('Seasonal maxima')


compo2 <- fam.sel.patt + 
  gene.sel.patt + 
  plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = 'A') & 
  trans.back  


ggsave(filename = 'results/figures/phosphorous_fams/pho_compar_16S_sel.pdf',
       plot = compo2, 
       width = 9,
       height = 9)


gene.other.df <- selection.df %>% 
  filter(cycle != 'phosporous',
         family %in% famsel)

gene.sel.other <- gene.other.df  %>% 
  ggplot(aes(season, relab.total)) + 
  geom_bar(stat = 'identity', aes(fill = gene_name), show.legend = T) +
  facet_wrap(~str_c(order, ', ', family), scales = 'free_y') +  
  scale_fill_manual(values = palette.it[unique(gene.other.df$gene_name)],
                    name = 'Gene name') +
  scale_y_continuous(labels = scales::percent) + 
  lil.strip  +
  theme(legend.text = element_markdown()) + 
  ylab('Total gene relative abundance') + 
  xlab('Seasonal maxima')