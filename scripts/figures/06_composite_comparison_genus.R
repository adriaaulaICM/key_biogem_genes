library(tidyverse)
library(cowplot)
library(patchwork)
library(geofacet)
library(ggtext)
library(ggforce)
source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')



# Data import -------------------------------------------------------------


fig.path <- 'results/figures/compar_mainord'
dir.create(fig.path)

lomb <- readRDS('data/intermediate_files/genes_seasonal.rds') 

ab.alr.filt <- read_rds('data/04_table_gen/abtblalr_filt0rem_8sams.rds')

envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')
annot <- read_rds('data/intermediate_files/annotation_key_genes.rds')
  
tax.extra <- readRDS('data/04_table_gen/taxonomy_contigwise.rds') 
tax.raw <- readRDS('data/04_table_gen/taxonomy_contigwise.rds') 
tax <- readRDS('data/04_table_gen/taxonomy_merged_contig_sgene_gtdb.rds')  %>% 
  filter(genename %in% lomb$genename)
tax.haro <- read_tsv('data/04_taxgenes/taxonomy_genevars_sags_haro.tsv')

contig.len <- read_tsv('data/stats/contig_length.tsv')
sel.season <- read_rds('data/intermediate_files/season_maxima.rds') %>% 
  rename( 'season_max' = season)


sel.genes <- tax.raw %>% 
  filter(family %in%  c('Pelagibacteraceae',
                        'Halieaceae',
                        'HIMB59',
                        'Flavobacteriaceae',
                        'Poseidoniaceae',
                        'Litoricolaceae',
                        'GCA-002718135',
                        'MB11C04',
                        'SAR86',
                        'TMED112',
                        'D2472',
                        'Puniceispirillaceae',
                        'Rhodobacteraceae')) %>%
  pull(genename)


mega <- ab.alr.filt %>% 
  filter(genename %in% sel.genes) %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) %>% 
  left_join(tax, by = 'genename')  %>% 
  left_join(annot %>% select(genename, gene_name),
            by = 'genename') %>% 
  left_join(sel.season, by = 'genename') %>% 
  mutate(season_max = as.character(season_max),
         season_max = case_when(
             season_max == 'Winter-Spring' ~ 'Winter',
             season_max == 'Spring-Summer' ~ 'Spring',
             season_max == 'Summer-Autumn' ~ 'Summer',
             season_max == 'Autumn-Winter' ~ 'Autumn',
             TRUE ~ season_max))
  
palette.cool <- tibble(gene_name = names(palette.it) %>% str_remove_all(pattern = '\\*'),
                       gene_name.it = names(palette.it), 
                       color = palette.it) %>% 
  mutate( label = str_c("<b style='color:", color,"'>", gene_name.it, "</b>"),
          gene_name = factor(gene_name, levels = gene.order)) %>% 
  arrange(gene_name) %>% 
  mutate(label = fct_inorder(label))
  

presence.genes.df <- tax.raw %>% 
  select(genename, family, genus) %>% 
  left_join(annot, by = 'genename') %>% 
  select(genename, family, genus, gene_name) %>% 
  left_join(tax.haro, by = 'genename') %>% 
  mutate(genus = case_when(Genomospecies == 'VII (gMED)' ~ Genomospecies,
                           TRUE ~ genus)) %>% 
  distinct(genus, family, gene_name) %>% 
  left_join(palette.cool, by = 'gene_name')


# Functions ---------------------------------------------------------------

gene_trend_plot <- function(df){
  
  df %>% 
    ggplot( aes(day_of_year, 2^count)) +
    geom_point(aes(color = gene_name),
               size = 0.5,
               alpha = 0.7, 
               show.legend = FALSE) +
    stat_smooth(aes(x = day_of_year,
                    group = genename,
                    color = gene_name),
                method = "gam",
                se = F,
                formula = y ~ s(x, k =12, bs = 'cc'),
                alpha = 0.8,
                show.legend = TRUE)  + 
    lil.strip + 
    scale_dayyear_shrt +
    scale_y_log10() +
    xlab('Month') + 
    ylab('Abundance ratio (log10)')  + 
    ggtitle('Seasonal genes') + 
    theme(panel.background = element_blank(), plot.background = element_blank()) + 
    scale_color_manual(values = palette.all, name = 'Gene') 
  
}

gene_pres_plot <- function(df, name){
  
  df %>%  
    mutate( gene_name = factor(gene_name,
                               levels = c('PR', gene.order),
                               labels = c('PR', ifelse(gene.order %in% gene.italics.sel,
                                                       str_c('*', gene.order, '*'),
                                                       gene.order)))) %>% 
    ggplot(aes(genus, label)) + 
    geom_point(aes(color = gene_name),
               show.legend = F, 
               size = 5.5) + 
    scale_x_discrete(position = "top") +
    scale_color_manual(values = palette.it, name = 'Gene')  + 
    theme_minimal() + 
    theme(axis.text.y = element_markdown(),
          legend.text = element_markdown(),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          plot.title = element_markdown()) + 
    ylab(NULL) + 
    xlab(NULL) + 
    ggtitle('Gene presence') + 
    coord_cartesian(clip = 'off')  
  
}

composite_plot <- function(plot1, plot2, main){
  
  plot2 + 
    plot1 + 
    plot_layout(widths = c(0.2, 0.8)) + 
    plot_annotation(title = main,
                    theme = theme(plot.title = element_markdown(size = 16,
                                                                hjust = 0.05)))
  
  
  }



# Pelagibacteraceae -------------------------------------------------------
pela.sel <- c('AG-414-E02', 'MED-G40', 'Pelagibacter A', 'VII (gMED)',
              'Pelagibacter', 'HIMB114')

pelago.dat <- mega %>% 
  filter(family %in% 'Pelagibacteraceae',
         genus %in% pela.sel,
         !is.na(gene_name)) %>% 
  left_join(tax.haro, by = 'genename') %>% 
  mutate(genus = case_when(Genomospecies == 'VII (gMED)' ~ Genomospecies,
                           TRUE ~ genus))

pelagos <- pelago.dat %>%  
  gene_trend_plot() + 
  guides(color = 'none') + 
  facet_wrap( ~genus, scales = 'fixed') 

genes.pelagos <-  presence.genes.df %>% 
  filter(genus %in% pela.sel,
         !is.na(gene_name)) %>% 
  gene_pres_plot(name = 'Pelagibacteraceae')

compo.pelago <- composite_plot(pelagos,
                               genes.pelagos,
                               main = 'Main *Pelagibacteraceae* genera')
  
ggsave(filename = str_c(fig.path, '/pelagos.pdf'),
       plot = compo.pelago, 
       width = 12, 
       height = 7)


pelag.sp <- mega %>% 
  filter(genus == 'Pelagibacter') %>% 
  filter(species %in% c('Pelagibacter sp003209915',
                        'Pelagibacter sp003212925',
                        'uc_Pelagibacter')) %>% 
  gene_trend_plot() + 
  facet_wrap( ~ species, scales = 'fixed') 


ggsave(filename = str_c(fig.path, '/pelagibacter_sp.pdf'),
       plot = pelag.sp, 
       width = 11, 
       height = 7)

# SAR86 -------------------------------------------------------------------
sar86.sel <- c('D2472')

mega %>% 
  filter(family %in%  sar86.sel,
         !is.na(gene_name)) %>%
  filter(!str_detect(genus, 'uc_')) %>% 
  ggplot( aes(day_of_year, 2^count)) +
  geom_point(aes(fill = gene_name),
             size = 1.5,
             alpha = 0.7,
             # color = NULL,
             stroke = 0.01,
             shape = 21,
             show.legend = TRUE) +
  facet_wrap( ~genus, scales = 'fixed') + 
  stat_smooth(aes(x = day_of_year,
                  group = season_max,
                  color = season_max),
              method = "gam",
              se = F,
              formula = y ~ s(x, k =12, bs = 'cc'),
              alpha = 0.8 ,
              show.legend = TRUE)  +
  lil.strip + 
  scale_dayyear_shrt +
  scale_y_sqrt() +
  xlab('Month') + 
  ylab('Abundance ratio')  + 
  ggtitle('Seasonal genes') + 
  scale_fill_manual(values = palette.all, name = 'Gene') + 
  # guides(fill = 'none') + 
  guides(fill =  guide_legend(override.aes = list(size = 3,
                                                  linetype = 'blank',
                                                  alpha = 1))) +
  theme(panel.background = element_blank(), plot.background = element_blank()) 

sar86 <- mega %>%
  filter(family %in% sar86.sel,
         !is.na(gene_name)) %>%
  gene_trend_plot() +
  guides(color = 'none') +
  facet_wrap( ~genus, scales = 'fixed')

genes.sar86 <-  presence.genes.df %>%
  filter(family %in% sar86.sel,
         !is.na(gene_name)) %>%
  gene_pres_plot(name = 'SAR86')


compo.sar86 <- composite_plot(sar86,
                             genes.sar86,
                             main = 'Main D2472 (SAR86) genera') 


ggsave(filename = str_c(fig.path, '/sar86.pdf'),
       plot = compo.sar86,
       width = 12,
       height = 7)


# Composite ---------------------------------------------------------------


final.composite <- plot_grid(compo.pelago,
                             compo.sar86,
                             ncol = 1,
                             # rel_heights = c(1.2,1),
                             labels = 'AUTO')


ggsave(filename = 'results/figures/compar_mainord/main_compar4.pdf',
       plot = final.composite,
       width = 10, 
       height = 12)

# # Litoricola --------------------------------------------------------------
# 
# # lito.sel <- c('Litoricolaceae') 
# # 
# # lito <- mega %>% 
# #   filter(family %in% lito.sel,
# #          !is.na(gene_name)) %>% 
# #   gene_trend_plot() +   
# #   guides(color = FALSE) + 
# #   scale_y_log10() +
# #   facet_wrap( ~ species, scales = 'fixed') 
# # 
# # genes.lito <-  presence.genes.df %>% 
# #   filter(family %in% lito.sel,
# #          !is.na(gene_name)) %>% 
# #   gene_pres_plot(name = 'lito')
# # 
# # compo.pelago <- lito + 
# #   genes.lito + 
# #   plot_layout(widths = c(0.8, 0.2))
# # 
# # ggsave(filename = str_c(fig.path, '/lito.pdf'),
# #        plot = compo.pelago, 
# #        width = 12, 
# #        height = 7)
# # 
# # 
# 
# 
# # Opitutales --------------------------------------------------------------
# 
# opi <- mega %>%
#   filter(family %in% 'MB11C04', 
#          !is.na(gene_name)) %>%
#   gene_trend_plot() +
#   guides(color = 'none') + 
#   scale_y_log10() +
#   facet_wrap( ~ genus, scales = 'free_y')
# 
# genes.opi <-  presence.genes.df %>%
#   filter(family %in% 'MB11C04',
#          !is.na(gene_name)) %>%
#   gene_pres_plot(name = 'SAR86')
# 
# 
# compo.opi <- composite_plot(opi,
#                               genes.opi,
#                               main = 'Main MB11C04 (Opitutales) genera') 
# 
# 
# ggsave(filename = str_c(fig.path, '/opitutales.pdf'),
#        plot = compo.opi,
#        width = 12,
#        height = 7)
# 
# 
# 
# # HIMB59 ------------------------------------------------------------------
# # HIMB59.sel <- c('HIMB59') 
# # 
# # himb59 <- mega %>% 
# #   filter(order %in% HIMB59.sel) %>% 
# #          # !is.na(gene_name)) %>% 
# #   gene_trend_plot() + 
# #   guides(color = FALSE) + 
# #   scale_y_log10() +
# #   facet_wrap( ~ species, scales = 'fixed') 
# # 
# # genes.himb59 <-  presence.genes.df %>% 
# #   filter(family %in% HIMB59.sel,
# #          !is.na(gene_name)) %>% 
# #   gene_pres_plot(name = 'HIMB59')
# # 
# # compo.himb59 <- himb59 + 
# #   genes.himb59 + 
# #   plot_layout(widths = c(0.8, 0.2))
# # 
# # ggsave(filename = str_c(fig.path, '/sar86.pdf'),
# #        plot = compo.pelago, 
# #        width = 12, 
# #        height = 7)
# 
# # Luminiphilus ------------------------------------------------------------
# 
# # lumi <- mega %>% 
# #   filter(family == 'Halieaceae') %>% 
# #   filter(genus == 'Luminiphilus') %>% 
# #   filter(species %in% c('Luminiphilus sp000227505',
# #                         'Luminiphilus sp002390485',
# #                         'Luminiphilus sp002689915',
# #                         'Luminiphilus sp002456975',
# #                         'uc_Luminiphilus') ) %>% 
# #   gene_trend_plot() + 
# #   facet_wrap( ~ species, scales = 'free_y') 
# # 
# # 
# # genepres.lumi <-  presence.genes.df %>% 
# #   filter(genus == 'Luminiphilus') %>% 
# #   gene_pres_plot(name = 'Luminiphilus')
# # 
# # ggsave(filename = str_c(fig.path, '/lumi_compar.pdf'),
# #        plot = lumi, 
# #        width = 9, 
# #        height = 7)
# 
# 
# # Poseidonaceae ------------------------------------------------------------
# 
# poseido.sel <- mega %>% 
#   filter(family == 'Poseidoniaceae' )
# 
# poseido <- poseido.sel %>% 
#   gene_trend_plot() + 
#   guides(color = 'none') + 
#   facet_wrap( ~genus, scales = 'free_y') 
# 
# 
# genes.poseido <-  presence.genes.df %>% 
#   filter(family == 'Poseidoniaceae') %>% 
#   gene_pres_plot() 
# 
# compo.poseido <- composite_plot(poseido,
#                              genes.poseido,
#                              main = 'Main *Poseidonaceae* genera') 
# 
# 
# ggsave(filename = str_c(fig.path, '/poseido_compar.pdf'),
#        plot = compo.poseido, 
#        width = 12, 
#        height = 7)
# 
# # we want to have some numbers from the contigs 
# 
# sel.l2.l3 <- poseido.sel %>% 
#   filter(genus %in% c('MGIIa-L2', 'MGIIa-L3')) %>% 
#   pull(genename) %>% 
#   unique()
# 
# tax.extra %>% 
#   filter(genename %in% sel.l2.l3) %>% 
#   mutate(contigname = str_extract(genename, 'BL[0-9]*\\.k127_[0-9]*') %>% 
#            str_replace(pattern = '\\.', replacement = '_')) %>%
#   # select(contigname)
#   left_join(contig.len, by = 'contigname') %>% 
#   rename( 'contig_length' = length) %>% 
#   select(genename,
#          contig_length,
#          n_protein_fragments,
#          n_fragments_labelled, n_fragments_coinciding, support_percent,
#          species, genus, family, order, class, phylum) 
#   
# mega  %>% 
#   filter(genename == 'BL101019.k127_1047786_12.pro_716808') %>% 
#   View()
# 
# 
# # Rhodobacteraceae --------------------------------------------------------
# 
# rhodo.sel <- c('Amylibacter', 'LFER01',
#                       'HIMB11', 'Planktomarina', 'MED-G52')
# 
# rhodo <- mega %>% 
#   filter(family == 'Rhodobacteraceae') %>% 
#   filter(genus %in% rhodo.sel) %>% 
#   gene_trend_plot() +
#   guides(color = 'none') + 
#   facet_wrap( ~genus, scales = 'free_y') 
# 
# 
# genes.rhodo <-  presence.genes.df %>% 
#   filter(genus %in% rhodo.sel) %>% 
#   gene_pres_plot(name = 'Rhodobacteraceae') + 
#   theme(plot.title = element_markdown())
# 
# 
# compo.rhodo <- composite_plot(rhodo , 
#                               genes.rhodo ,
#                               main = 'Main *Rhodobacteraceae* genera')
# 
# ggsave(filename = str_c(fig.path, '/rhodo_compar.pdf'),
#        plot = compo.rhodo, 
#        width = 12, 
#        height = 7)
# 

# Proteos Flavos ----------------------------------------------------------
flavo.sel <- c('Flavobacteriaceae')

flavo.dat <- mega %>% 
  filter(family %in% flavo.sel, 
         gene_name == 'PR green (M105)') %>% 
  mutate(season_max = factor(season_max, levels = c('Winter', 'Spring',
                                            'Summer', 'Autumn')))

names(sea.col) <- c('Winter', 'Spring', 'Summer', 'Autumn')

flavos <- flavo.dat %>%  
  ggplot( aes(day_of_year, 2^count)) +
  geom_point(aes(color = season_max),
             size = 0.5,
             alpha = 0.7, 
             show.legend = FALSE) +
  stat_smooth(aes(x = day_of_year,
                  group = genename,
                  color = season_max),
              method = "gam",
              se = F,
              formula = y ~ s(x, k =12, bs = 'cc'),
              alpha = 0.8,
              show.legend = TRUE)  + 
  lil.strip + 
  scale_dayyear_shrt +
  scale_y_log10() +
  scale_color_manual(values = sea.col, name = 'Season maximum') + 
  xlab('Month') + 
  ylab('Abundance ratio (log10)')  + 
  ggtitle('Seasonality of proteorhodopsin M105 gene variants') + 
  theme(panel.background = element_blank(), plot.background = element_blank()) + 
  facet_wrap( ~genus, scales = 'fixed') + 
  leg.bottom

ggsave(filename = 'results/figures/compar_mainord/flavo_proteo.pdf',
       plot = flavos,
       width = 9,
       height = 9)

# # Puniceispirillales------------------------------------------------------------
# 
# puni.sel <- c('HIMB100', 'MED-G116', 'UBA8309')
# 
# puni <- mega %>% 
#   filter(family == 'Puniceispirillaceae',
#          genus %in% puni.sel) %>% 
#   gene_trend_plot() + 
#   guides(color = 'none') + 
#   facet_wrap( ~genus, scales = 'free_y') 
# 
# 
# genes.puni <-  presence.genes.df %>% 
#   filter(family == 'Puniceispirillaceae',
#          genus %in% puni.sel) %>% 
#   gene_pres_plot(name = 'Punisher')
# 
# compo.puni <- composite_plot(puni,
#                              genes.puni,
#                              main = 'Main *Puniceispirillaceae* genera') 
# 
# 
# ggsave(filename = str_c(fig.path, '/puniceispiri_compar.pdf'),
#        plot = compo.puni, 
#        width = 12, 
#        height = 7)