library(tidyverse)
library(mgcv) 
library(patchwork)
# library(pheatmap)
library(broom)
source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')

# Plotting the annual trends of the 26 selected genes 
# Checking if it is seasonal or not 
path.fig <- 'results/figures/gene_agg_seasonality/'
dir.create(path = path.fig)

ab.agg <- read_rds('data/04_table_gen/alragg.tbl.rds') 
descriptions <- readxl::read_xlsx('results/gene_info/database_genes.xlsx')
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv') 
seasonality <- read_tsv('results/stats/seasonality_agg.tsv') %>% 
  mutate( seasonal = ifelse(peak >= 8 & pval <= 0.05 & interval.start <= 1,
                            TRUE, FALSE))

sea.genes <- seasonality %>% filter(seasonal) %>% pull(gene_name)

outscg <- descriptions %>% filter(str_detect(cycle, 'scg')) %>% pull(gene_name)
ab.agg <- ab.agg %>% filter(!gene_name %in% outscg)

# Merging everything together 
megafile.agg <-  ab.agg %>% 
  pivot_longer( names_to = 'sample', values_to = 'logratio', -gene_name) %>% 
  left_join(descriptions, by = 'gene_name') %>% 
  left_join(seasonality, 
            by = c('gene_name')) %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) 

# In this graph we will present PR in global and the PR variations below 
gen.agg.plot <- megafile.agg %>%
  filter(!str_detect(gene_name, 'PR '),
         !str_detect(gene_name, 'Other'),
         sample != 'BL120313') %>% 
  mutate(gene_name = factor(gene_name, levels = gene.order.main, 
                            labels = ifelse(gene.order.main %in% gene.italics.sel,
                                            str_c('*', gene.order.main, '*'),
                                            gene.order.main))) %>% 
  ggplot( aes(day_of_year, 2^ logratio)) + 
  geom_point(shape = 21, size = 1) + 
  stat_smooth(aes(x = day_of_year,
                  group = gene_name,
                  color = peak), 
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'), 
              show.legend = T) + 
  facet_wrap(~gene_name,
             scale = 'free_y') +
  scale_colour_stepsn(colours = viridis::inferno(n = 5),
                      limits = c(1.5, 24.2),
                      name = 'PN value (recurrence strength)') + 
  scale_dayyear_shrt + 
  lil.strip + 
  theme(strip.text.x = ggtext::element_markdown(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.8, 0.07),
        legend.direction = 'horizontal') + 
  xlab('Day of the year (labelled by month)') + 
  ylab('Ratio (gene read count / geomean(single copy genes))') 

#only PRs
gen.agg.pr <- megafile.agg %>%
  filter(str_detect(gene_name, 'PR ') | str_detect(gene_name, 'Other PR'),
         sample != 'BL120313') %>% 
  mutate(gene_name = factor(gene_name, levels = gene.order, 
                            labels = ifelse(gene.order %in% gene.italics.sel,
                                            str_c('*', gene.order, '*'),
                                            gene.order))) %>% 
  ggplot( aes(day_of_year, 2^ logratio)) + 
  geom_point(shape = 21, size = 1) + 
  stat_smooth(aes(x = day_of_year,
                  group = gene_name,
                  color = peak), 
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'), 
              show.legend = FALSE) + 
  facet_wrap(~gene_name,
             scale = 'free_y', nrow = 1) +
  scale_colour_stepsn(colours = viridis::inferno(n = 5), limits = c(1.5, 24.2)) + 
  scale_dayyear_shrt + 
  lil.strip + 
  theme(strip.text.x = ggtext::element_markdown(),
        panel.grid.minor = element_blank()) + 
  xlab(NULL) + 
  ylab(NULL)


layout <- '
AAAAA
AAAAA
AAAAA
AAAAA
AAAAA
AAAAA
AAAAA
AAAAA
BBBB#
'
composite.plot <- wrap_plots(A = gen.agg.plot, B = gen.agg.pr, design = layout) + 
  plot_annotation(tag_levels = 'A') & 
  trans.back  

# composite.plot <- gen.agg.plot +
#   gen.agg.pr + 
#   ylab(NULL) + 
#   plot_layout(ncol = 1, heights = c(0.9, 0.1)) + 
#   plot_annotation(tag_levels = 'A')

ggsave(str_c(path.fig, 'whole_genes_sea_pr.pdf'), bg = 'transparent',
       plot = composite.plot, 
       width = 10, height = 11)


# Heatmap -----------------------------------------------------------------
# gam.gene.model <- megafile.agg %>%
#   group_by(gene_name) %>%
#   nest()  %>%
#   mutate( model = map(data,
#                       ~gam(logratio ~ s(day_of_year, k =12, bs = 'cc'),
#                            data = .)),
#           predict = map(model, ~predict(.))) %>%
#   unnest(c(data,predict))
# 
# pred.matrix <- gam.gene.model %>%
#   select(gene_name, sample, predict) %>%
#   pivot_wider(names_from = 'gene_name', values_from = 'predict') %>%
#   select(-sample) %>%
#   as.matrix()
# 
# rownames(pred.matrix) <- gam.gene.model$sample %>% unique()
# 
# 
# genes.df <- tibble(genes = colnames(pred.matrix))
# 
# newnames <- genes.df %>% 
#   mutate( genes.italics = lapply(
#     colnames(pred.matrix),
#     function(x) bquote(italic(.(x))))
#   ) %>% 
#   mutate( genes.bold.italics = lapply(
#     colnames(pred.matrix),
#     function(x) bquote(bolditalic(.(x))))
#   )  %>% 
#   mutate( genes.other = lapply(
#     colnames(pred.matrix),
#     function(x) bquote(.(x)))
#   ) %>% 
#   mutate( newnames = case_when(
#     genes %in% sea.genes ~ genes.bold.italics,
#     genes %in% gene.italics.sel[!gene.italics.sel %in% sea.genes] ~ genes.italics,
#     TRUE ~ genes.other
#   ) 
#   ) %>% 
#   pull(newnames)
# 
# newmat <- pred.matrix %>%
#   as_tibble(rownames = 'sample') %>% 
#   left_join(envdata %>%
#               select(Sample_name, month),
#             by = c('sample' = 'Sample_name')) %>%
#   pivot_longer(names_to = 'gene',
#                values_to = 'logratio',
#                cols = c(-sample, -month)) %>%
#   group_by(month,gene) %>%
#   summarize(av.logratio = median(logratio)) %>%
#   ungroup() %>%
#   pivot_wider(names_from = gene,
#               values_from = av.logratio,
#               id_cols = month) %>% select(-month) %>% 
#   as.matrix()
# 
# pheatmap(mat = t(scale(newmat)),
#          color = viridis::plasma(150),
#          cellheight = 18, cellwidth = 18,
#          labels_col = month.order %>% str_to_title(),
#          labels_row = as.expression(newnames),
#          cluster_cols = F,
#          filename = str_c(path.fig,'heatmap_seasonal.pdf'))
# 
# 
# # seasonality test for psbA 
# comm <- megafile.agg %>% 
#   filter(gene_name == 'psbA') %>% 
#   select(logratio) %>% 
#   as.matrix()
# 
# seasonality.test(comm.tab = comm, n = 1000)
