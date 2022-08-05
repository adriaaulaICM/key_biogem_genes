# A plot with the yearly pattern of abundance both by the bacterial count 
# and the SCG

library(tidyverse)
library(ggtext)
library(patchwork)

source('scripts/utils/backbone_params-graphs.r')
source('scripts/utils/backbone_functions.r')
source('scripts/utils/params_biogem.r')


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# data import -------------------------------------------------------------

ab.all <- read_rds('data/04_table_gen/abtblraw_all.rds') %>% 
  group_by(sample) %>% 
  summarize(total = sum(count))

ab.scg <- read_rds('data/04_table_gen/abtblraw_allonlyscg.rds')
annot <- read_rds('data/intermediate_files/annotation_all.rds')
descriptions <- readxl::read_xlsx('results/gene_info/database_genes.xlsx')
envdata <- read_tsv('data/metadata-raw/metadata_blanes_compact_may2017.tsv')


megafile.scg <-  ab.scg %>% 
  left_join(annot %>% select(genename, annotation), by = 'genename') %>% 
  left_join(descriptions, by = 'annotation') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name'))

totals.kegg <- megafile.scg %>% 
  group_by(sample,gene_name) %>% 
  summarize( total = sum(count)) %>% 
  left_join(envdata, by = c('sample'= 'Sample_name'))  %>% 
  filter(total < 150000, !gene_name %in% c('LeuS','ValS'))

geomean.kegg <- totals.kegg %>% 
  group_by(sample) %>% 
  summarise( geomean = gm_mean(total)) %>% 
  left_join(envdata, by = c('sample'= 'Sample_name')) 

# we take away some outliers since we are focusing in the main trend
count.trend <- ggplot(data = totals.kegg, 
                      aes(month, total)) + 
  geom_jitter(alpha = 0.7, size = 2, aes(color = gene_name)) + 
  geom_point(data = geomean.kegg, aes(month, geomean),
             color = 'red') + 
  stat_smooth(aes(x = month),
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = T, show.legend = F) + 
  ylab('Read count (aggregated by SCG gene)') + 
  scale_color_pander(name = 'SCG gene name') + 
  scale_x_continuous(breaks = 1:12,
                     name = 'Month',
                     labels = str_to_title(month.order))


# Plotting the cytometry counts 
cytometry_counts <- envdata %>% 
  filter(year %in% c(2009:2015)) %>% 
  mutate(Bact_Syne = Bacteria_joint + Synechococcus) %>% 
  ggplot( aes(month, Bact_Syne)) + 
  geom_point(alpha = 0.7, size =2) + 
  stat_smooth(aes(group = 1),
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = T, size = 1.2,
              alpha = 0.8) + 
  theme(axis.title.y.left =  element_markdown()) + 
  scale_x_continuous(breaks = 1:12,
                     name = 'Month',
                     labels = str_to_title(month.order))  +
  xlab('Month') + 
  ylab('(cells ml<sup>-1</sup>)')


composite <- count.trend +   cytometry_counts + 
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = 'A')

ggsave(plot = composite,
       filename = 'results/figures/general/scg_cytometry_comparison.pdf',
       width = 8, height = 7)


# Relationships between geomean and other values --------------------------
total.sample <- read_tsv('data/stats/seq_depth_sample.txt') %>% 
  mutate(sample = str_remove(file, pattern = '.*/') %>%
           str_remove(pattern = '_.*')) %>% 
  select(sample, num_seqs) %>% 
  distinct()



# is there a relationship between the geomean and the total read count
corr.geomean.kegg.plot <- geomean.kegg %>% 
  left_join(total.sample, by = 'sample')  %>% 
  ggplot( aes(geomean, num_seqs)) + 
  geom_point(size = 2 ) + 
  geom_smooth(method = 'lm') + 
  ylab('Total read count') + 
  xlab('SCG geometric mean')


composite.2  <- count.trend +   corr.geomean.kegg.plot + 
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = 'A') & 
  trans.back

ggsave(plot = composite.2,
       filename = 'results/figures/general/scg_readcount_compar.pdf',
       width = 8, height = 7)