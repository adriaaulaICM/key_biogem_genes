library(tidyverse)
source('scripts/analysis/sourcefiles/backbone_params-graphs.R')
source('scripts/analysis/sourcefiles/backbone_functions.R')
source('scripts/analysis/sourcefiles/params_biogem.R')

# This scripts prints the temporal trend of the genes aggregating all the variants 
#  

theme_set(theme_bw(base_size = 20))

# Data import -------------------------------------------------------------
ab.scg <- read_rds('data/04_table_gen/abtblraw_allonlyscg.rds')
descriptions <- readxl::read_xlsx('results/summary_charact/database_genes.xlsx')
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')


megafile.scg <-  ab.scg %>% 
  left_join(descriptions, by = 'annotation') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name'))


thetops <- megafile.scg %>%
  filter(!annotation %in% c("K01873", "K06942")) %>%
  group_by(annotation,genename) %>% 
  summarize( totals = sum(count)) %>% 
  group_by(annotation) %>% 
  top_n(wt = totals, n = 10) %>% 
  pull(genename)

write_tsv(x = data.frame(genename = thetops),
          path  = 'data/04_table_gen/selected_scg_norm.tsv')

trends <- megafile.scg %>% 
  filter(genename %in% thetops) %>% 
  group_by(genename) %>% 
  mutate(scaled = count)
  # mutate(scaled = scale(count))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geomean <- trends %>% 
  group_by(month,sample) %>% 
  summarise( geomean = gm_mean(count + 1))


opt1 <- ggplot(data = trends %>% filter(scaled < 4000), aes(month, scaled)) + 
  geom_jitter(alpha = 0.7) + 
  stat_smooth(aes(x = month,
                  group = genename), # continuous x-axis
              method = "lm",
              formula = y ~ poly(x, 3),
              se = F, show.legend = F) + 
  geom_point(data = geomean, aes(group =1, y = geomean), color = 'red') +
  stat_smooth(data = geomean, aes(x = month,
                                  y = geomean,
                                  group = 1), # continuous x-axis
              method = "lm",
              formula = y ~ poly(x, 3),
              show.legend = F, color = 'red') + 
  scale_x_month + 
  ggtitle(' Top 10 gene variants SCG KEGG',
          subtitle = 'Red points: geomeans x sample') + 
  ylab('Read count')


geomean.agg <- trends %>% 
  group_by(month, sample, annotation) %>% 
  summarize( total = sum(count)) %>% 
  group_by(month,sample) %>% 
  summarise( geomean = gm_mean(total + 1))

opt2 <- trends  %>%
  group_by(month, sample, annotation) %>% 
  summarize( total = sum(count)) %>% 
  filter(total < 85000) %>% 
  ggplot(aes(month, total)) + 
  geom_jitter() + 
  stat_smooth(aes(x = month,
                  # color = ,
                  group = 1), # continuous x-axis
              method = "lm",
              formula = y ~ poly(x, 3),
              se = F, show.legend = F) + 
  geom_point(data = geomean.agg, aes(group =1, y = geomean), color = 'red') +
  stat_smooth(data = geomean.agg, aes(x = month,
                                  y = geomean,
                                  group = 1), # continuous x-axis
              method = "lm",
              formula = y ~ poly(x, 3),
              show.legend = F, color = 'red') + 
  scale_x_month + 
  ggtitle(' Each SCG kegg agregated',
          subtitle = 'Red points: geomeans x sample') + 
  ylab('Read count') 

opt2

library(patchwork)

opt3 <- ggplot(envdata, aes(month, Bacteria_joint)) + 
  geom_jitter(shape = 21, alpha = 0.9) + 
  stat_smooth(aes(group =1), method = "loess",
              show.legend = F)  + 
  scale_x_month + 
  ggtitle('Flow cytometry counts') + 
  ylab('(cells * ml-1)')
opt3

opt1 + opt2 + opt3 + plot_layout(ncol = 2)

ggsave('results/figures/scg_trends/comparison_geomeans.png',
       width = 15, height = 10)


trends  %>% 
  group_by(month, sample, annotation) %>% 
  summarize( total = sum(count)) %>% 
  ggplot(aes(month, total)) + 
  geom_jitter() + 
  stat_smooth(aes(x = month,
                  # color = ,
                  group = annotation), # continuous x-axis
              method = "lm",
              formula = y ~ poly(x, 3),
              se = F, show.legend = F) + 
  scale_x_month + 
  ylab('Count') + 
  facet_wrap(~annotation, scales = 'free')

ggsave('results/figures/scg_var_annot.png', height = 9, width = 19)
