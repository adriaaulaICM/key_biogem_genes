library(tidyverse)
library(ggtext)
library(patchwork)

source('scripts/utils/params_biogem.R')
source('scripts/utils/backbone_params-graphs.R')

ab.alr <- read_rds('data/04_table_gen/abtblraw_all.rds') 
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')
descriptions <- readxl::read_xlsx('results/gene_info/database_genes.xlsx') %>% 
  mutate( cycle = ifelse(cycle == 'phosporous', 'phosphorous', cycle))
#TODO clean this bad name
annot <- read_rds('data/intermediate_files/annotation.rds') %>%
  mutate( cycle = ifelse(cycle == 'phosporous', 'phosphorous', cycle))
  
tax <- read_rds('data/04_table_gen/taxonomy_merged_contig_sgene_gtdb.rds')

# Number of annotated genes  ----------------------------------------------
annot %>% pull(genename) %>% length()

# without hao and nifH and the SCG 
annot %>% 
  filter(!str_detect(cycle, 'scg')) %>% 
  filter(!gene_name %in% c('hao', 'nifH', 'chiA')) %>%
  pull(genename) %>% length()


# Number of variants  -----------------------------------------------------

palette.it <- palette.all[gene.order]
names(palette.it) <- gene.italics
      

total.rich.df <- annot %>% 
  filter(genename %in% unique(ab.alr$genename)) %>% 
  filter(!gene_name %in% c('hao', 'nifH', 'chiA', 'dddP')) %>%
  filter(!str_detect(cycle, 'scg')) %>% 
  group_by(cycle,gene_name) %>% 
  summarize(counts = n()) %>% 
  arrange(-counts) %>% 
  mutate( gene_name = factor(gene_name,
                             levels = gene.order,
                             labels = gene.italics)) %>% 
  mutate(cycle = ifelse(cycle %in% c('sulphur', 'iron'),
                        '(Sulphur + iron)',
                        cycle %>% str_to_title()),
         cycle = factor(cycle, levels = c('Carbon', 'Nitrogen', 'Phosphorous',
                                          '(Sulphur + iron)')))

pr.rich.df <- total.rich.df %>% 
  filter(str_detect(gene_name, 'PR'))

pr.total <- tibble(cycle = factor('Carbon'),
                   gene_name = 'PR',
                   counts = sum(pr.rich.df$counts))

total.gg <- total.rich.df %>%
  filter(!str_detect(gene_name, 'PR')) %>% 
  bind_rows(pr.total) %>% 
  ggplot(aes(y = reorder(gene_name,-counts), x = counts)) + 
  geom_bar(stat = 'identity', aes(fill = gene_name))  + 
  geom_text(aes(label = counts),nudge_x = -0.3, color = 'white') + 
  ggforce::facet_col(vars(cycle), scales = 'free_y', space = 'free') + 
  # facet_wrap(~cycle, ncol = 1, scales = 'free_y') + 
  scale_fill_manual(values = c('PR' = '#3d405b', palette.it),
                    name = 'Gene name') + 
  ylab('Gene name') + 
  xlab(NULL) + 
  guides(fill = 'none') + 
  scale_x_log10(expand = c(0,0)) + 
  lil.strip + 
  theme(axis.text.y = element_markdown(),
        panel.grid.minor = element_blank() )

total.pr <- pr.rich.df %>% 
  ggplot(aes(y = reorder(gene_name,-counts), x = counts)) + 
  geom_bar(stat = 'identity', aes(fill = gene_name))  + 
  geom_text(aes(label = counts),nudge_x = -0.3, color = 'white') + 
  scale_fill_manual(values = palette.it, name = 'Gene name') + 
  ylab(NULL) + 
  xlab('Number of variants (log scale)') +
  guides(fill = 'none') + 
  scale_x_log10(expand = c(0,0), limits = c(1, 14683)) + 
  lil.strip + 
  theme(axis.text.y = element_markdown(),
        panel.grid.minor = element_blank() )
  

# Richness  trend  --------------------------------------------------------

# How the richness behaves. 

rich.trends <- ab.alr %>% 
  left_join(annot, by = 'genename') %>% 
  select(sample, genename, count, gene_name, cycle) %>% 
  group_by(gene_name,sample) %>% 
  summarise(richness = n()) %>% 
  left_join(descriptions, by = 'gene_name') %>%
  left_join(envdata, by = c('sample' = 'Sample_name')) %>% 
  mutate( rich.scaled = base::scale(richness)[,1] , 
          cycle = ifelse(cycle %in% c('sulphur', 'iron'),
                         '(Sulphur + iron)',
                         cycle %>% str_to_title()),
          cycle = factor(cycle, levels = c('Carbon', 'Nitrogen', 'Phosphorous',
                                           '(Sulphur + iron)')))  %>% 
  ungroup()

write_tsv(x = rich.trends,
          file = 'data/intermediate_files/richness_gene_trends.tsv')

rich.trends %>% 
  group_by(gene_name, month) %>% 
  summarize(meantrend = mean(rich.scaled))


trend.gg <- ggplot(rich.trends, aes(month, rich.scaled)) + 
  stat_smooth(aes(x = day_of_year,
                  group = gene_name, 
                  color = gene_name), 
              method = "gam",
              se = FALSE,
              formula = y ~ s(x, k =12, bs = 'cc'), 
              show.legend = FALSE) + 
  facet_wrap(~cycle) + 
  scale_color_manual(values = palette.all, name = 'Gene name') + 
  # guides(color = guide_legend(override.aes = list(size = 2))) + 
  scale_dayyear + 
  lil.strip + 
  theme(panel.grid.minor = element_blank()) + 
  xlab('Month') +
  ylab('Scaled richness')


# Composition ----------------------------------------------------------

compo.gg <- (((total.gg  / total.pr ) + plot_layout(heights = c(0.9, 0.1))) | trend.gg ) + 
  plot_layout(widths = c(0.3,0.7)) + 
  plot_annotation(tag_levels = 'A') & 
  trans.back  

ggsave(filename = 'results/figures/general/richness_composition.pdf',
       plot = compo.gg, 
       width = 13, 
       height = 8)


# Specific genes of interest ----------------------------------------------

# let's look deeply some of these genes as we know they are distributed differently
specific.trends.gg <- rich.trends %>% 
  filter(gene_name %in% c('dmdA', 'psbA', 'narB', 'nasA')) %>% 
                          # "PR green (L105)", "PR green (M105)")) %>% 
  mutate( gene_name = factor(gene_name,
                             levels = gene.order,
                             labels = gene.italics)) %>% 
  ggplot(aes(day_of_year, richness)) + 
  geom_point( aes(fill = as.character(year)), shape = 21, size = 2) +
  stat_smooth(aes(x = day_of_year,
                  group = gene_name), 
              method = "gam",
              se = TRUE,
              formula = y ~ s(x, k =12, bs = 'cc'), 
              show.legend = FALSE) + 
  facet_wrap(~gene_name, scales = 'free_y') + 
  scale_fill_brewer(name = 'Year') + 
  leg.bottom +
  scale_dayyear_shrt + 
  lil.strip + 
  theme(strip.text.x = element_markdown(),
        panel.grid.minor = element_blank()) + 
  xlab('Month') +
  ylab('Number of variants')

ggsave(filename = 'results/figures/general/richness_specifc_trends.pdf',
       plot = specific.trends.gg, 
       width = 9,
       height = 7)



# Taxonomic dominance -----------------------------------------------------

# psbA --------------------------------------------------------------------


dataset <- ab.alr %>% 
  left_join(annot, by = 'genename') %>% 
  filter(gene_name == 'psbA')   %>% 
  select(sample, genename, count, gene_name, cycle) %>% 
  left_join(tax, by = 'genename') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) 


sam.more.40 <- dataset %>% 
  group_by(sample) %>% 
  summarise(count = n()) %>% 
  filter(count > 20) %>% 
  pull(sample)




specifially.syne <- dataset %>%
  filter(family %in% c('Nostocaceae',
                       'Ketobacteraceae',
                       'Gloeomargaritaceae',
                       'Cyanobiaceae')) %>%
  mutate(separ = ifelse(genus %in% 'Synechococcus_C',
                        '*Synechococcus*',
                        'All other groups'),
         sample.more = ifelse(sample %in% sam.more.40,
                              '> 20 variants',
                              '<= 20 variants')) %>%
  ggplot(aes(str_sub(sample, start = 3))) + 
  geom_bar(stat = 'count',
           # position = 'fill',
           aes(fill = family), show.legend = T) + 
  geom_hline(yintercept = 0.5, linetype = 2, size = 0.5) +
  lab.flip + 
  xlab('Samples') + 
  ylab('Number of variants') + 
  lil.strip  + 
  leg.bottom + 
  ggtitle('*psbA* richness with a bimodal distribution') + 
  guides(fill=guide_legend(ncol=2)) + 
  # scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_discrete(name = 'Family') + 
  theme(plot.title = element_markdown(),
        panel.grid.minor = element_blank()) + 
  facet_wrap(~sample.more, scales = 'free_x', ncol = 1)

specifially.syne



compo.gg <- specific.trends.gg +
  specifially.syne +
  plot_annotation(tag_levels = 'A')

ggsave(filename = 'results/figures/general/richness_specifc_trends_compo.pdf',
       plot = compo.gg, 
       width = 13,
       height = 7)


# dmdA --------------------------------------------------------------------

dataset <- ab.alr %>% 
  left_join(annot, by = 'genename') %>% 
  filter(gene_name == 'dmdA') %>% 
  select(sample, genename, count, gene_name, cycle) %>% 
  left_join(tax, by = 'genename') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) 


sam.more.70 <- dataset %>% 
  group_by(sample) %>% 
  summarise(count = n()) %>% 
  filter(count > 90) %>% 
  pull(sample)



# dataset %>% 
#   mutate( sample.more = ifelse(sample %in% sam.more.40,
#                               '> 70 variants',
#                               '<= 70 variants')) %>%
#   group_by(sample.more, sample, family, genus, species) %>% 
#   summarize( richness = n()) %>% 
#   group_by(sample) %>% 
#   mutate( relab.richness = richness / sum(richness)) %>% 
#   filter( family == 'Pelagibacteraceae') %>%
#   # ggplot(aes(sample.more, y = relab.richness)) +
#   # geom_violin() + 
#   # geom_jitter() + 
#   # facet_wrap(~family)
#   
#   ggplot(aes(str_sub(sample, start = 3), y = relab.richness)) +
#   geom_bar(stat = 'identity',
#            aes(fill = species), show.legend = T) +
#   geom_hline(yintercept = 0.5, linetype = 2, size = 0.5) +
#   lab.flip +
#   xlab('Samples') +
#   ylab('Number of variants (%)') +
#   lil.strip  +
#   leg.bottom +
#   guides(fill=guide_legend(ncol=4)) +
#   scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
#   theme(plot.title = element_markdown()) +
#   facet_wrap(~sample.more, scales = 'free_x', ncol = 1)

  

  
specifially.dmda <- dataset %>%
  mutate( sample.more = ifelse(sample %in% sam.more.40,
                              '> 70 variants',
                              '<= 70 variants')) %>%
  ggplot(aes(str_sub(sample, start = 3))) + 
  geom_bar(stat = 'count',
           # position = 'fill',
           aes(fill = order), show.legend = T) + 
  geom_hline(yintercept = 0.5, linetype = 2, size = 0.5) +
  lab.flip + 
  xlab('Samples') + 
  ylab('Number of variants (%)') + 
  lil.strip  + 
  leg.bottom + 
  guides(fill=guide_legend(ncol=4)) + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
  scale_fill_discrete(name = 'Family') + 
  theme(plot.title = element_markdown(),
        panel.grid.minor = element_blank()) + 
  facet_wrap(~sample.more, scales = 'free_x', ncol = 1)

specifially.dmda
