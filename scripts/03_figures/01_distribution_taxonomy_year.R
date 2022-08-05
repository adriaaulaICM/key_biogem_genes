library(tidyverse)
library(ggtext)
library(patchwork)
source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')


fig.path <- 'results/figures/taxonomy_distr_year/'
dir.create(fig.path)

ab.raw.filt <- read_rds('data/04_table_gen/abtblraw_filt0rem_8sams.rds') 
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')
annot <- read_rds('data/intermediate_files/annotation.rds')
  

tax <- readRDS('data/04_table_gen/taxonomy_merged_contig_sgene_gtdb.rds')

# tax <- readRDS('data/04_table_gen/taxonomy_subset.rds')
tax.uniref <- read_tsv('data/04_taxgenes/unirefResult.tsv',
                col_names = c('genename', 'genomeid', 'last_rank',
                              'rank_name', 'taxonomy')) %>% 
  filter(genename %in% ab.raw.filt$genename) %>% 
  separate(col = taxonomy,
           into = c("species","genus","family",
                    "order","class","phylum","superkingdom"), sep = ";") %>% 
  mutate( label = case_when(is.na(phylum) ~ rank_name,
                            !is.na(phylum) ~ str_c(phylum, family,sep = ", "),
                            TRUE ~ rank_name)) %>% 
  distinct(genename, .keep_all = T)

ab.mega <- ab.raw.filt %>% 
  left_join(annot, by = 'genename') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) 

pr.data <- ab.mega %>% 
  filter(str_detect(gene_name, "PR")) %>% 
  mutate(gene_name = 'PR')

ab.mega <- bind_rows(ab.mega, pr.data) 
# Formatting the data -----------------------------------------------------
gene.data.persample <- ab.mega %>%
  select(sample, genename, gene_name, count, month) %>% 
  left_join(tax, by = 'genename') %>% 
  mutate(family = case_when( family %in% family.selection ~ family,
                             last_rank %in% c("no rank", 'superkingdom')  ~ 'unclassified // no assignation',
                             TRUE ~ 'Other'))  %>% 
  mutate(family = case_when( family == 'Other' &
                               class == 'Alphaproteobacteria' ~ 'Other Alpha',
                             family == 'Other' &
                               class == 'Gammaproteobacteria' ~ 'Other Gamma',
                             TRUE ~ family)) %>%
  # mutate(family = factor(family,
  #                        levels = family.order,
  #                        labels = family.labels)) %>% 
  group_by(gene_name, sample) %>% 
  mutate(relab = count / sum(count)) 

# we need to change the Other Gammaproteobacteria to Eukaryotes 

theeukas.photo <- gene.data.persample %>%
  filter(gene_name %in% c('psbA', 'rbcL I'), 
         family %in% c('Other Gamma', 'Other')) %>%
  select(genename) %>%
  left_join(tax.uniref, by = 'genename') %>%
  filter(superkingdom == 'Eukaryota') %>%
  pull(genename)

gene.data.persample <- gene.data.persample  %>% 
  mutate(family = ifelse( genename %in% theeukas.photo, 'Eukaryota', family)) %>% 
  mutate(family = factor(family,
                          levels = family.order.euk,
                          labels = family.labels.euk)) 
  

gene.data <- gene.data.persample %>% 
  group_by(gene_name, family, month) %>% 
  summarize(relab.month = sum(relab))

# We complete the data adding the 0 values to the samples without presence of a group
gene.data.c <- complete(gene.data %>% ungroup(),
                        month, nesting(gene_name, family),
                        fill = list(relab.month = 0)) %>% 
  mutate( gene_name = factor(gene_name,
                             levels = c('PR', gene.order),
                             labels = c('PR', ifelse(gene.order %in% gene.italics.sel,
                                             str_c('*', gene.order, '*'),
                                             gene.order))))


# we want to have an idea of how the relative abundance is compared between 
# the other groups and the selected ones. 
gene.data %>%
  mutate( relab.month = relab.month / 7) %>%
  mutate( type = ifelse(str_detect(family, 'Other'), 'other', 'groups')) %>% 
  group_by(type,gene_name) %>% 
  summarise(m.relab = median(relab.month))  %>% 
  # group_by(type) %>% 
  # summarise(m.relab = median(m.relab))  %>% 
  print()
  



# Plotting ----------------------------------------------------------------
  
areaplot.all <- gene.data.c %>% 
  filter(!gene_name %in% c("PR blue", "PR green (M105)",
                          "PR green (L105)", "Other PR")) %>% 
  ggplot( aes(month, relab.month, fill = family)) + 
  geom_area( stat = 'identity', position = 'fill') + 
  facet_wrap(~gene_name) + 
  lil.strip +
  scale_fill_manual(values = family.colors.euk, name = 'Family') + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
  scale_x_continuous(breaks = 1:12,
                     labels = month.order %>%
                       str_to_title() %>%
                       str_sub(1,1),
                     expand = c(0,0),
                     name = 'Month') + 
  xlab(NULL) + 
  ylab('Relative abundance (samples aggregated by month)') +
  theme(strip.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.position = c(0.7, 0.08)) + 
  guides(fill = guide_legend(nrow = 5, title.position = 'left'))


areaplot.pr <- gene.data.c %>% 
  filter(gene_name %in% c("PR blue", "PR green (M105)",
                          "PR green (L105)", "Other PR")) %>% 
  ggplot( aes(month, relab.month, fill = family)) + 
  geom_area( stat = 'identity', position = 'fill', show.legend = FALSE) + 
  facet_wrap(~gene_name, nrow = 1) + 
  lil.strip +
  scale_fill_manual(values = family.colors.euk, name = 'Family') + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
  scale_x_continuous(breaks = 1:12,
                     labels = month.order %>%
                       str_to_title() %>%
                       str_sub(1,1),
                     expand = c(0,0),
                     name = 'Month') + 
  ylab(NULL) +
  theme(strip.text.x = element_markdown(),
        legend.text = element_markdown())
  
layout <- '
AAAAA
AAAAA
AAAAA
AAAAA
AAAAA
AAAAA
AAAAA
BBBB#
'
composite.plot <- wrap_plots(A = areaplot.all, B = areaplot.pr, design = layout) + 
  plot_annotation(tag_levels = 'A') & 
  trans.back  

ggsave(plot = composite.plot,
       filename = str_c(fig.path, 'varplot_all.pdf'),
       width = 10,
       height = 11)


# Plotting separately each sample -----------------------------------------
megaplot_bars <- gene.data.persample %>% 
  group_by(gene_name, family, month, sample) %>% 
  summarize(relab = sum(relab)) %>% 
  mutate( gene_name = factor(gene_name,
                             levels = gene.order,
                             labels = ifelse(gene.order %in% gene.italics.sel,
                                             str_c('*', gene.order, '*'),
                                             gene.order))) %>% 
  mutate(month = factor(month,
                        levels = 1:12,
                        labels = month.order %>% str_to_title())) %>% 
  ggplot( aes(reorder(sample, month), relab)) + 
  geom_bar(aes(fill = family),
           stat = 'identity',
           position = 'fill') + 
  scale_fill_manual(values = family.colors.euk, name = 'Family') + 
  facet_grid(gene_name~month, scales = 'free_x') + 
  lil.strip +
  lab.flip + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
  ylab('Relative abundance (samples aggregated by month)') +
  theme(strip.text.y = element_markdown(),
        legend.text = element_markdown())

ggsave(plot = megaplot_bars, 
       filename = str_c(fig.path, 'varplot_bysample.pdf'),
       width = 13,
       height = 25)
