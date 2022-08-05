# this script calculates proportionality for a gene set 

library(tidyverse) 
library(propr) 
source('scripts/analysis/sourcefiles/backbone_params-graphs.R')
source('scripts/analysis/sourcefiles/backbone_functions.R')
source('scripts/analysis/sourcefiles/params_biogem.R')

# Data import  ------------------------------------------------------------

ab.raw <- read_rds('data/04_table_gen/abtblraw_all.rds')
ab.alr.filt <- read_rds('data/04_table_gen/abtblalr_filt0rem_25sams.rds') 
ab.scg <- read_rds('data/04_table_gen/abtblraw_allonlyscg.rds')
alr_denom <- read_tsv('data/04_table_gen/selected_scg_norm.tsv') %>% 
  pull(genename)

annot <- read_rds('data/intermediate_files/annotation.rds')
descriptions <- readxl::read_xlsx('results/summary_charact/database_genes.xlsx')
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')

tax <- read_tsv('data/04_taxgenes/taxonomyResult.tsv',
                col_names = c('genename', 'genomeid', 'last_rank',
                              'rank_name', 'taxonomy')) %>% 
  separate(col = taxonomy,
           into = c("species","genus","family",
                    "order","class","phylum","superkingdom"), sep = ";") %>% 
  mutate( label = case_when(is.na(phylum) ~ rank_name,
                                 !is.na(phylum) ~ str_c(phylum, family,sep = ", "),
                                 TRUE ~ rank_name)) %>% 
  #erase eventually 
  distinct(genename, .keep_all = T)  

seasonal <- readRDS('data/intermediate_files/genes_seasonal.rds') %>% 
  rename( 'genename' = asv) %>% 
  filter(peak >= 10) %>% 
  filter(qval <= 0.01)

top5_genvars <- readRDS('data/intermediate_files/top5vars.rds')

ab.raw.filt <- ab.raw %>% filter(genename %in% seasonal$genename)

totals <- ab.raw.filt %>% 
  group_by(genename) %>% 
  summarize( total = sum(count))
  
ab.mega <- ab.alr.filt %>% 
  left_join(annot, by = 'genename') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) %>% 
  left_join(seasonal %>% select(genename, peak), by = 'genename') %>% 
  mutate( seasonal = ifelse( !is.na(peak), 'seasonal', 'not seasonal'),
          count = ifelse(count < 0, count + 0,count))

# Proportionality calculation ---------------------------------------------

# in order to have all samples for every genename we will switch fast between 
# long and wide format 
ab.raw.filt <- ab.raw.filt %>% 
  pivot_wider(names_from = genename,
              values_from = count, values_fill = list(count = 0) ) %>% 
  pivot_longer(names_to = 'genename', values_to = 'count', cols = -sample)

ab.scg <- ab.scg %>% 
  pivot_wider(names_from = genename,
              values_from = count, values_fill = list(count = 0) ) %>% 
  pivot_longer(names_to = 'genename', values_to = 'count', cols = -sample)

obtain_proportionality <- function(selection){
  # From the whole dataset, this function calculates 
  # the rho proportionality measurement
  dataset <- ab.raw.filt %>% 
    left_join(annot %>% select(genename, gene_name), by = 'genename') %>% 
    filter( gene_name == selection) %>% 
    select(-gene_name) %>% 
    group_by(genename) %>%
    pivot_wider(names_from = genename,
                values_from = count, values_fill = list(count = 0) ) 
  
  
  denom <- ab.scg %>% 
    filter( genename %in% alr_denom) %>% 
    distinct(sample, genename, .keep_all = T) %>% 
    pivot_wider(names_from = genename,
                values_from = count,
                values_fill = list(count = 0))
  
  merged.df <- bind_cols(dataset, select(denom, -sample)) %>% 
    column_to_rownames(var = 'sample')
  
  
  propor <- propr(counts = merged.df, metric = 'rho',
                  symmetrize = T, ivar = alr_denom)
  
  propor@matrix[lower.tri(propor@matrix)] <- NA
  
  propr.long <- as_tibble(propor@matrix) %>% 
    mutate( genename2 = colnames(.)) %>% 
    pivot_longer(names_to = 'genename',
                 values_to = 'rho', cols = -c(genename2)) %>% 
    filter(genename != genename2) %>% 
    filter(!is.na(rho)) %>% 
    filter( !genename2 %in% alr_denom, !genename %in% alr_denom)
  
  return(propr.long)
}

propers <- annot %>% filter(genename %in% seasonal$genename) %>%
  pull(gene_name) %>%
  unique() %>% 
  set_names(map(., ~obtain_proportionality(.x)), .)

besties <- propers %>% 
  map(~filter(.x, rho >= 0.8))

bind_rows(besties, .id = 'gene') %>% 
  group_by(gene) %>% 
  filter(n() > 10) %>% 
  ggplot(aes(rho)) + 
  geom_histogram() + 
  facet_wrap(~gene, scales = 'free_y')


# Network creation and clustering  ----------------------------------------

library(ggraph) 
library(igraph)
library(tidygraph)

create_graph <- function(this){
  graph <- as_tbl_graph(this) %>% 
    to_undirected() %>% 
    left_join(totals, by = c('name' = 'genename')) %>% 
    left_join(tax, by = c('name' = 'genename'))
  
  print(graph)
  
  graph <- graph  %>%
    activate(nodes) %>%
    mutate(components = group_components(),
           clusters = group_edge_betweenness(weights = rho) %>% as.character())
  
  V(graph)$components %>% table() %>% sort() %>% print()
  
  graph <-  graph %>%
    activate(nodes) %>%
    filter(components %in% c(1,2))
  
  V(graph)$components %>% table() %>% sort() %>% print()
  return(graph)
  
}


mod_clusters <- function(df){
  df %>% 
    activate(nodes) %>% 
    mutate(clusters = ifelse( clusters %in% c(1,2,3,4), clusters, NA)) %>% 
    return()
}

thegraphs <-  tibble(genes = names(besties), df = besties) %>% 
  rowwise() %>% 
  dplyr::filter(nrow(df) > 0) %>% 
  ungroup() %>% 
  mutate( tidygraph = map(df, ~create_graph(.x) ))

thegraphs <- mutate(thegraphs, tidymod = map(tidygraph, ~mod_clusters(.x)))


  


 
# Plotting  ---------------------------------------------------------------

plot_sea_graph <- function(graph){
  
  graph <- graph %>% 
    activate(edges) %>% 
    filter(rho >= 0.8) %>% 
    activate(nodes) %>% 
    left_join(top5_genvars %>% select(genename, varnum), by = c('name' = 'genename')) 
  
  ggraph(graph, layout = 'kk') + 
    geom_edge_fan(show.legend = FALSE, color = 'dodgerblue') + 
    geom_node_point(aes(fill = as.character(clusters),
                        shape = as.character(clusters),
                        size = log10(total))) + 
    scale_shape_manual(values = c(`1` = 21, `2` = 22, `3` = 23, `4` = 24 ),
                       na.value = 25, name = 'Clusters') + 
    guides(
      fill = guide_legend(override.aes = list(size = 3,
                                              stroke = 0.5,
                                              shape = 21)),
      shape = guide_legend(override.aes = list(fill = "grey"))) + 
    theme_void() + 
    scale_fill_manual(values = clust_varparcall,
                      na.value = 'black', name = 'clusters') 
  
}

# Palette
clust_varparcall <-  c("#00478E", "#00BFB3", "#DBE442", "#FFB1BB", "#E24585")
names(clust_varparcall) <- 1:5

thegraphs <- thegraphs %>% 
  mutate(plots = map(tidymod, ~plot_sea_graph(.x)))
  
save_mult_plot(thegraphs$plots,
               vec_files = thegraphs$genes,
               common = 'results/figures/network/network_',
               endfile = '.pdf', width = 9, height = 9)


save_toprep_clusters <- function(graphs){
  
  graphs %>% 
    activate(nodes) %>% 
    select(name, total, clusters) %>% 
    as_tibble() %>% 
    group_by(clusters) %>% 
    top_n(n = 1, wt = total) %>% 
    ungroup() %>% 
    return()
  
}

repclusters <- thegraphs %>% 
  mutate(toprep = map(tidymod, ~save_toprep_clusters(.x))) %>% 
  pull(toprep) %>% 
  bind_rows() %>% 
  rename('genename'  = name)

plot_repclust_trends <- function(genenam){
  ab.repclust <- ab.mega %>% 
    filter(genename %in% repclusters$genename) %>% 
    left_join(repclusters %>% select(genename, clusters), by = 'genename') %>% 
    filter(gene_name == genenam) 
  
  taxinfo <-  tax[tax$genename %in% ab.repclust$genename,] %>% 
    as_tibble() %>% 
    select(genename, rank_name, label) %>% 
    left_join(ab.repclust %>% select(genename, clusters), by = 'genename') %>% 
    distinct(genename, .keep_all = T) %>% 
    arrange(clusters) %>% 
    mutate(y = rep(c(6,5), length.out = length(clusters)),
           x = seq(1,11, length.out = length(clusters))  )
  
  ab.repclust %>% 
    ggplot(aes(month, count, color = as.character(clusters)) ) + 
    geom_jitter(shape = 21, show.legend = F) + 
    stat_smooth(aes(x = month,
                    group = as.character(clusters)), # continuous x-axis
                method = "lm",
                formula = y ~ poly(x, 3), se = F, size = 2, show.legend = F) + 
    ggrepel::geom_label_repel(data = taxinfo,
                              aes( x = x, y = y, label = label,
                                   fill = as.character(clusters)),
                              show.legend = F,
                              inherit.aes = F,
                              size = 3,
                              color = 'white',
                              vjust = "inward", hjust = "inward") + 
    ylab('Log ratio (gene / scg)')  +
    scale_color_manual(values = clust_varparcall, na.value = 'black') +
    scale_x_month  +
    scale_fill_manual(values = clust_varparcall, na.value = 'black', name = 'Clusters') + 
    theme_minimal(base_size = 19)
  
}


library(patchwork)

thegraphs <- thegraphs %>% 
  mutate( reptrends = map(genes, ~plot_repclust_trends(.x) + ggtitle(.x))) 

# Pulling the magic together  ---------------------------------------------

join_plots <- function(plot1, plot2){
 plot1 + plot2 + plot_layout(ncol = 2 ) 
  
  }
plots <- map2(thegraphs$reptrends, thegraphs$plots, join_plots)

save_mult_plot(plots = plots,
              vec_files = thegraphs$genes,
              common = "results/figures/main/netwhole",
              endfile = '.pdf',width = 15, height = 9)

