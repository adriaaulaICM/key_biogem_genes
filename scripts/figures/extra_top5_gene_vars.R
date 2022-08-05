library(tidyverse)
library(ggtext)
library(patchwork)
source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')

fig.path <- 'results/figures/top5_vars'
dir.create(fig.path)

ab.raw.all <- read_rds('data/04_table_gen/abtblraw_all.rds') 

ab.raw.filt <- read_rds('data/04_table_gen/abtblraw_filt0rem_8sams.rds') 
ab.alr <- read_rds('data/04_table_gen/abtblalr_filt0rem_8sams.rds')  
# ab.agg <- read_rds('data/04_table_gen/alragg.tbl.rds') %>% 
#   pivot_longer(names_to = 'sample', values_to = 'logratio', cols = -gene_name) %>% 
#   mutate(logratio = 2 ^ logratio)
  
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')
annot <- read_rds('data/intermediate_files/annotation.rds') %>% 
  select(genename, cycle, gene_name)
  

lomb <- readRDS('data/intermediate_files/genes_seasonal.rds') 
tax <- readRDS('data/04_table_gen/taxonomy_merged_contig_sgene_gtdb.rds')

relab.total <- readRDS('data/intermediate_files/total_relab_genev.rds')

# Megafile ----------------------------------------------------------------

mega.df <- ab.raw.filt %>% 
  left_join(annot, by = 'genename') %>% 
  select(sample, genename, count, gene_name) 

# Variant top order -------------------------------------------------------

n.vars <- 10
 
var.n <- mega.df %>% 
  group_by(gene_name, sample) %>% 
  mutate(relab = count / sum(count))  %>% 
  group_by(gene_name, genename) %>% 
  summarise(m.relab = mean(relab)) %>% 
  group_by(gene_name) %>% 
  arrange(gene_name, -m.relab) %>% 
  top_n(n = n.vars, wt = m.relab) %>% 
  mutate(variant.n = 1:n.vars)  %>% 
  ungroup() %>% 
  select(genename, variant.n) %>% 
  mutate( variant.n = factor(variant.n))
  

trend.df <- ab.alr %>% 
  left_join(annot, by = 'genename') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name'))  %>% 
  left_join(var.n, by = 'genename') %>% 
  left_join(tax %>% select(genename, class, genus, family), by = 'genename') %>% 
  mutate(seasonal = ifelse(genename %in% lomb$genename, 'seasonal', 'no-sea') ) %>% 
  filter(!is.na(variant.n)) 
  
  
tax.vars.df <- var.n %>% 
  left_join(annot %>% select(genename, gene_name), by = 'genename') %>% 
  left_join(tax %>% select(genename, class, order, family, genus, species), 
            by = 'genename')


# Main fig ------------------------------------------------------------------

colors.vars <-  palette_pander(n = n.vars)
names(colors.vars) <- 1:n.vars
scale.var.col <- scale_color_manual(values = colors.vars)
scale.var.fill <- scale_fill_manual(values = colors.vars)

seasonal_plot <- function(gene.sel){
  trend.df %>% 
    filter(gene_name == gene.sel) %>% 
    ggplot( aes(day_of_year, 2^count)) +
    geom_point(aes(color = variant.n, 
                   shape = seasonal),
               size = 2,
               alpha = 0.8, 
               show.legend = FALSE) +
    stat_smooth(aes(x = day_of_year,
                    group = genename,
                    color = variant.n),
                method = "gam",
                se = F,
                formula = y ~ s(x, k =12, bs = 'cc'),
                alpha = 0.8,
                show.legend = FALSE) +
    facet_wrap(~str_c(variant.n, ', ',
                      str_replace_na(class),
                      ', ',
                      str_replace_na(family),
                      ', ',
                      str_replace_na(genus)), ncol = 2) + 
    scale_shape_manual(values = c( seasonal = 1, `no-sea` = 2)) +
    scale_dayyear_shrt + 
    scale.var.col + 
    ylab('Abundance ratio') + 
    xlab('Months (day of the year)') + 
    ggtitle(gene.sel, subtitle = "circles = seasonal ; triangles = non-seasonal")
}

gen.nams <-  trend.df$gene_name %>% unique()

daplots <- gen.nams %>% 
  map(~seasonal_plot(.x))



# Barplot totals ----------------------------------------------------------

rel.df <- relab.total %>% 
  left_join(var.n, by = 'genename') %>% 
  filter(genename %in% trend.df$genename)

bartotal <- function(gene.name.sel = 'amoA'){
  
  rel.df %>% 
    filter(gene_name == gene.name.sel) %>% 
    ggplot(aes(x = 1, y = relab.total)) + 
    geom_bar(stat = 'identity', aes(fill = variant.n), show.legend = F) + 
    scale_y_continuous(limits = c(0,1),labels =  scales::percent,
                       expand = c(0,0)) + 
    scale.var.fill + 
    ylab('Total relative abundance') + 
    xlab(NULL) + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  
}



bar.totals <- gen.nams %>% 
  map(~bartotal(.x))



# Composite ---------------------------------------------------------------

compoplot <- map2(daplots, bar.totals,
                  .f = ~(.x + .y + plot_layout(widths = c(0.9, 0.1))))


filenames <- str_c( fig.path, '/', gen.nams, '.pdf')

map2(filenames, compoplot, ggsave, width = 9 , height = 9 )
# Table -------------------------------------------------------------------

interactive.tab <- DT::datatable(tax.vars.df)


DT::saveWidget(interactive.tab, 'results/tables/top_5_table.html')

