library(tidyverse)

source('scripts/analysis/sourcefiles/backbone_params-graphs.R')
source('scripts/analysis/sourcefiles/params_biogem.R')
theme_set(theme_bw(base_size = 20))

# Data import -------------------------------------------------------------
ab.agg.tbl <- read_rds('data/04_table_gen/scglengthnorm_aggkegg.tbl.rds') %>% 
  as_tibble() %>% 
  filter(!keggid %in% scg)  

colnames(ab.agg.tbl) <-  str_replace_all(colnames(ab.agg.tbl),
                                     c( "BL100413" = "BL100414",
                                        "BL110704" = "BL110705",
                                        "BL120518" = "BL120511",
                                        "BL131204" = "BL131215"))


# Gene information + envirionmental data 
descriptions <- readxl::read_xlsx('results/summary_charact/database_genes.xlsx')
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')

# Join all together 
megafile.agg <-  ab.agg.tbl %>% 
  gather(key = 'sample', value = 'logratio', -keggid) %>% 
  left_join(descriptions, by = 'keggid') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name'))


plot_genes <- function( df){
  ggplot(data = df, aes( x = month, y = logratio)) +
    geom_jitter( aes(fill = gene_name),
                 shape = 21, alpha = 0.6, size = 3.2) + 
    stat_smooth(aes(x = month,
                    color = gene_name,
                    group = gene_name), # continuous x-axis
                method = "lm",
                formula = y ~ poly(x, 3)) + 
    # facet_wrap(~gene_name, scales = 'free_y' ) +
    # ylab('Log10 ratio (gene / single copy genes)') + 
    ylab(NULL) + 
    ggtitle( unique(df$gene_name)) + 
    scale_color_manual(values = palette.all) +
    scale_fill_manual(values = palette.all) +
    guides(fill = FALSE, color = FALSE) + 
    lil.strip + 
    scale_x_month + 
    theme(axis.text.y = element_text(size =19),
          plot.title = element_text(size = 22, face = 'italic'),
          strip.text = element_text(size = 21,face = 'bold',hjust = 0),
          plot.background = element_blank())
  
  }

trends <- megafile.agg %>%
  filter(gene_name %in% normal.genes) %>% 
  filter(!gene_name %in% c('chiA', 'nifH')) %>% 
  split(.$gene_name) %>% 
  map(~plot_genes(.) )

trends[c(1,5,9,11)] <-  trends[c(1,5,9,11)]  %>% 
  map(~.x + 
    ylab('Log10 ratio (gene / single copy genes)'))

trends[[2]]

filenames <- str_c('results/figures/',
                   str_sub(names(trends),1,4),
                   '_distribution_genes.pdf')

map2(filenames, trends, ggsave, width =6 , height =5 )

trends.double <- megafile.agg %>% 
  filter(gene_name %in% double.genes) %>% 
  filter(!gene_name %in% c('chiA', 'nifH')) %>% 
  split(.$cycle) %>% 
  map(~plot_genes(.)) 
        
filenames.dou <- str_c('results/figures/',
                   str_sub(names(trends.double),1,4),
                   '_distribution_genes.pdf')

map2(filenames.dou, trends.double, ggsave, width =6 , height =5 )
        

# together <- ab.raw %>% 
#   gather(key = 'sample', value = 'count', -keggid, -genename) %>% 
#   group_by(sample, keggid) %>% 
#   summarize( total  = sum(count)) %>% 
#   spread(key = 'sample', value = 'total') %>% 
#   select(-keggid) %>% 
#   as.data.frame()
  
  

# Plot for SCG  -----------------------------------------------------------

scgdf <- readRDS('data/04_table_gen/scglengthnorm.tbl.rds')

colnames(scgdf) <-  str_replace_all(colnames(scgdf),
                                         c( "BL100413" = "BL100414",
                                            "BL110704" = "BL110705",
                                            "BL120518" = "BL120511",
                                            "BL131204" = "BL131215"))


scgwhole <- scgdf %>% 
  gather(key = 'sample', value = 'geomean', -genename, -keggid) %>% 
  left_join(envdata, by = c('sample' = 'Sample_name'))

sctotals <- readRDS('data/04_table_gen/sctotals.tbl.rds') %>% 
  filter(normcount > 0)

sctotals$sample <- str_replace_all(sctotals$sample, 
                                   c( "BL100413" = "BL100414",
                                      "BL110704" = "BL110705",
                                      "BL120518" = "BL120511",
                                      "BL131204" = "BL131215"))


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

scdf <- sctotals %>% 
  filter(!keggid %in% c("K01873", "K06942")) %>%
  group_by(sample, keggid) %>% 
  summarize( totals = sum(normcount)) %>% 
  group_by(sample) %>% 
  summarise( geomean = gm_mean(totals))  %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) 

  
ggplot(data = scdf, aes( x = month, y = geomean)) +
  geom_jitter( shape = 21, alpha = 0.6, size = 3.2) + 
  stat_smooth(aes(x = month,
                  group = 1), # continuous x-axis
              method = "lm",
              formula = y ~ poly(x, 3))   + 
  scale_x_month + 
  ylab('Geometric mean (read counts)')  + 
  ggtitle( 'Single copy genes trend')

ggsave('results/figures/SCG_geomean.pdf',  width =7 , height =5 )
