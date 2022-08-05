library(tidyverse)
library(lubridate)
source('scripts/analysis/sourcefiles/backbone_params-graphs.R')
source('scripts/analysis/sourcefiles/backbone_functions.R')
source('scripts/analysis/sourcefiles/params_biogem.R')

# I want to plot each gene with the 5th variants and the taxonomy to the lowest 
# level for each variant 

# Once i do that, I might explore with the most abundant variants if they cluster 
# based in taxonomy



theme_set(theme_bw(base_size = 20))

# Data import -------------------------------------------------------------
ab.alr.filt <- read_rds('data/04_table_gen/abtblalr_filt0rem_25sams.rds') 
ab.raw.filt <- read_rds('data/04_table_gen/abtblraw_filt0rem_25sams.rds') 
ab.rawcount.tbl <- read_rds('data/03_mapping/abtable_raw_famli.rds')  %>% 
  rename( 'genename' = sseqid)
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')
annot <- read_rds('data/intermediate_files/annotation.rds')
  
tax <- read_tsv('data/04_taxgenes/taxonomyresult.tsv',
                col_names = c('genename', 'genomeid', 'last_rank',
                              'rank_name', 'taxonomy')) %>% 
  separate(col = taxonomy,
           into = c("species","genus","family",
                    "order","class","phylum","superkingdom"), sep = ";") %>% 
  mutate( label = case_when(is.na(phylum) ~ rank_name,
                                 !is.na(phylum) ~ str_c(phylum, family,sep = ", "),
                                 TRUE ~ rank_name))

lomb <- readRDS('data/intermediate_files/genes_seasonal.rds') 

ab.mega <- ab.alr.filt %>% 
  left_join(annot, by = 'genename') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) %>% 
  left_join(lomb %>% select(genename, peak), by = 'genename') %>% 
  mutate( seasonal = ifelse( !is.na(peak), 'seasonal', 'not seasonal'),
          count = ifelse(count < 0, count + 0,count))
  


# Variant calling (top 5) -------------------------------------------------

top5_genvars <- ab.mega %>% 
  filter(seasonal == 'seasonal') %>% 
  group_by(annotation, genename) %>% 
  summarise( mean.sam = mean(count)) %>% 
  top_n(n = 5, wt = mean.sam)  %>% 
  arrange(annotation,-mean.sam)  %>% 
  mutate(varnum = factor(genename,
                         levels = unique(genename),
                         labels = str_c('Var. ', 1:length(genename)) 
                         ))  %>% 
  ungroup()

write_rds( top5_genvars, 'data/intermediate_files/top5vars.rds')

pal_varcall <- c("#00478E", "#00BFB3", "#DBE442", "#FFB1BB", "#E24585")
names(pal_varcall) <- str_c('Var. ', 1:5)

plot_variants_trends <- function(genenam){
  ab.topsub5 <- ab.mega %>% 
    filter(genename %in% top5_genvars$genename) %>% 
    left_join(top5_genvars %>% select(genename, varnum), by = 'genename') %>% 
    filter(gene_name == genenam) 
  
  ab.topsub5 %>% 
    ggplot(aes(month, count, color = varnum) ) + 
    geom_jitter(shape = 21, show.legend = F) + 
    stat_smooth(aes(x = month,
                    group = varnum), # continuous x-axis
                method = "lm",
                formula = y ~ poly(x, 3), se = F, size = 2, show.legend = F) + 
    ylab('Log ratio (gene / SCG)') + 
    scale_x_month +
    scale_color_manual(values = pal_varcall) + 
    scale_fill_manual(values = pal_varcall) + 
    theme_minimal(base_size = 19)
  
}

  

plot_varlabels <- function(this){
  
  taxinfo <-  tax %>%
    filter(genename %in% top5_genvars$genename) %>% 
    left_join(annot, by = 'genename') %>% 
    filter(gene_name == this) %>% 
    as_tibble() %>%
    select(genename, label) %>%
    left_join(top5_genvars %>% select(genename, varnum), by = 'genename') %>%
    distinct(genename, .keep_all = T) %>%
    arrange(varnum) %>%
    mutate( y = 1:length(varnum) ) 
  
  ggplot(taxinfo, aes(x = 1, y = y, label = label, fill = varnum)) + 
    geom_label(show.legend = F,
               size = 4,
               color = 'white',
               vjust = "inward", hjust = "inward") + 
    scale_color_manual(values = pal_varcall) + 
    scale_fill_manual(values = pal_varcall) + 
    theme_void(base_size = 19)
  
}

label5plot <- ab.mega %>% pull(gene_name) %>% unique() %>% 
  set_names(map(., ~plot_varlabels(.x)), .) 

top5plot <- ab.mega %>% pull(gene_name) %>% unique() %>% 
  set_names(map(., ~plot_variants_trends(.x)), .) 

write_rds(top5plot, 'data/intermediate_files/variant_callingplot.rds')
write_rds(label5plot, 'data/intermediate_files/variant_labelplot.rds')


save_mult_plot(plots =  top5plot, 
              vec_files = names(top5plot), 
              common = 'results/figures/variant_calling/var_call_tax_',
              endfile = '.pdf', width = 9, height = 6)
  

# Seasonality information  ------------------------------------------------

# we want to see the maximum abundance month for each variant 
# and plot it in a descriptive way


# Summary table  ----------------------------------------------------------
total.counts <- ab.rawcount.tbl %>% 
  left_join(annot, by = 'genename') %>% 
  select(gene_name, genename) %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarise(total = n()) 


month.max <- ab.mega %>% 
  filter(seasonal == 'seasonal') %>% 
  group_by(gene_name, genename, month) %>% 
  summarize( medianmax = mean(count)) %>% 
  filter(medianmax == max(medianmax)) %>% 
  group_by(gene_name,month) %>% 
  summarise( count = n())  %>% 
  group_by(gene_name) %>% 
  top_n( n = 1, wt = count)

month.max



# Month maximum variant graph ---------------------------------------------

summary.maxm <- ab.mega %>% 
  filter(seasonal == 'seasonal') %>% 
  group_by(gene_name,genename, peak, month) %>% 
  summarize( medianmax = mean(count)) %>% 
  group_by(genename) %>% 
  top_n(wt = medianmax, n = 2) %>% 
  # filter(medianmax == max(medianmax))  %>% 
  ungroup() %>% 
  rename( 'month.max' = month) %>% 
  mutate( month.max = factor(month.max, levels = month.num,
                             labels = str_to_title(month.order)))


plot_maxtrend_distr <- function(gene){
   
  dataset <- summary.maxm  %>%
    left_join(top5_genvars, by = 'genename') %>%
    filter(gene_name == gene)
  
  line.data <-   dataset %>%
    filter(!is.na(varnum)) %>%
    group_by(genename) %>%
    mutate( month.max.num = as.numeric(month.max)) %>%
    mutate( diff = max(month.max.num) - min(month.max.num)) %>% 
    filter(diff < 4)
  
  ggplot(dataset , 
         aes(x = fct_rev(month.max),
             y = peak)) +
    geom_point(alpha = 0.5, size = 3)  +
    geom_point(data = dataset %>% filter(!is.na(varnum)),
               aes(color = varnum), size = 3.1, show.legend = F)  +
    geom_line(data = line.data,
               aes(color = varnum, group = varnum), lineend = "round",
              size = 3, show.legend = F)  +
    ylab('Recurrence index') + 
    xlab('Month maximum') + 
    scale_color_manual(values = pal_varcall) + 
    scale_x_discrete(drop = FALSE) + 
    theme_minimal() + 
    coord_flip()
  
}

plot_bartrend <- function(gene){
  
  dataset <- summary.maxm  %>%
    left_join(top5_genvars, by = 'genename') %>%
    filter(gene_name == gene)
  
  ggplot(dataset,
         aes(x = (month.max))) +
    geom_bar(aes(fill = varnum),stat = 'count', show.legend = F) + 
    xlab('Month maximum') + 
    ylab('Count') + 
    scale_x_discrete(drop = FALSE) + 
    scale_fill_manual(values = pal_varcall, na.value = 'grey')  +
    theme_minimal(base_size = 15)
  
}

bartrendsplot <- summary.maxm %>% pull(gene_name) %>% unique() %>% 
  set_names(map(., ~plot_bartrend(.x)), .) 

maxtrendsplot <- summary.maxm %>% pull(gene_name) %>% unique() %>% 
  set_names(map(., ~plot_maxtrend_distr(.x)), .) 


write_rds(maxtrendsplot, 'data/intermediate_files/maxtrendsplots.rds')

save_mult_plot(plots =  maxtrendsplot, 
              vec_files = names(maxtrendsplot), 
              common = 'results/figures/variant_calling/seasonal_tendencies_',
              endfile = '.pdf', width = 9, height = 6)
  
write_rds(bartrendsplot, 'data/intermediate_files/bartrendsplots.rds')

save_mult_plot(plots =  bartrendsplot, 
              vec_files = names(bartrendsplot), 
              common = 'results/figures/variant_calling/monthmax_',
              endfile = '.pdf', width = 9, height = 6)
  

heatmap.df <- summary.maxm %>% 
  group_by(gene_name, month.max) %>%
  summarize( count = n()) %>% 
  ungroup()


this <- pivot_wider(heatmap.df,
            names_from = gene_name,
            values_from = count,
            values_fill = list(count = 0)) %>% 
  column_to_rownames(var = 'month.max')

library(pheatmap)
pheatmap(this, cluster_rows = F, cluster_cols = T)

