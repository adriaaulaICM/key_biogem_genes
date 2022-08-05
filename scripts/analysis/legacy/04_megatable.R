library(tidyverse)



# Data import -------------------------------------------------------------
ab.alr.filt <- read_rds('data/04_table_gen/abtblalr_filt0rem_25sams.rds') 
ab.agg <- read_rds('data/04_table_gen/alragg.tbl.rds') 
descriptions <- readxl::read_xlsx('results/summary_charact/database_genes.xlsx')
annot <- read_rds('data/intermediate_files/annotation.rds')
tax <- read_tsv('data/04_taxgenes/taxonomyresult.tsv',
                col_names = c('genename', 'genomeid', 'last_rank',
                              'rank_name', 'taxonomy')) %>% 
  separate(col = taxonomy,
           into = c("species","genus","family",
                    "order","class","phylum","superkingdom"), sep = ";")

envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')

lomb <- readRDS('data/intermediate_files/genes_seasonal.rds')

ab.mega <- ab.alr.filt %>% 
  left_join(annot, by = 'genename') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) %>% 
  left_join(lomb %>% select(genename, peak), by = 'genename') %>% 
  mutate( seasonal = ifelse( !is.na(peak), 'seasonal', 'not seasonal'),
          count = ifelse(count < 0, count + 0,count))

eval <- ab.mega %>% 
  select(genename,cycle, gene_name, seasonal) %>% 
  distinct(genename, .keep_all = T)  %>% 
  group_by(gene_name,cycle) %>% 
  summarize(evaluated = n()) %>% 
  ungroup()

seasonalvar <- ab.mega %>% 
  select(genename, gene_name, seasonal) %>% 
  distinct(genename, .keep_all = T)  %>% 
  filter(seasonal == 'seasonal') %>%
  group_by(gene_name) %>% 
  summarise( nseasonal = n()) 

table <- eval %>% 
  left_join(seasonalvar, by = 'gene_name') %>% 
  left_join(total.counts, by = 'gene_name') %>% 
  select(gene_name, cycle, total, evaluated, nseasonal) %>% 
  mutate( relab_seasonal = ((nseasonal * 100) / evaluated) %>% round(digits = 1))
  


# Preparing the first datasets --------------------------------------------

# Seasonality of  the aggregated counts 
sea.agg <- read_tsv('results/stats/seasonality_agg.tsv') %>% 
  mutate(agg_seasonality  = ifelse( pval <= 0.01, 'yes', 'no')) %>% 
  select(gene_name, agg_seasonality)

# Seasonality of the variants 
lomb.filt <- lomb %>% 
  filter(peak >= 10) %>% 
  filter(qval <= 0.01)

# Season that the aggregated values present the max value if seasonal
maxseason <- ab.agg %>% 
  left_join(descriptions %>% select(annotation, gene_name, cycle),
            by = 'annotation') %>% 
  filter(!str_detect(string = cycle, pattern = 'scg')) %>% 
  select(-annotation, -cycle) %>% 
  pivot_longer(names_to = 'sample',
               values_to = 'logratio',
               cols = -c(gene_name)) %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) %>% 
  group_by(gene_name,season, month) %>% 
  summarize( mean.trend = median(logratio)) %>% 
  left_join(sea.agg, by = 'gene_name') %>% 
  group_by(gene_name) %>% 
  top_n(n = 2, wt = mean.trend) %>% 
  group_by(gene_name, agg_seasonality) %>% 
  summarise(max_trend_agg  = str_c(unique(season), collapse = ','))  %>% 
  mutate( max_trend_agg =  ifelse(agg_seasonality == 'yes', max_trend_agg, '-')) %>% 
  mutate( max_trend_agg =  ifelse(max_trend_agg == 'autumn,summer',
                                  'summer,autumn', max_trend_agg)) 
  
# recurreence index value median 
peak.mean <- lomb.filt %>% 
  left_join(annot, by = 'genename') %>% 
  group_by(gene_name) %>% 
  summarize(mean.peak = median(peak) %>% round(digits = 1)) 

predtax <- tax %>% 
  distinct(genename, .keep_all = T) %>% 
  filter(genename %in% lomb.filt$genename)  %>% 
  left_join(annot %>% select(genename, gene_name), by = 'genename') %>% 
  group_by(gene_name, phylum, class,order,family) %>% 
  summarize(count = n()) %>% 
  group_by(gene_name) %>% 
  filter(sum(count) > 5) %>% 
  mutate( relab = count / sum(count)) %>% 
  ungroup()
  

classpred <- predtax %>% 
  group_by(gene_name, class) %>% 
  summarize( relab = sum(relab) * 100 %>% round(digits = 2)) %>% 
  filter(relab >= 25) %>% 
  arrange(gene_name, -relab)  %>% 
  mutate(class_pred = str_c( class, str_c(round(relab, digits = 1), '%'),
                             sep = ',', collapse = ' // ')) %>% 
  distinct(gene_name, class_pred)

fampred <- predtax %>% 
  group_by(gene_name, family) %>% 
  summarize( relab = sum(relab) * 100 %>% round(digits = 2)) %>% 
  filter(relab >= 25) %>% 
  arrange(gene_name, -relab)  %>% 
  mutate(fam_pred = str_c( family, str_c(round(relab, digits = 1), '%'),
                           sep = ',', collapse = ' // ')) %>% 
  distinct(gene_name, fam_pred )

# Major taxonomic group (25% min relab)
predmajor <- classpred %>% 
  left_join(fampred, by = 'gene_name')

# Clusters of seasonal variants 
clusters_season <- read_tsv('results/stats/possible_clusters.txt')

# Joining all together ----------------------------------------------------

finaltable <- table %>% 
  left_join(descriptions %>% select(gene_name, `basic function`),
            by = 'gene_name') %>% 
  left_join(maxseason, by = 'gene_name') %>% 
  left_join(peak.mean, by = 'gene_name') %>% 
  left_join(predmajor, by = 'gene_name') %>% 
  left_join(clusters_season, by = 'gene_name') %>% 
  select(gene_name, cycle, `basic function`, total, evaluated, agg_seasonality,
         max_trend_agg, nseasonal, relab_seasonal, mean.peak, class_pred, fam_pred,
         n_clusters, season_maximums) 
  
write_rds(finaltable, 'data/intermediate_files/finaltable.rds')

  