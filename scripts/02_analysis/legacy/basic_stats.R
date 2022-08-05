library(tidyverse)

source('scripts/analysis/sourcefiles/backbone_params-graphs.R')
source('scripts/analysis/sourcefiles/params_biogem.R')


# ------------------------------------------------------------
ab.tbl <- read_rds('data/04_table_gen/scglengthnorm.tbl.rds') %>% 
  filter(!keggid %in% scg) 

colnames(ab.tbl) <-  str_replace_all(colnames(ab.tbl),
                                     c( "BL100413" = "BL100414",
                                        "BL110704" = "BL110705",
                                        "BL120518" = "BL120511",
                                        "BL131204" = "BL131215"))

# Dataset with the abundance counts 
ab.raw <- read_rds('data/04_table_gen/lengthnorm.tbl.rds') %>% 
  filter(!keggid %in% scg)  
# Descriptions and env variables 
descriptions <- readxl::read_xlsx('results/summary_charact/database_genes.xlsx')
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')



# Number of gene variants  ------------------------------------------------
nvariants.gene <- ab.tbl %>% 
  group_by(keggid) %>%
  tally() %>% 
  left_join(descriptions, by = 'keggid') %>% 
  select(cycle, gene_name, n ) %>% 
  arrange(cycle, gene_name, n)

gt::gt(nvariants.gene)

# Count distribution  -----------------------------------------------------
counts.gene <- ab.raw %>% 
  gather(key = 'sample', value = 'counts', -keggid, -genename) %>% 
  group_by(genename,keggid) %>% 
  summarize( total = sum(counts)) 
  


ggplot(counts.gene, aes(x = total)) + 
  geom_histogram() +  
  scale_x_log10(labels = function(x) format(x, scientific = FALSE)) + 
  ggtitle( 'Count distribution for gene variants', 
           subtitle = str_c('Due to the logscale, ',
                            sum(counts.gene$total == 0),
                            ' genes are not present'))

# For every gene, how many of them are with 0 counts? 
counts.gene %>%  
  left_join(descriptions, by = 'keggid') %>% 
  filter(total == 0) %>% 
  group_by( gene_name) %>% 
  tally()
