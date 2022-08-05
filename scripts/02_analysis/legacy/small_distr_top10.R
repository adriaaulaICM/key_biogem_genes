library(tidyverse)
source('scripts/analysis/sourcefiles/backbone_params-graphs.R')

#This script prints the distribution of the top10 variants to look if the tendency 
#of abundance is shared between the top or instead there is a rank curve in which 
#the most abundant is really super abundant
#We will also check if this is due to gene length or not

# It will also check how many genes are present in all the timeseries and how 
# many are only in a small amount of samples 

ab.raw <- read_rds('data/04_table_gen/abtable_raw_0rem_25samrem_noscg.rds') 
descriptions <- readxl::read_xlsx('results/summary_charact/database_genes.xlsx')
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')

megafile.agg <-  ab.raw  %>% 
  gather(key = 'sample', value = 'count', -keggid, -genename)

onlytops <- megafile.agg %>% 
  group_by(genename) %>% 
  filter(sum(count) > 1000)

sel.vars <- onlytops %>% 
  group_by(keggid, genename) %>% 
  summarise( mean.sam = mean(count)) %>% 
  group_by( keggid) %>% 
  top_n(n = 8, wt = mean.sam) %>% 
  pull(genename)


genvar.df <- megafile.agg %>% 
  left_join(descriptions, by = 'keggid') %>% 
  filter( genename %in%  sel.vars) 

# Adding up the gene length  ----------------------------------------------

library(seqinr)
fasta <- read.fasta('data/annotation/putative_genes/sel0vars.fasta')
distribution <- tibble(genename =  names(fasta),
                       gene.length = getLength(fasta)) 

plot <- genvar.df %>%  
  left_join( distribution, by = 'genename') %>% 
  group_by(gene_name,gene.length, genename) %>% 
  summarise(total = sum(count))  %>% 
  filter(!is.na(gene_name)) %>% 
  ggplot( aes( fct_reorder(genename, total, .desc = TRUE), total)) + 
  geom_col() + 
  ggrepel::geom_label_repel( aes( label = gene.length), nudge_x = 0, size = 2) + 
  facet_wrap(~gene_name, scales = 'free') + 
  xlab( 'Gene number') + 
  lil.strip + 
  ggtitle( 'Read counts per gene variants' , 
           subtitle = 'Label on top of the bar indicates the gene length of the variant') + 
  theme(axis.text.x = element_blank())

ggsave(plot = plot, 'results/figures/suppl_abundance_genelenghtdistr.pdf',
       width = 12, height = 9)



