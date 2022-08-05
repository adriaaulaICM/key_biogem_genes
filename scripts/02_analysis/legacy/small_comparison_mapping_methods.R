library(tidyverse)
library(seqinr)


res.diamond <- read_rds('data/03_mapping/abtable_raw_famli.rds')
res.bowtie <-  read_rds('data/03_mapping/bowtie/lengthnorm.tbl.rds') 


names.correct <- data.frame( sample = as.character(1:84),
            sample_names = colnames(res.bowtie)[-c(1,2)])

res.diamond <- res.diamond %>% 
  left_join(names.correct, by = 'sample') %>% 
  select(-sample) %>% 
  rename( 'sample' = sample_names)

comparison <- res.bowtie %>% 
  pivot_longer(names_to = 'sample', cols =  -c(keggid, genename),
  values_to = 'count') %>% 
  left_join(res.diamond, by = c('genename' = 'sseqid', 'sample'))  %>% 
  rename('bowtie' = count.x, 
         'diamond' =  count.y)
   
 
# Big picture 
totals = data.frame(  bowtie = sum(comparison$bowtie),
                      diamond = sum(comparison$diamond, na.rm = T))

totals %>% 
  mutate( rel.bowtie = bowtie / diamond) %>% 
  knitr::kable( caption = 'Diamond protocol doubles down the amount of reads')

  

comparison %>% 
  mutate(diff = bowtie - diamond) %>% 
  filter(!is.na(diff)) %>% 
  pull(diff) %>% 
  skimr::skim()

only.bowtie <- comparison %>% 
  filter(is.na(diamond), bowtie > 0)

only.diamond <- comparison %>% 
  filter(bowtie == 0, diamond > 0)

ggplot(only.diamond, aes( x = diamond)) + 
  geom_histogram() + 
  scale_y_log10()

ggplot(only.bowtie, aes( x = bowtie)) + 
  geom_histogram() + 
  scale_y_log10()

fasta <- read.fasta('data/annotation/putative_genes/sel0vars.fasta')
distribution <- tibble(genename =  names(fasta),
                       gene.length = getLength(fasta)) 
descriptions <- readxl::read_xlsx('results/summary_charact/database_genes.xlsx')

compartop10 <- results %>%  
  pivot_longer(cols = c(bowtie,diamond),
               names_to = 'strategy', values_to = 'counts')  %>% 
  left_join( distribution, by = 'genename') %>% 
  left_join( descriptions, by = 'keggid') %>% 
  filter(!is.na(counts), counts >0) %>% 
  group_by(strategy,gene_name,gene.length, genename) %>% 
  summarise(total = sum(counts))  %>% 
  ungroup() %>% 
  group_by(strategy, gene_name) %>% 
  top_n(n = 15, wt = total) %>% 
  filter(!is.na(gene_name)) 



compartop10 %>% 
  group_by(genename) %>% 
  filter(n() == 2) %>% 
  ggplot( aes( fct_reorder(genename, total, .desc = TRUE), total)) + 
  geom_col(aes(fill = strategy)) + 
  # ggrepel::geom_label_repel( aes( label = gene.length), nudge_x = 0, size = 2) + 
  facet_wrap(~gene_name, scales = 'free') + 
  xlab( 'Gene number') + 
  ggtitle( 'Read counts per gene variants' , 
           subtitle = 'Label on top of the bar indicates the gene length of the variant') + 
  theme(axis.text.x = element_blank())

ggsave('results/figures/comparison_bowtie_diamond.pdf', height = 7, width = 7)

