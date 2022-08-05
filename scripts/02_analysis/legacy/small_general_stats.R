library(tidyverse)
library(seqinr)

# a small script for long calculations from our beautiful dataset

# Frequency in diamond and bowtie  ----------------------------------------

ab.rawcount.tbl <- read_rds('data/03_mapping/abtable_raw_famli.rds')  %>% 
  rename( 'genename' = sseqid)

ab.rawcount.tbl %>% 
  group_by(genename) %>% 
  summarise(ocurrence =  sum(count > 0) ) %>% 
  pull(ocurrence) %>% 
  table() %>% 
  data.frame() %>% 
  write_tsv('results/stats/freqdist_per_gene_diamond.tsv')

ab.rawcount.bowtie <- read_rds('data/03_mapping/bowtie/lengthnorm.tbl.rds') %>% 
  pivot_longer(names_to = 'sample', values_to = 'count', cols = -c(genename, keggid))

ab.rawcount.tbl %>% 
  group_by(genename) %>% 
  summarise(ocurrence =  sum(count > 0) ) %>% 
  pull(ocurrence) %>% 
  table() %>% 
  data.frame() %>% 
  write_tsv('results/stats/freqdist_per_gene_diamond.tsv')
ab.rawcount.bowtie <- read_rds('data/03_mapping/bowtie/lengthnorm.tbl.rds') %>% 
  pivot_longer(names_to = 'sample', values_to = 'count', cols = -c(genename, keggid))

shit <- bind_rows( ab.rawcount.tbl  %>% 
             group_by(genename) %>% 
             summarise(ocurrence =  sum(count > 0), abundance = sum(count)) ,
           
           ab.rawcount.bowtie  %>% 
             group_by(genename) %>% 
             summarise(ocurrence =  sum(count > 0), abundance = sum(count)), 
           .id = 'method') %>% 
  filter( ocurrence > 20)


shit %>% 
  ggplot( aes(ocurrence, abundance, color = method))  + 
  geom_point(alpha = 0.7) + 
  facet_wrap(~method) + 
  scale_y_log10() 

freq.diamond <- read_tsv('results/stats/freqdist_per_gene_diamond.tsv')
freq.bowtie <- read_tsv('results/stats/freqdist_per_gene_bowtie.tsv')

freq.compar <- bind_rows(freq.diamond, freq.bowtie, .id = 'method')
colnames(freq.compar) <- c('method', 'ocurrence', 'freq')

ggplot(freq.compar %>% filter(ocurrence > 10),
       aes(x = ocurrence, y = freq, fill = method)) + 
  geom_bar(stat = 'identity') + 
  facet_wrap(~method)

freq.compar %>% 
  group_by(method) %>% 
  summarise(total_genes = sum(freq))


# Gene length stat  -------------------------------------------------------

fasta <- read.fasta('data/annotation/putative_genes/final_min250.fna')
lengthdist <- tibble(genename =  names(fasta),
                       gene.length = getLength(fasta)) 

write_tsv(lengthdist, 'results/stats/length_dist.tsv')