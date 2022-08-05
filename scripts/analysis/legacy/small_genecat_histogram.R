library(seqinr)

# This script generates an histogram showing the  various genes predicted, the 
# distribution of these along the various gene lengths 

fasta <- read.fasta('~/icm/tunnel/nuc_genes.fna')
head(names(fasta))

distribution <- tibble(names =  names(fasta),
                       gene.length = getlength(fasta)) %>% 
  mutate(origin = ifelse( str_detect(names, 'mgm'), 'mgmark', 'prodigal')) 
  


ggplot(distribution, aes(x = gene.length, fill = origin)) + 
  geom_histogram()+ 
  geom_vline(xintercept = 310, linetype = 2, color = 'black')  + 
  scale_x_log10() + 
  theme(legend.position = c(0.8, 0.8), legend.background = element_blank()) + 
  scale_fill_fivethirtyeight() + 
  ggtitle('Gene length distribution',
          subtitle = 'Mgmark predictions present many genes ~300 pb')


