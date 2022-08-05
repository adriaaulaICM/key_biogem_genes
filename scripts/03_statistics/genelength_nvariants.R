library(tidyverse)
library(seqinr)

source('scripts/utils/backbone_params-graphs.r')
source('scripts/utils/backbone_functions.r')
source('scripts/utils/params_biogem.r')


# data import -------------------------------------------------------------

gene.count <- read_rds('data/intermediate_files/')
annot <- read_rds('data/intermediate_files/annotation.rds')
descriptions <- readxl::read_xlsx('results/gene_info/database_genes.xlsx')
fasta <- read.fasta(file = 'data/04_table_gen/genenames_all_more.fasta.gz')


gene.lengths <- tibble( genename = names(fasta),
                       length = getLength(fasta)) %>% 
  left_join(annot %>% select(genename, annotation,
                             gene_name, `length (aa)`,
                             `basic function`),
            by = 'genename')

View(gene.lengths)


ratio_genelength_nvariants <- gene.lengths %>% 
  group_by(annotation,gene_name) %>% 
  summarise(n_variants = n(),
            mean_gene_length = mean(length)) %>% 
  mutate( ratio = mean_gene_length / n_variants)


ggplot(ratio_genelength_nvariants %>% filter( !gene_name %in% c('chiA','nifH')),
       aes(ratio, gene_name)) + 
  geom_bar(stat= 'identity')
  

ratio_genelength_nvariants %>% filter(is.na(gene_name))
