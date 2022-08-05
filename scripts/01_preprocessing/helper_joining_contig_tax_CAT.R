library(tidyverse)

ab.alr <- read_rds('data/04_table_gen/abtblraw_all.rds')
tax <- readRDS('data/04_table_gen/taxonomy_subset.rds')

path <- '~/icm/tunnel/CAT_tax/'
contigs.file <- list.files(path)

tax.contig.df <- tibble(files = contigs.file) %>% 
  # read the tax files
  mutate( dataset = map(contigs.file,
                        ~read_tsv(str_c(path, .x), skip = 1,
                                  col_names =  c("contig_id", "classification",
                                                 "reason", "lineage",
                                                 "lineage scores", "superkingdom",
                                                 "phylum", "class", "order",
                                                 "family", "genus", "species") 
                        ))) %>% 
  # obtain the sample origin for our contigs
  mutate(sam_origin = str_remove(files, pattern = '.tax.txt')) %>% 
  select(sam_origin, dataset) %>% 
  unnest(cols = dataset) %>% 
  select(sam_origin, contig_id, phylum, class, order, family, genus, species)


# correspondence for each genename
contigs.df <- ab.alr %>%
  select(genename) %>% 
  distinct() %>% 
  mutate(contigs = str_extract(genename, 'BL[0-9]*\\.k127_[0-9]*')) %>%
  separate(col = contigs,
           into = c('sam_origin', 'contig_id'),
           sep = '\\.')

# Joining both dataset to obtain a complementary taxonomy 
tax.final.df <- contigs.df %>% 
  left_join(tax.contig.df, by = c('sam_origin', 'contig_id')) %>% 
  select(-sam_origin, -contig_id)


write_rds(x = tax.final.df,
          path = 'data/04_table_gen/taxonomy_contigwise_CAT.rds')

