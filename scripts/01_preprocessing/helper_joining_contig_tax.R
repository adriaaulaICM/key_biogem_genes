library(tidyverse)

ab.alr <- read_rds('data/04_table_gen/abtblraw_all.rds')
tax <- readRDS('data/04_table_gen/taxonomy_subset.rds')

path <- '~/icm/tunnel/contigs_mmseqs2_only_lca/'
contigs.file <- list.files(path, pattern = 'lca')

tax.contig.df <- tibble(files = contigs.file) %>% 
  # read the tax files
  mutate( dataset = map(contigs.file,
                        ~read_tsv(str_c(path, .x), 
                                  col_names = c('contig_id', 'genomeid',
                                                'last_rank', 'rank_name',
                                                'n_protein_fragments', 'n_fragments_labelled', 'n_fragments_coinciding', 'support_percent',
                                                'taxonomy')
                                  ))) %>% 
  # obtain the sample origin for our contigs
  mutate(sam_origin = str_remove(files, pattern = '_lca.tsv')) %>% 
  # select(sam_origin, dataset) %>% 
  unnest(cols = dataset) %>% 
  # select(sam_origin, contig_id, last_rank, rank_name, taxonomy) %>% 
  # separating the content
  separate(col = taxonomy,
           into = c("species","genus","family",
                    "order","class","phylum","superkingdom"), sep = ";") 


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
          path = 'data/04_table_gen/taxonomy_contigwise.rds')

# Check and compare a little the results  ---------------------------------

tax.at.gene.lev <- tax %>% 
  filter(genename %in% tax.final.df$genename) 

tax.final.df %>% 
  filter(genename %in% coxL.gene) %>% 
  View()

tax.final.df %>% 
  filter(genename %in% gammas.psba) %>% 
  View('Tax contigs')


tax.at.gene.lev %>% 
  filter(genename %in% gammas.psba) %>% 
  View('Tax genes')

library(naniar)


vis_miss(tax.final.df, warn_large_data = F)
vis_miss(tax.at.gene.lev, warn_large_data = F)
