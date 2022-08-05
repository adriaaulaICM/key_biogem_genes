library(tidyverse) 

# Relationship between one and the other db 
summar.char <- read_tsv('results/summary_charact/database_genes.txt')  %>% 
  filter(!is.na(cogs))

# Kegg annotation 
genes.kegg <- read_delim('data/annotation/kofam_search/putative_function.txt',
                       col_names = c('gene', 'keggid'), delim = ' ') 
# Results above 70% identity 
cog.output <- read_tsv('data/annotation/cog_recheck/cogannot_70id.txt',
                       col_names = F) %>% 
  separate(X2, into = c('gi', 'protein_id', 'ref', 'genbankid'), sep = "\\|")

# COG categories and protein  relationship
cog.ass <- read_csv('data/annotation/cog_recheck/cog2003-2014.csv',
                    col_names =  c("domain_id", "genome_name", "protein_id",
                                   "protein_length", "domain_start",
                                   "domain_end", "COG_id", "membership_class")) %>% 
  mutate( protein_id = as.character(protein_id))


# Keeping phosporous genes  -----------------------------------------------
wholelist.genes <- filter(genes.kegg, keggid %in% summar.char$keggid) %>% 
  pull(gene) %>% 
  unique()

output.joined <- cog.output %>% 
  left_join(cog.ass, by = 'protein_id')
  

process <- output.joined  %>% 
  select(X1,X3, COG_id) %>%  
  left_join(summar.char, by = c("COG_id" = 'cogs')) %>% 
  filter(!is.na(gene_name)) %>% 
  group_by(X1) %>% 
  top_n(n = 3, wt = X3) 

genes.keep <- process  %>%  
  select(X1, COG_id) %>% 
  distinct() %>% 
  group_by(X1) %>% 
  filter( n() <=  1) %>% 
  ungroup()
  
theonesout <- setdiff( wholelist.genes, genes.keep$X1)


# The genes which were a FP for phosporous genes
length(theonesout)
# The whole list of phosporous genes 
length(wholelist.genes)
# The whole cog output at 70% identity
length(unique(cog.output$X1))
# This ones summing the ones out should be wholelist.genes
length(genes.keep$X1)
  

write_tsv(data.frame( genes = theonesout),
          'data/annotation/cog_recheck/genes_out_putative.tsv', 
          col_names = F)  
length(setdiff( wholelist.genes, unique(cog.output$X1)) )

# Most of the wholelist genes is not present in the cog.output, which is ...
# weird? Most of them should be there also, since the annotation has taken place
# with different methods but 
