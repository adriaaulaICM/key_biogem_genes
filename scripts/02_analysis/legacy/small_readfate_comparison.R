library(tidyverse)

pufms_ab <- read_rds('data/04_table_gen/abtblalr_all.rds') %>%
  left_join(annot, by = 'genename') %>%
  filter(gene_name == 'pufM') %>%
  distinct(genename, .keep_all = T)  %>% 
  pull(genename)

pufms_annot <- read_rds('data/intermediate_files/annotation.rds') %>% 
  filter(gene_name == 'pufM') %>%
  distinct(genename, .keep_all = T)  %>% 
  pull(genename)
  
setdiff(pufms_annot, pufms_ab) 


# From these 73 variants, I take one and one sample. Specifiacally I have taken
# the sample from where the variant has come (and therefore hypotetically there 
# will be more reads than in other

# I go to the diamond mapping of the reads and I retreive all the results 
# that present this variant. I do the same with the final result coming from famli
# in which each read is having a hard fate to one variant 

reads.matching.variant <- read_tsv('data/intermediate_files/reads_in_pufm_BL150609.k127_600175.mgm_788968.txt',
                       col_names = c('qseqid','sseqid','pident','length','mismatch',
                                     'gapopen','qstart','qend','sstart',
                                     'send','evalue','bitscore','qlen','slen')) %>% 
  distinct(qseqid, sseqid, .keep_all = T)


reads.fate <- read_tsv('data/intermediate_files/reads_fate.txt',
                       col_names = c('qseqid','sseqid','pident','length','mismatch',
                                     'gapopen','qstart','qend','sstart',
                                     'send','evalue','bitscore','qlen','slen')) %>% 
  distinct(qseqid, sseqid, .keep_all = T)

reads.variant.good <- read_tsv('data/intermediate_files/reads_BL110517.k127_1088171.mgm_1360378_final.txt',
                       col_names = c('qseqid','sseqid','pident','length','mismatch',
                                     'gapopen','qstart','qend','sstart',
                                     'send','evalue','bitscore','qlen','slen')) %>% 
  distinct(qseqid, sseqid, .keep_all = T)

bothdf <- bind_rows(reads.matching.variant, reads.fate)

# read fate if we only would consider the maximum identity
readstop <- bothdf %>% 
  group_by(qseqid) %>% 
  top_n(n = 1, wt = pident) 

readstop  %>% 
  pull(sseqid)  %>% 
  table()

therealdeal <- readstop %>%  
  ungroup() %>% 
  group_by(qseqid) %>% 
  filter(n() != 2) %>% 
  filter(sseqid == 'BL150609.k127_600175.mgm_788968')

# real read fate 
reads.fate %>% 
  pull(sseqid) %>% 
  table()

library(patchwork)

plot_identity_gene <- function(df){
  
  ggplot(df, aes( y = qseqid, color = pident)) +
    geom_segment(aes(x = sstart,xend =send, yend = qseqid)) + 
    facet_wrap(~sseqid) + 
    theme(axis.text.y = element_blank(), axis.ticks = element_blank())
  
}

plot_shared_reads <- function(df){
  
  ggplot(df, aes( y = qseqid, color = topsies)) +
    geom_segment(aes(x = sstart,xend =send, yend = qseqid)) + 
    facet_wrap(~sseqid) + 
    theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + 
    scale_color_manual(values = c('shared' = 'firebrick', 'nonshared' = 'grey'))
}

prepare_readdf <- function(df, identity = 95, compar.reads = readstop$qseqid){
  df %>%  
    filter(pident >= identity) %>% 
    arrange(sseqid, sstart) %>% 
    mutate( topsies = ifelse(qseqid %in% compar.reads, 
                             'shared', 'nonshared')) %>% 
    mutate(qseqid = fct_inorder(qseqid))   
}



prepare_readdf(reads.variant.good) %>% 
  plot_identity_gene()  +
prepare_readdf(reads.variant.good) %>% 
  plot_shared_reads() +
prepare_readdf(reads.matching.variant) %>% 
  plot_identity_gene() +
prepare_readdf(reads.matching.variant) %>% 
  plot_shared_reads()
  
prepare_readdf(reads.variant.good, identity = 100) %>% 
  plot_identity_gene()  +
prepare_readdf(reads.variant.good, identity = 100) %>% 
  plot_shared_reads() +
prepare_readdf(reads.matching.variant, identity = 100) %>% 
  plot_identity_gene() +
prepare_readdf(reads.matching.variant, identity = 100) %>% 
  plot_shared_reads()
  
prepare_readdf(reads.variant.good, identity = 100) %>% 
  plot_identity_gene()  +
prepare_readdf(reads.variant.good, identity = 100) %>% 
  plot_shared_reads() +
prepare_readdf(reads.matching.variant, identity = 100) %>% 
  plot_identity_gene() +
prepare_readdf(reads.matching.variant, identity = 100) %>% 
  plot_shared_reads()
  

prepare_readdf(reads.matching.variant) %>% 
  plot_identity_gene() +
prepare_readdf(reads.matching.variant) %>% 
  mutate( topsies = ifelse(qseqid %in% therealdeal$qseqid, 'unique', 'non')) %>% 
  plot_shared_reads() + 
  scale_color_manual(values = c('gray75','black'))
 

