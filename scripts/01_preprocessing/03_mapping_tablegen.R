library(tidyverse) 

dir.famli <- '~/scratch/biogem_7years_bigoutputs/mapping_parse_famli/'
files_aln <- list.files(dir.famli, pattern = 'aln', full.names = T)

count_reads <- function(df){
  
  df %>% 
    distinct(qseqid, .keep_all = T) %>%
    filter(pident >= 95) %>% 
    group_by(sseqid) %>% 
    summarise( count  = n()) %>% 
    return()
  
}

res <-  files_aln %>% 
  map(~read_tsv(.x, 
                col_names = c("qseqid","sseqid", 
                                 "pident","length","mismatch","gapopen","qstart",
                                 "qend","sstart","send","evalue", "bitscore",
                                 "qlen","slen"))) %>% 
  map(~count_reads(.x))

names(res) <- list.files(dir.famli, pattern = '.aln')  %>% 
  str_replace(., ".aln", "")

allsamplesdf <- bind_rows(res, .id = 'sample')

write_rds(allsamplesdf, path = 'data/03_mapping/abtable_raw_famli.rds')

