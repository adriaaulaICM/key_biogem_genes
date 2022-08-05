library(tidyverse)
library(seqinr)
source('scripts/utils/backbone_params-graphs.R')

# This script imports the data, generates some general statistics regarding the 
# amount of variants and zeroes and filters and corrects some wrong things in
# the document 

ab.rawcount.tbl <- read_rds('data/03_mapping/abtable_raw_famli_reduced.rds')  %>% 
  rename( 'genename' = sseqid)

#I have decided that I will kill 4 samples due to a  small amount of reads
#proceeding into doing it 
ab.rawcount.tbl <- ab.rawcount.tbl %>% filter(!sample %in% c('BL100322',
                                                             'BL100622',
                                                             'BL101116',
                                                             'BL090804')) 

# Import annotation ------------------------------------------------------
kegg <- read_delim('data/annotation/kofam_search/putative_function.txt', 
                      col_names = c('genename', 'annotation'), delim = ' ')

pr <-  read_delim('data/annotation/PR_search/putative_pr.txt',
                          col_names = c('genename', 'annotation'), delim = ' ')

cog <- read_tsv('data/annotation/cog_annot/putative_function.txt',
                          col_names = c('genename', 'annotation'))

descriptions <- readxl::read_xlsx('results/gene_info/database_genes.xlsx')

annot <- bind_rows(kegg,pr,cog)  %>% 
  # Some SCG and narB/nasA have double annotations, we have to erase them 
  distinct(genename, .keep_all = TRUE)  %>% 
  left_join(descriptions, by = 'annotation')


# Filter out genes considered unusable ------------------------------------

# stats 
discarded.genes <- annot %>% 
  filter(gene_name %in% c('hao', 'nifH', 'chiA', 'amoA')) %>% 
  pull(genename) %>% 
  unique()

ab.rawcount.tbl %>% 
  filter(genename %in% discarded.genes) %>% 
  left_join(annot, by = 'genename') %>% 
  group_by(gene_name) %>% 
  summarize(  min = min(count),
              mean = mean(count),
              max = max(count))

ab.rawcount.tbl %>% 
  filter(genename %in% discarded.genes) %>% 
  left_join(annot, by = 'genename') %>% 
  group_by(gene_name) %>% 
  summarize(total = sum(count))

ab.rawcount.tbl %>% 
  filter(genename %in% discarded.genes) %>% 
  left_join(annot, by = 'genename') %>% 
  distinct(genename, gene_name) %>% 
  group_by(gene_name) %>%
  summarize( count = n())



# we will erase from the dataset nifH, hao, dddP and chiA 
annot <- annot %>% 
  filter(!gene_name %in% c('hao', 'nifH', 'chiA', 'dddP')) 


# Import specific differentations  ----------------------------------------
dir.pr <- 'data/annotation/PR_search/motif_finding/'
pr.spec <- read_tsv(str_c(dir.pr,'spec_tunning_aa_checked.tsv')) %>% 
  select(genename, pr.description)
  
codh1 <- read_tsv('data/annotation/kofam_search/codh1_class.tsv', 
                       col_names = c('genename', 'codh_type'))

annot <- annot %>% 
  left_join(pr.spec, by = 'genename') %>% 
  mutate( gene_name = ifelse( gene_name == 'PR',
                               pr.description, gene_name)) %>% 
  filter(gene_name != 'coxL' | gene_name == 'coxL' & genename %in% codh1$genename) %>% 
  filter(gene_name != 'PR')

write_rds(annot, 'data/intermediate_files/annotation_all.rds')
# SCG  --------------------------------------------------------------------

# We will generate a separate file for the scg genes 
scg <- c("K06942", "K01889", "K01887", "K01875", "K01883",
         "K01869", "K01873", "K01409", "K03106", "K03110")

# I have selected an specific subset of SCG gene variants to perform the alr
# normalization 
scg_genes <- annot %>% filter(annotation %in% scg, 
                 !gene_name %in% c('LeuS','ValS')) %>% 
  pull(genename) 

count.scg <- ab.rawcount.tbl %>% 
  filter(genename %in% scg_genes) %>% 
  left_join(annot, by = 'genename')


# Write small annot  ------------------------------------------------------

#...without the SCG basically 
annot %>% 
  filter(!annotation %in% scg) %>% 
  write_rds('data/intermediate_files/annotation.rds')

# Geometric mean SCG ------------------------------------------------------

#before all calculations we will stick to only the genes from the annotation
ab.rawcount.tbl <- ab.rawcount.tbl %>% 
  filter(genename %in% unique(annot$genename))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geomean.sample <- count.scg %>% 
  group_by(sample, gene_name) %>% 
  summarize( count.kegg = sum(count))  %>% 
  group_by(sample) %>% 
  summarise( geomean = gm_mean(count.kegg)) 
  

ab.wider <- ab.rawcount.tbl %>% 
  pivot_wider(names_from = sample,
              values_from = count,
              values_fill = list(count = 0)) 

genenames <- ab.wider %>% select(genename)

ab.alr.tbl  <- map2(ab.wider %>% select(-genename),
                    geomean.sample$geomean, ~ log2((.x + 0.01)  / .y)) %>% 
  bind_rows() %>% 
  bind_cols(genenames, . )


ab.agg.wider <-  ab.rawcount.tbl %>% 
  left_join(annot, by = 'genename') %>% 
  group_by(gene_name, sample) %>% 
  summarize( total = sum(count)) %>% 
  pivot_wider(names_from = sample,
              values_from = total,
              values_fill = list(total = 0))  %>% 
  filter(!is.na(gene_name)) %>% 
  ungroup()

ab.agg.pr.wider <- ab.agg.wider %>%  
  filter(str_detect(gene_name, 'PR ')) 


pr.agg <- map_df(ab.agg.pr.wider[,-1], ~sum(.x)) 


genenames <- ab.agg.wider %>% select(gene_name)
ab.alr.keggagg <- map2(ab.agg.wider %>% select(-gene_name),
                       geomean.sample$geomean, ~log2((.x + 0.01) / .y)) %>% 
  bind_rows() %>% 
  bind_cols(genenames, .)

ab.alr.keggagg.nolog <- map2(ab.agg.wider %>% select(-gene_name),
                       geomean.sample$geomean, ~(.x  / .y)) %>% 
  bind_rows() %>% 
  bind_cols(genenames, .)

ab.alr.pr <- map2(pr.agg, 
                  geomean.sample$geomean, ~log2((.x + 0.01) / .y)) %>% 
  bind_rows() %>% 
  mutate( gene_name = 'PR') %>% 
  select(gene_name, everything())

ab.alr.keggagg <- bind_rows(ab.alr.keggagg, ab.alr.pr)

# Keep all variants with > 0 reads  ---------------------------------------
variants_abov0 <- ab.rawcount.tbl %>% 
  group_by(genename) %>% 
  filter( sum(count) > 0)

# Keep all variants happening in at least 10% samples ---------------------
variants_20samples <- ab.rawcount.tbl %>% 
  group_by(genename) %>% 
  filter( sum(count > 0) >= 8) # 80 samples, 8 its  a 10%

count.dist <- ab.rawcount.tbl %>% 
  group_by(genename) %>% 
  summarize(total.read = sum(count),
            ocurrence =  sum(count > 0) ) 

count.plot <- ggplot(count.dist, aes(ocurrence)) + 
  geom_histogram() + 
  scale_y_log10() + 
  ylab('Count') + 
  xlab('Ocurrence (n. of samples)')

ggsave(plot = count.plot,
       'results/figures/ocurrence_gene_distr.pdf', width = 7, height = 5)

dist.plot <- ggplot(count.dist %>%
                      filter(ocurrence >= 20),
                    aes(total.read, ocurrence)) + 
  geom_hex()  + 
  scale_x_log10() +
  xlab('Total read') + 
  ylab('Ocurrence (n. of samples)')

ggsave(plot = dist.plot,
       'results/figures/ocurrence_ab_gene_distr.pdf',
       width = 7, height = 5)

# Intersection of both datasts 
vars.keep <- intersect(unique(variants_abov0$genename), 
                       unique(variants_20samples$genename))


# Data import and creation of datasets  -----------------------------------

correct_samnames <- function(df){
  
  df$sample <-  str_replace_all(df$sample,
                                c( "BL100413" = "BL100414",
                                   "BL110704" = "BL110705",
                                   "BL120518" = "BL120511",
                                   "BL131204" = "BL131215"))
  return(df)
}

save_abtbl_rds <- function(df, keep, scg_genes, name){
  
  df %>% 
    correct_samnames() %>% 
    filter( genename %in% keep)  %>% 
    filter( !genename  %in% scg_genes) %>% 
    write_rds(str_c('data/04_table_gen/abtbl', name, '.rds'))
  
  df %>% 
    correct_samnames() %>% 
    filter( genename  %in% scg_genes) %>% 
    write_rds(str_c('data/04_table_gen/abtbl', name, 'onlyscg.rds'))
  
}

scg_genes <- annot %>% filter(annotation %in% scg) %>% pull(genename)

save_abtbl_rds(ab.rawcount.tbl, 
               keep = ab.rawcount.tbl$genename %>% unique(),
               scg_genes = scg_genes, name = 'raw_all')

save_abtbl_rds(ab.rawcount.tbl,
               keep = vars.keep,
               scg_genes = scg_genes, name = 'raw_filt0rem_8sams')

save_abtbl_rds(ab.alr.tbl %>% 
                 pivot_longer(names_to = 'sample',
                              values_to = 'count', cols = -c(genename)) ,
               keep = vars.keep,
               scg_genes = scg_genes, name = 'alr_filt0rem_8sams')

save_abtbl_rds(ab.alr.tbl %>% 
                 pivot_longer(names_to = 'sample',
                              values_to = 'count', cols = -c(genename)) %>% 
                 group_by(sample) %>%
                 filter( count > min(count)),
               keep = unique(ab.rawcount.tbl$genename),
               scg_genes = scg_genes, name = 'alr_all')

colnames(ab.alr.keggagg) <- str_replace_all(colnames(ab.alr.keggagg),
                                            c( "BL100413" = "BL100414",
                                               "BL110704" = "BL110705",
                                               "BL120518" = "BL120511",
                                               "BL131204" = "BL131215"))

write_rds(ab.alr.keggagg,
          'data/04_table_gen/alragg.tbl.rds')


colnames(ab.alr.keggagg.nolog) <- str_replace_all(colnames(ab.alr.keggagg.nolog),
                                            c( "BL100413" = "BL100414",
                                               "BL110704" = "BL110705",
                                               "BL120518" = "BL120511",
                                               "BL131204" = "BL131215"))

write_rds(ab.alr.keggagg.nolog ,
          'data/04_table_gen/alragg.nolog.tbl.rds')

# Taxonomy formatting -----------------------------------------------------
tax <- read_tsv('data/04_taxgenes/taxonomyresult.tsv',
                col_names = c('genename', 'genomeid', 'last_rank',
                              'rank_name', 'taxonomy')) %>% 
  separate(col = taxonomy,
           into = c("species","genus","family",
                    "order","class","phylum","superkingdom"), sep = ";") %>% 
  mutate( label = case_when(is.na(phylum) ~ rank_name,
                            !is.na(phylum) ~ str_c(phylum, family,sep = ", "),
                            TRUE ~ rank_name)) %>% 
  distinct(genename, .keep_all = T) 

saveRDS(tax, 'data/04_table_gen/taxonomy_formatted.rds')

annot_noscg <-  annot %>% 
  filter(!annotation %in% scg) %>% 
  pull(genename)

tax %>% 
  filter(genename %in% annot_noscg) %>% 
  saveRDS(., 'data/04_table_gen/taxonomy_subset.rds')
