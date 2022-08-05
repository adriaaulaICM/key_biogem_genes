library(tidyverse)
library(phyloseq)
library(corncob)

scg <- c("K06942", "K01889", "K01887", "K01875", "K01883",
         "K01869", "K01873", "K01409", "K03106", "K03110")

ab.raw <- read_rds('data/04_table_gen/lengthnorm.tbl.rds') %>% 
  filter(!keggid %in% scg) 

colnames(ab.raw) <-  str_replace_all(colnames(ab.raw),
                                     c( "BL100413" = "BL100414",
                                        "BL110704" = "BL110705",
                                        "BL120518" = "BL120511",
                                        "BL131204" = "BL131215"))
descriptions <- readxl::read_xlsx('results/summary_charact/database_genes.xlsx')

envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv') %>% 
  filter(Sample_name %in% colnames(ab.raw)) %>% 
  distinct(Sample_name, .keep_all = T) %>% 
  column_to_rownames(var = 'Sample_name') 

otuliketab <- ab.raw %>% 
  gather(key = 'sample', value = 'count', -keggid, -genename) %>% 
  group_by(sample, keggid) %>% 
  summarize( total  = sum(count)) %>% 
  spread(key = 'sample', value = 'total') %>% 
  column_to_rownames(var = 'keggid') %>% 
  ceiling() %>% 
  as.data.frame()

phyobj <- phyloseq(otu_table(as.matrix(otuliketab), taxa_are_rows = T),
                   sample_data(envdata))

da_analysis <- differentialTest(formula = ~ season,
                                phi.formula = ~ season,
                                formula_null = ~ 1,
                                phi.formula_null = ~ season,
                                test = "Wald", boot = FALSE,
                                data = phyobj,
                                fdr_cutoff = 0.01)

plot(da_analysis)
taxa_names(phyobj)
descriptions[descriptions$keggid %in% da_analysis$significant_taxa,]

