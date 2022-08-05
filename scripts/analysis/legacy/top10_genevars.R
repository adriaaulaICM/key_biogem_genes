library(tidyverse)
source('scripts/analysis/sourcefiles/backbone_params-graphs.R')

theme_set(theme_bw(base_size = 20))

scg <- c("K06942", "K01889", "K01887", "K01875", "K01883",
         "K01869", "K01873", "K01409", "K03106", "K03110")

ab.tbl <- read_rds('data/04_table_gen/scglengthnorm.tbl.rds') %>% 
  filter(!keggid %in% scg)  

colnames(ab.tbl) <-  str_replace_all(colnames(ab.tbl),
                                     c( "BL100413" = "BL100414",
                                        "BL110704" = "BL110705",
                                        "BL120518" = "BL120511",
                                        "BL131204" = "BL131215"))

ab.raw <- read_rds('data/04_table_gen/lengthnorm.tbl.rds') 

descriptions <- readxl::read_xlsx('results/summary_charact/database_genes.xlsx')

envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')

megafile.agg <-  ab.raw  %>% 
  gather(key = 'sample', value = 'count', -keggid, -genename)

onlytops <- megafile.agg %>% 
  group_by(genename) %>% 
  filter(sum(count) > 1000)

sel.vars <- onlytops %>% 
  group_by(keggid, genename) %>% 
  summarise( mean.sam = mean(count)) %>% 
  group_by( keggid) %>% 
  top_n(n = 5, wt = mean.sam) %>% 
  pull(genename)


# Drawing distributions  --------------------------------------------------

selvars.df <-  ab.tbl %>% 
  gather(key = 'sample', value = 'logratio', -keggid, -genename) %>% 
  filter(genename %in% sel.vars) %>% 
  left_join(descriptions, by = 'keggid') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name'))

plots.cool <- unique(selvars.df$gene_name) %>% 
  map(~ggplot(selvars.df %>% filter(gene_name == .x)%>%
                mutate(genename = factor(genename,
                                         levels = unique(genename),
                                         labels = str_c('Var. ',
                                                        1:length(unique(genename))))),
              aes( month, logratio)) + 
        geom_jitter( aes(fill = genename),
                     shape = 21, alpha = 0.6, size = 3.2) + 
        stat_smooth(aes(x = month,
                        color = genename,
                        group = genename), # continuous x-axis
                    method = "lm",
                    formula = y ~ poly(x, 3))  + 
        guides(fill = guide_legend(title = 'nasA variant',
                                   override.aes= list(size =7, alpha  =1)),
               color = FALSE) + 
        scale_color_stata() + 
        scale_fill_stata() + 
        ggtitle(.x) + 
        ylab('Log10 ratio (gene / single copy genes)') + 
        lil.strip + 
        scale_x_month + 
        leg.bottom  +
        theme(axis.text.y = element_text(size =19),
              plot.title = element_text(size = 25, face = 'italic'),
              plot.background = element_blank() ))

plots.cool[[2]]
ggsave(filename = 'results/figures/variant_calling_nasA.pdf', 
       plot = plots.cool[[2]], width =9.5 , height =7.5)


library(propr)


raw.sel <- ab.raw[ab.raw$genename %in% selvars.df$genename, -c(2)]  %>% 
  column_to_rownames( var = 'genename') %>% 
  t()
  
                  
respropr <- propr(counts = raw.sel, metric = 'rho', updateCu)

correspondence <- selvars.df %>% 
  select(genename, keggid) %>% 
  distinct()


coolname <- selvars.df %>%
  filter(gene_name == 'nasA') %>%
  mutate(varnum = factor(genename,
                           levels = unique(genename),
                           labels = str_c('Var. ',
                                          1:length(unique(genename))))) %>% 
  select(genename, varnum) %>% 
  distinct()

prop <- respropr@matrix %>%  
  as.data.frame() %>% 
  rownames_to_column(var = 'Partner') %>% 
  gather(key = 'Pair', value = 'rho', -Partner) %>% 
  left_join(correspondence, by = c('Partner' = 'genename')) %>% 
  left_join(correspondence, by = c('Pair' = 'genename')) %>% 
  left_join(descriptions, by = c('keggid.x' = 'keggid')) %>% 
  left_join(descriptions, by = c('keggid.y' = 'keggid')) %>% 
  left_join(coolname, by = c('Partner' = 'genename')) %>% 
  left_join(coolname, by = c('Pair' = 'genename')) %>% 
  filter( keggid.x == keggid.y) 


prop %>% 
  filter( gene_name.x == 'nasA')  %>%
  ggplot( aes(varnum.x, varnum.y, fill = rho)) + 
  geom_tile() + 
  ylab(NULL) + 
  xlab(NULL) + 
  theme_minimal(base_size = 21) + 
  theme( axis.text.y = element_text(size = 21)) + 
  scale_fill_viridis_c()

ggsave('results/figures/heatmap_nasA.pdf', width = 9, height = 7)

hist.rho <- prop %>% 
  filter( keggid.x == keggid.y) %>% 
  filter(rho != 1) %>% 
  select(rho) %>% 
  ggplot( aes(rho)) + 
  geom_histogram() + 
  ggtitle('Distribution of rho value within same functional group',
          subtitle = 'High absolute values indicate covariance') + 
  ylab( 'Count' )  + 
  
  xlab( 'Rho value')

hist.rho
ggsave('results/figures/rho_histogram.pdf', height = 7, width = 9)


prop %>%
  filter(!is.na(Pair)) %>% 
  ggplot(aes(Partner, Pair)) + 
  geom_tile(aes(fill = propr))  + 
  scale_fill_viridis_c()
