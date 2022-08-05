library(tidyverse)
source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')

lomb <- readRDS('data/intermediate_files/genes_seasonal.rds') 

ab.raw.filt <- read_rds('data/04_table_gen/abtblalr_filt0rem_8sams.rds')  %>% 
  filter(genename %in% lomb$genename)
  
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')

sel.season <- ab.raw.filt %>% 
  left_join(envdata, by = c('sample' = 'Sample_name'))  %>% 
  # for each gene, group by year
  group_by(genename, year) %>% 
  # we take for each year the month with the max count
  filter(count == max(count)) %>% 
  # we group again, this time by season
  group_by(genename, season) %>% 
  summarize( ntimes.max = n()) %>% 
  # we only keep strong maximums, repeating more than one year
  filter(ntimes.max >= 3) %>% 
  select(-ntimes.max)


pair.dataset <- sel.season %>%
  group_by(genename) %>%
  filter(n() >= 2) %>%
  arrange(genename, season) %>% 
  split(.$genename)  %>% 
  map( ~(.x$season) ) %>% 
  map_chr(~str_c(str_to_title(.x[1]),
                 str_to_title(.x[2]),
                 sep = '-')) %>% 
  as_tibble(rownames = 'genename') %>% 
  mutate(season = case_when(
    value == 'Autumn-Spring' ~ 'unclear',
    value == 'Autumn-Summer' ~ 'Summer-Autumn',
    value == 'Spring-Winter' ~ 'Winter-Spring',
    TRUE ~ value )) %>% 
  select(-value)


sel.season <- sel.season %>% 
  filter(!genename %in% pair.dataset$genename) %>% 
  bind_rows(pair.dataset) %>% 
  mutate(season = ifelse(season %in% season.order,
                         str_to_title(season),
                         season)) %>% 
  mutate(season = factor(season, levels = c('Winter', 'Winter-Spring',
                                            'Spring', 'Spring-Summer',
                                            'Summer', 'Summer-Autumn',
                                            'Autumn', 'Autumn-Winter'))) %>% 
  ungroup()


write_rds(x = sel.season, 
          path =  'data/intermediate_files/season_maxima.rds')


# we will do the same with the months in a simpler fashion 

sel.months <- ab.raw.filt %>% 
  left_join(envdata, by = c('sample' = 'Sample_name'))  %>% 
  # for each gene, group by year
  group_by(genename, year) %>% 
  # we take for each year the month with the max count
  filter(count == max(count)) %>% 
  # we group again, this time by season
  group_by(genename, month) %>% 
  summarize( ntimes.max = n()) %>% 
  # we only keep strong maximums, repeating more than one year
  filter(ntimes.max >= 3)  %>% 
  select(-ntimes.max)


write_rds(x = sel.months, 
          path =  'data/intermediate_files/month_maxima.rds')
