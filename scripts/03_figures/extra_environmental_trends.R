library(tidyverse)
library(mgcv) 
library(broom) 


envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv') 
rich.df <- read_tsv('data/intermediate_files/richness_gene_trends.tsv')

# select 
env.trend <- envdata %>% 
  select(day_of_year, month, day, season,Day_length,
         Temperature, PO4, NH4, NO2, NO3, Chla_total, Synechococcus)  %>% 
  pivot_longer(names_to = 'env', values_to = 'values', cols = -c(day_of_year, month, day, season)) %>% 
# create a nested df 
  group_by(env) %>%
  mutate(scaled.values = scale(values)) %>% 
  nest()  %>%
  mutate( model = map(data,
                      ~gam(values ~ s(day_of_year, k =12, bs = 'cc'),
                           data = .)),
          predict = map(model, ~predict(.))) %>%
  unnest(c(predict), ) %>% 
  print()


# Relationship HNA and LNA ------------------------------------------------

# can we say something about the HNA / LNA distinction? 
envdata %>% 
  select(day_of_year, HNA, LNA, Bacteria_joint) %>% 
  pivot_longer(names_to = 'type_bac', values_to = 'counts', cols = -day_of_year) %>% 
  ggplot( aes(day_of_year, counts)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~type_bac)


envdata %>% 
  ggplot( aes(HNA, LNA)) + 
  geom_point()  + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_smooth()



# Blooms!  ----------------------------------------------------------------
envdata %>%
  filter(year >= 2009) %>% 
  ggplot(aes(Date, Chla_total)) + 
  geom_line() + 
  geom_point(aes(color = season), size = 3) 

# Relationship between psbA, dmdA and blooms 

rich.df %>% 
  filter(gene_name %in% c('psbA', 'dmdA')) %>% 
  ggplot( aes(Chla_total, richness)) + 
  geom_point() + 
  facet_wrap(~gene_name)
