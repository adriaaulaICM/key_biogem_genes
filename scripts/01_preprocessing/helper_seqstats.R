library(tidyverse)
library(skimr)


histodf <- read_tsv('~/icm/tunnel/helper_seqstats.txt',
                    col_names = c('type', 'pos', 'count')) %>% 
  filter( count > 0)


histodf %>% 
  arrange(pos)
  
histodf %>% 
  spread( key = type, value = count ) %>% 
  # group_by(type) %>%  
  skim() 
  
    

ggplot( histodf %>% filter(pos < 2000), aes( pos, count, fill = type)) + 
  geom_bar( stat = 'identity') + 
  #facet_wrap(~type) + 
  # scale_y_log10() + 
  theme_bw() + 
  theme(legend.position = c(0.8, 0.8), legend.background = element_blank()) + 
  ylab('N variants') + 
  xlab('Gene length')


histodf %>% 
  group_by(type) %>%  
  summarise( total.count = sum(count))
