library(tidyverse)
library(cowplot)
library(patchwork)
library(ggtext)
library(ggforce)
source('scripts/utils/backbone_params-graphs.R')
source('scripts/utils/backbone_functions.R')
source('scripts/utils/params_biogem.R')

set.seed(42)


# Data import -------------------------------------------------------------
fig.path <- 'results/figures/ordination_variants'
dir.create(fig.path)


lomb <- readRDS('data/intermediate_files/genes_seasonal.rds') 

ab.raw <- read_rds('data/04_table_gen/abtblalr_filt0rem_8sams.rds')  
ab.raw.filt <- ab.raw  %>% 
  filter(genename %in% lomb$genename)
  
envdata <- read_tsv('data/metadata-raw/metadata_Blanes_compact_may2017.tsv')
annot <- read_rds('data/intermediate_files/annotation.rds')
  
tax <- readRDS('data/04_table_gen/taxonomy_contigwise.rds')  %>% 
  filter(genename %in% ab.raw.filt$genename) 

sel.season <- read_rds('data/intermediate_files/season_maxima.rds')

# Parameters -------------------------------------------------------------------


# Selection and subsetting ------------------------------------------------

genes.wtax <- tax %>% 
  filter(genename %in% lomb$genename) %>% 
  filter(family %in% family.order) %>%
  pull(genename)


# Plotting by families  --------------------------------------------------------

palette.it <- palette.all[gene.order]
names(palette.it) <- gene.italics

coord.fam <- read_rds('data/intermediate_files/umap_coord_fam.rds')

coord.dat <- coord.fam %>% 
  left_join(tax %>% select(genename, species, genus, family), by = 'genename') %>% 
  mutate(family = factor(family,
                         levels = family.order,
                         labels = family.labels)) %>% 
  mutate( gene_name = factor(gene_name,
                             levels = gene.order,
                             labels = gene.italics)) %>% 
  left_join(sel.season, by = 'genename') 



# Example -----------------------------------------------------------------

season.max <- ptol_pal()(8)
names(season.max) <- season.trans

sel <- c('BL131204.k127_3527521_1.pro_1927797',
         # 'BL131204.k127_770738_1.pro_506787', autumn
         # 'BL101116.k127_1091631.mgm_1234259', possibly summer? 
         'BL140505.k127_2488324.mgm_2715285',
         'BL090113.k127_1040684.mgm_1210741',
         'BL130115.k127_2918295.mgm_3479586'
         )

example.dataset <-  ab.raw %>% 
  filter(genename %in% sel) %>% 
  left_join(sel.season, by = 'genename') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) %>% 
  mutate(season.x = fct_explicit_na(season.x, na_level = 'Non-seasonal'))

plot_examples <- function(df){
  
  df %>% 
    ggplot(aes(day_of_year, 2^count)) + 
    geom_point( color = 'grey60',
                size = 0.8,
               alpha = 0.8, show.legend = F) + 
    stat_smooth(aes(x = day_of_year,
                    group = genename,
                    color = season.x),
                method = "gam",
                size =1,
                se = F,
                formula = y ~ s(x, k =12, bs = 'cc'), 
                show.legend = F) +
    scale_color_manual(values = c(season.max, `Non-seasonal` = 'black')) + 
    ylab('Gene variant ratio') + 
    xlab('Time') + 
    theme_bw(base_size = 8) + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.background = element_blank(),
          title = element_text(size = 8))
  
}

winter <- example.dataset %>% 
  filter(genename %in% 'BL131204.k127_3527521_1.pro_1927797') %>% 
  plot_examples() + 
  theme(axis.title.x = element_text(color = 'transparent'))
  
spring <- example.dataset %>% 
  filter(genename %in% 'BL140505.k127_2488324.mgm_2715285') %>% 
  plot_examples() 
  
nonseasonal <- example.dataset %>% 
  filter(genename %in% 'BL130115.k127_2918295.mgm_3479586') %>% 
  plot_examples() + 
  ggtitle('Non-seasonal')  



# Colors by family --------------------------------------------------------

colors.by.fam <- coord.dat %>% 
  filter(!is.na(season),
         !is.na(family),
         !family %in% c('HIMB59')) %>% 
  ggplot(aes( V1, V2)) + 
  geom_point( aes(color = season), size = 1, show.legend = F) + 
  facet_wrap(~family, scales = 'fixed') + 
  scale_color_ptol()  + 
  theme(strip.background = element_blank(),
        strip.text.x =element_text(size = 8,
                                   margin = margin(0, 0.2, 0, 0, "cm"))) + 
  xlab(NULL) + 
  ylab(NULL) + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        strip.text.x = element_markdown(),
        legend.text = element_markdown()) 


# Main plot ---------------------------------------------------------------


label.pos <- data.frame( season = factor(season.trans, levels = season.trans),
                         V1 = c(2.5, 2.3, 2, 0, 2, -1.5, 3, -1),
                         V2 = c(4, 2.2, 0.5, -2.3, -2.5, -5.5, 7, 6))


distribution <- coord.dat %>% 
  filter(!is.na(season)) %>% 
  ggplot(aes( V1, V2, fill = season)) + 
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level..),
                  show.legend = F) +
  geom_point( size = 3, shape = 21, alpha = 0.7,
              show.legend = F)  + 
  geom_label( data = label.pos,
              aes(V1,V2,
                  label = season,
                  fill = season),
              color = 'white',
              inherit.aes = F,
              show.legend = F) +
  coord_cartesian(xlim = c(-4.3, 4), ylim = c(-6, 7)) + 
  theme_void() +
  leg.bottom + 
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(-7, 9)) +
  scale_fill_ptol() +
  scale_color_ptol()


# Barplot number of variants seasonal -------------------------------------

patt.whole <- lomb %>% 
  left_join(sel.season, by = 'genename') %>% 
  filter(!is.na(season)) %>% 
  group_by(season) %>% 
  summarise(total = n()) %>% 
  ggplot(aes(y = fct_relevel(season, season.trans[c(7,8,1,2,3,4,5,6)]) %>%
               fct_rev(),
             x = total)) + 
  geom_bar(stat = 'identity', aes(fill = season), show.legend = F) +
  geom_text(aes(x = 20,
                label = total), 
            colour = 'white',
            hjust = 0,
            size = 3.5) + 
  coord_cartesian(xlim = c(0, 1500)) + 
  lil.strip +  
  ylab(NULL) + 
  xlab('Count (n. variants)') + 
  ggtitle('Number of seasonal gene variants') + 
  scale_fill_ptol() + 
  scale_x_continuous(expand = c(0,10)) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


# Composition -------------------------------------------------------------


arrowA <- data.frame(x1 = 15, x2 = 49, y1 = 90, y2 = 71)
arrowB <- data.frame(x1 = 15, x2 = 52, y1 = 75, y2 = 55)

coolthing <- ggdraw(xlim = c(0, 100), ylim = c(0, 100)) + 
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = arrowA, 
               arrow = arrow(type = 'closed', length = unit(0.2, "cm") ),
               linetype = 3,
               color = 'grey60',
               size  = 0.4,
               lineend = "round") +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = arrowB, 
               arrow = arrow(type = 'closed', length = unit(0.2, "cm") ),
               linetype = 3,
               color = 'grey60',
               size  = 0.4,
               lineend = "round") + 
  draw_plot(distribution, x = 5, y = 0, width = 95, height = 100)  + 
  draw_plot(winter, x = 0, y = 85, width = 16, height = 16) + 
  draw_plot(nonseasonal, x = 5, y = 3, width = 16, height = 16) + 
  draw_plot(spring, x = 0, y = 70, width = 16, height = 16) 

whole.compo <- (coolthing +
    (patt.whole / colors.by.fam)) + 
  plot_layout( widths = c(0.7, 0.3)) + 
  plot_annotation(tag_levels = 'A')


ggsave(filename = 'results/figures/ordination_variants/ordination_composition.pdf',
       plot = whole.compo,
       width = 12, 
       height = 9)

  