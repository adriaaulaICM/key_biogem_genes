library(tidyverse)
library(patchwork)
library(umap)
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
  
tax <- readRDS('data/04_table_gen/taxonomy_merged_contig_sgene_gtdb.rds')  %>% 
  filter(genename %in% ab.raw.filt$genename) 

sel.season <- read_rds('data/intermediate_files/season_maxima.rds')

# Param -------------------------------------------------------------------
# theme ordination 
theme_ordination <- list( lil.strip,
                          xlab(NULL), 
                          ylab(NULL),
                          theme(axis.ticks = element_blank(),
                                axis.text = element_blank(),
                                strip.text.x = element_markdown(),
                                legend.text = element_markdown()))


# Selection and subsetting ------------------------------------------------

genes.wtax <- tax %>% 
  filter(genename %in% lomb$genename) %>% 
  filter(family %in% family.order) %>%
  pull(genename)



# Plotting by families  ----------------------------------------------------------------

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

# Main plot -----------------------------------------------------------

season.max <- ptol_pal()(8)
names(season.max) <- season.trans

# Example -----------------------------------------------------------------

sel <- c('BL131204.k127_3527521_1.pro_1927797',
         'BL131204.k127_770738_1.pro_506787',
         'BL101116.k127_1091631.mgm_1234259',
         'BL140505.k127_2488324.mgm_2715285',
         'BL090113.k127_1040684.mgm_1210741'
         )

example.dataset <-  ab.raw %>% 
  filter(genename %in% sel) %>% 
  left_join(sel.season, by = 'genename') %>% 
  left_join(envdata, by = c('sample' = 'Sample_name')) %>% 
  mutate(season.x = fct_explicit_na(season.x, na_level = 'Non-seasonal'))

examplot <- example.dataset %>% 
  ggplot(aes(day_of_year, 2^count)) + 
  geom_point(aes(color = season.x), show.legend = F) + 
  stat_smooth(aes(x = day_of_year,
                  group = genename,
                  color = season.x),
              method = "gam",
              size =3,
              se = F,
              formula = y ~ s(x, k =12, bs = 'cc'), 
              show.legend = F) +
  facet_wrap(~season.x, ncol = 1, scales = 'free_y') + 
  scale_color_manual(values = c(season.max, `Non-seasonal` = 'black')) + 
  theme_void()

examplot

label.pos <- data.frame( season = factor(season.trans, levels = season.trans),
                         V1 = c(2.5, 2.3, 2, 0, -2, 2, 2.5, -1),
                         V2 = c(4, 2.2, 0.5, -2.5, -5, -2, 6, 6))


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
                  color = season), 
              inherit.aes = F,
              show.legend = F) +
  theme_void() + 
  theme_ordination +
  leg.bottom + 
  scale_fill_ptol() +
  scale_color_ptol()

colors.by.fam <- coord.dat %>% 
  filter(!is.na(season)) %>% 
  ggplot(aes( V1, V2)) + 
  geom_point( aes(color = season), size = 1, show.legend = F) + 
  facet_wrap(~family, scales = 'free') + 
  scale_color_ptol()  + 
  theme_ordination +
  theme(legend.position = c(0.85, 0.17))

# Counts per family per season --------------------------------------------

plot_gene_season <- function(df){
  df %>% 
    ggplot(aes(count, fct_rev(season), fill = gene_name)) + 
    geom_bar(stat = 'identity', show.legend = T) + 
    facet_wrap(~family, scales = 'free_x') + 
    ylab('Season') + 
    xlab('Count (n. variants)') + 
    lil.strip + 
    scale_fill_manual(values = palette.it, name = 'Gene name') + 
    theme(strip.text.x = element_markdown(),
          legend.text = element_markdown(),
          legend.position = c(0.75, 0.15),
          legend.direction = 'horizontal') +
    guides(fill = guide_legend(title.position = 'top',
                               ncol = 3))
  
}

coord.fam <- coord.dat %>% 
  filter(!is.na(season)) %>% 
  group_by(family,season , cycle, gene_name) %>% 
  summarize(count = n())

phosphorous.plot <- coord.fam  %>% 
  filter(cycle == 'phosporous') %>% 
  plot_gene_season()
 

other.genes.plot <- coord.fam  %>% 
  filter(cycle != 'phosporous') %>% 
  plot_gene_season()

composite <- distribution + 
  colors.by.fam + 
  phosphorous.plot + 
  other.genes.plot + 
  plot_layout(ncol = 2) + 
  plot_annotation(tag_levels = 'A' )

ggsave(filename = str_c(fig.path, '/',
                        'ordination_composite_old.pdf'),
      plot = composite,
      width = 15,
      height = 15)
  


# Composition alternative  ------------------------------------------------

# areas <- c(area(10, 10, 50, 10),
#            area(10, 20, 60, 60))
# 
# # plot(areas)
# 
# exam.final <- examplot + 
#   theme_classic() + 
#   lil.strip + 
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank()) + 
#   ylab('Abundance ratio') + 
#   xlab('Day of year')
# 
# exam.final + 
#   distribution + 
#   plot_layout(design = areas)
# 
# blueline <- ggplot() + 
#   geom_abline(slope = -1,
#               intercept = 1,
#               color = 'blue',
#               linetype = 2) + 
#   theme_void()

# Interactive version -----------------------------------------------------

# Making an interactive version to check some properties of our data. 
# library(plotly)
# coord.plot <-  coord.dat %>% 
#   ggplot(aes( V1, V2)) + 
#   geom_point( aes(color = gene_name,
#                   shape = str_to_title(season),
#                   text = paste(genename,
#                                "\ngenus: ", genus,
#                                "\nspecies: ", species))) +
#   facet_wrap(~family) + 
#   scale_color_manual(values = palette.it, name = 'gene') + 
#   scale_shape_manual(values = shape.val, name = 'season maxima') + 
#   leg.bottom
# 
# 
# pint <- ggplotly(coord.plot, tooltip = 'text')
# 
# library(htmlwidgets)
# saveWidget(pint, "results/figures/ordination_variants/interactive_ordination.html", 
#            selfcontained = TRUE)


# Plotting each gene ------------------------------------------------------

splitted.plots <- coord.dat %>% 
  filter(!gene_name %in% c('PR blue', 'PR green (L105)',
                          'PR green (M105)', 'Other PR')) %>% 
  mutate(gene_name = as.character(gene_name)) %>% 
  split(.$gene_name) %>% 
  map(~ggplot(., aes( V1, V2)) + 
        geom_point( aes(color = gene_name, shape = str_to_title(season))) +
        facet_wrap(~family) + 
        scale_color_manual(values = palette.it, name = 'Gene') + 
        scale_shape_manual(values = shape.val, name = 'Season maxima') + 
        lil.strip + 
        theme_ordination + 
        leg.bottom
  ) 

dir.create(str_c(fig.path, '/splitted_genes/'))

filenames <- str_c(fig.path,
                   '/splitted_genes',
                   '/umap_ord_',
                   names(splitted.plots) %>% str_replace_all('\\*', ''),
                   '.pdf')

map2(filenames, 
     splitted.plots, 
     ggsave,
     width = 9, 
     height = 9)


# Specific distributions --------------------------------------------------

pr.plot <-  coord.dat %>% 
  filter(gene_name %in% c('PR blue', 'PR green (L105)',
                          'PR green (M105)', 'Other PR')) %>% 
  ggplot(aes( V1, V2)) + 
  geom_point( aes(color = gene_name, shape = str_to_title(season))) +
  facet_wrap(~family) + 
  scale_color_manual(values = palette.it, name = 'Gene') + 
  scale_shape_manual(values = shape.val, name = 'Season maxima') + 
  theme_ordination +
  leg.bottom

ggsave(plot = pr.plot,
       str_c(fig.path, '/umap_ord_proteorhodopsins.pdf'),
       width = 11, height = 9)

# Let's have more detail for this one 
coord.dat %>% 
  filter(gene_name %in% c('PR blue', 'PR green (L105)',
                          'PR green (M105)', 'Other PR')) %>% 
  filter(family == '*Puniceispirillaceae*') %>%
  # left_join(tax %>% select(genename, genus), by = 'genename') %>% 
  ggplot(aes( V1, V2)) + 
  ggrepel::geom_label_repel( aes(label = genus)) + 
  geom_point( aes(color = gene_name, shape = str_to_title(season))) +
  facet_wrap(~family) + 
  scale_color_manual(values = palette.it, name = 'Gene') + 
  scale_shape_manual(values = shape.val, name = 'Season maxima') + 
  theme_ordination +
  leg.bottom

ggsave('results/figures/ordination_variants/umap_proteo_puni.pdf')
 

# and what about fecA and Flavobacteriaceae
flavo.fecA <- coord.dat %>% 
  filter(gene_name %in% c('*fecA*')) %>% 
  filter(family == '*Flavobacteriaceae*') %>%
  # left_join(tax %>% select(genename, genus,species), by = 'genename') %>% 
  ggplot(aes( V1, V2)) + 
  # ggrepel::geom_label_repel( aes(color = genus)) + 
  geom_point( aes(color = gene_name, shape = str_to_title(season))) +
  facet_wrap(~genus) + 
  scale_color_manual(values = palette.it, name = 'Gene') + 
  scale_shape_discrete(name = 'Season maxima') + 
  theme_ordination + 
  leg.bottom

ggsave('results/figures/ordination_variants/umap_proteo_flavo_genus.pdf', 
       plot = flavo.fecA, 
       width = 8,
       height = 8)



# Plotting the whole dataset ----------------------------------------------
# coord <- read_rds('data/intermediate_files/umap_coord.rds')
# 
# coord.2d.plot <- coord %>% 
#   left_join(sel.season, by = 'genename') %>% 
#   filter(!is.na(season)) %>% 
#   ggplot(aes( V1, V2)) + 
#   stat_density_2d(geom = "polygon",
#                   aes(alpha = ..level.., fill = season)) + 
#   facet_wrap(~season)
# 
# ggsave(plot = coord.2d.plot,
#        str_c(fig.path, '/umap_ord_2d.pdf'),
#        width = 11, height = 9)

# Plotting by seasons -----------------------------------------------------


# distribution <- coord.dat %>% 
#   filter(!is.na(season)) %>% 
#   ggplot(aes( V1, V2)) + 
#   geom_point( aes(color = season), size = 3, alpha = 0.7  )  + 
#   facet_wrap( ~ season) + 
#   scale_color_ptol()
# 
# 
# 
# ggsave(plot = distribution, 
#        str_c(fig.path, '/umap_ord_seasonal_distribution.pdf'),
#        width = 9, height = 9)


# coord.plot <-  coord.dat %>% 
#   ggplot(aes( V1, V2)) + 
#   geom_point( aes(color = gene_name, shape = str_to_title(season))) +
#   facet_wrap(~family) + 
#   scale_color_manual(values = palette.it, name = 'Gene') + 
#   scale_shape_manual(values = shape.val, name = 'Season maxima') + 
#   leg.bottom
# 
# ggsave(plot = coord.plot,
#        str_c(fig.path, '/umap_ord_seasonal_genes.pdf'),
#        width = 13, height = 9)

# same, but for genes --------------------------------------------

# gene.count.plot <- coord.dat %>% 
#   filter(!is.na(season)) %>% 
#   group_by(gene_name, season) %>% 
#   summarize(count = n()) %>% 
#   mutate(season = fct_rev(season)) %>% 
#   ggplot(aes(count, season, fill = season)) + 
#   geom_bar(stat = 'identity', show.legend = F) + 
#   facet_wrap(~gene_name) + 
#   scale_fill_ptol()  + 
#   ylab('Season') + 
#   xlab('Count (n. variants)') + 
#   lil.strip + 
#   theme(strip.text.x = element_markdown())
# 
# gene.count.plot
# 
# ggsave(filename = str_c(fig.path, '/',
#                         'ordination_counts.pdf'),
#        plot =  gene.count.plot, 
#        width = 8,
#        height = 10)
