
gene.order <- c("PR blue", 
                'PR green (M105)',
                'PR green (L105)',
                'Other PR',
                "pufM", "pufL","rbcL I", 'psbA',
                'coxL',
                "tauA", "chiA",  'narB', 'nasA',
                'nifH', 'hao','amoA', 'ureC', 'phnD', 'phnM', 'pstS', 'phoD',
                'phoX', 'ppx', 'ppk1', 'plcP', 'dmdA', "fecA")


gene.order.main <- c("PR", 
                     "pufM", "pufL","rbcL I", 'psbA',
                     'coxL',
                     "tauA", "chiA",  'narB', 'nasA',
                     'nifH', 'hao','amoA', 'ureC', 'phnD', 'phnM', 'pstS', 'phoD',
                     'phoX', 'ppx', 'ppk1', 'plcP', 'dmdA', "fecA")

gene.italics <- c("PR blue", 
                'PR green (M105)',
                'PR green (L105)',
                'Other PR',
                "*pufM*", "*pufL*","*rbcL I*",'*psbA*',
                '*coxL*',
                "*tauA*", "*chiA*", '*narB*', '*nasA*',
                '*nifH*', '*hao*','*amoA*', '*ureC*', '*phnD*', '*phnM*',
                '*pstS*', '*phoD*', '*phoX*', '*ppx*', '*ppk1*',
                '*plcP*', '*dmdA*', "*fecA*")

gene.italics.sel <- c("amoA", "dddP", "dmdA", "fecA", "narB", "nasA",
                  "phnD", "phnM", "phoD", "phoX", "plcP", "ppk1", 'coxL',
                  "ppx", "psbA", "pstS", "pufL", "pufM", "rbcL I",
                  "tauA", "ureC")


palette.carbon <- c(`Other PR`= "#49525A",
                    `PR green (M105)` =  "#307473",
                    coxL = "#726C5B",
                    pufM = "#B79047", 
                    pufL = "#F3D3BD", 
                    tauA = "#6247AA", 
                    chiA = '#062726',
                    `PR green (L105)` = "#12664F",
                    `PR blue` = "#2DC2BD",
                    `rbcL I` = '#1E4772', 
                    fecA = '#291720' )


palette.nitro <- c(  narB = "#88AD9D",
                     nasA = "#F4E87A",
                     hao = '#DDD535',
                     nifH = 'darkolivegreen',
                     amoA = "#E9B154",
                     ureC = "#E0CD4E")

palette.phos <-   c( phnD = "#C7BDC0",
                     phnM = "#F2B38C",
                     psbA = "#B73F1B",
                     pstS = "#F0774D",
                     phoD = "#DE514E",
                     phoX = "#E56399",
                     ppx = '#C00A0F',
                     ppk1 = "#C37E78",
                     plcP = "#C2556E")
  
palette.sul <- c(dmdA = "#4062BB", 
                 dddP = "#7FD1B9")

palette.all <- c(palette.carbon, palette.nitro, palette.phos, palette.sul)


palette.it <- palette.all[gene.order]
names(palette.it) <- gene.italics

normal.genes <- c("amoA", "coxL", "dmdA", "hao", "phoD",
                  "plcP", "ppk1", "ppx", "pstS", "pufM",
                  "rbcL I", "ureC")

double.genes <- data.frame( genes = c('phnD', 'phnM', 'narB', 'nasA'), 
                            gene_name = c('phnD + phnM',
                                          'phnD + phnM',
                                          'narB + nasA',
                                          'narB + nasA'))

scg <- c("K06942", "K01889", "K01887", "K01875", "K01883",
         "K01869", "K01873", "K01409", "K03106", "K03110")

family.selection <- c('Nitrosopumilaceae',
                      'Rhodobacteraceae',
                      'Puniceispirillaceae',
                      'Pelagibacteraceae',
                      'Flavobacteriaceae',
                      'HIMB59', 'D2472',
                      'Halieaceae',
                      'Cyanobiaceae', 'Litoricolaceae',
                      'Unclassified', 'Other')

family.order <- c( 'Pelagibacteraceae',
                   'HIMB59',
                   'Rhodobacteraceae',
                   'Puniceispirillaceae',
                   'Other Alpha',
                   'Halieaceae',
                   'Litoricolaceae',
                   # 'HTCC2089',
                   'D2472',
                   'Other Gamma',
                   'Flavobacteriaceae',
                   'Nitrosopumilaceae',
                   'Cyanobiaceae',
                   'Other',
                   'unclassified // no assignation')

family.order.euk <- c( 'Pelagibacteraceae',
                       'HIMB59',
                       'Rhodobacteraceae',
                       'Puniceispirillaceae',
                       'Other Alpha',
                       'Halieaceae',
                       'Litoricolaceae',
                       # 'HTCC2089',
                       'D2472',
                       'Other Gamma',
                       'Flavobacteriaceae',
                       'Nitrosopumilaceae',
                       'Cyanobiaceae',
                       'Eukaryota',
                       'Other',
                       'unclassified // no assignation')

# family.colors <- c("#89023E",
#                    "#AF42AE",
#                    "#CC7178",
#                    "#FFA5AB",
#                    "#D0A3BF", #alpha
#                    "#2C8C99",
#                    "#634B66", #tochange
#                    "#141B41", #tochnage
#                    "#175676", #tochange
#                    "#E8C547",
#                    '#1B2021', #gammnew
#                    "#3BB273",
#                    "#95A3A4",
#                    '#404E5C' ) #uncla
  

family.colors.euk <- c("#89023E",
                       "#AF42AE",
                       "#CC7178",
                       "#FFA5AB",
                       "#D0A3BF", #alpha
                       "#2C8C99",
                       "#634B66", 
                       "#141B41", 
                       "#175676", 
                       "#E8C547",#gammnew
                       '#F5AB00', 
                       "#3BB273",
                       "#f8bd7f",
                       "#95A3A4",#other
                       '#404E5C' ) #uncla

family.labels <- c( '*Pelagibacteraceae*',
                    'HIMB59',
                    '*Rhodobacteraceae*',
                    '*Puniceispirillaceae*',
                   'Other Alpha',
                   '*Halieaceae*',
                   # 'HTCC2089<br>(Pseudomonadales)',
                   '*Litoricolaceae*',
                   'D2472 (SAR86)',
                   'Other Gamma',
                   '*Flavobacteriaceae*',
                   '*Nitrosopumilaceae*',
                   '*Cyanobiaceae*',
                   'Other',
                   'unclassified // no assignation')


family.labels.euk <- c( '*Pelagibacteraceae*',
                   'HIMB59',
                   '*Rhodobacteraceae*',
                   '*Puniceispirillaceae*',
                   'Other Alpha',
                   '*Halieaceae*',
                   '*Litoricolaceae*',
                   # 'HTCC2089<br>(Pseudomonadales)',
                   'D2472 (SAR86)',
                   'Other Gamma',
                   '*Flavobacteriaceae*',
                   '*Nitrosopumilaceae*',
                   '*Cyanobiaceae*',
                   'Eukaryota',
                   'Other', 
                   'unclassified // no assignation')
                   
# names(family.colors) <- family.labels
names(family.colors.euk) <- family.labels.euk
  

other.interesting.families <- c('Thalassoarchaeaceae', #PhnD
                                'AAA536-G10', #Punicieispi
                                'TMED127', #until order
                                'Marinoscillaceae',
                                'UBA10066', #PR blue
                                'Thioglobaceae', #PR green
                                'Litoricolaceae')

# Shapes for the seasons

shape.val <- c(15,0,16,1,17,2,18,5)
season.trans <- c('Winter', 'Winter-Spring',
                  'Spring', 'Spring-Summer',
                  'Summer', 'Summer-Autumn',
                  'Autumn', 'Autumn-Winter')

  names(shape.val) <- str_to_title(season.trans)
