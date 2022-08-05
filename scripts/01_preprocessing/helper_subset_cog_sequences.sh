library(tidyverse)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)

pathdna<- args[1] #where the fastq files are located
pathaa <-args[2] #name of analysis directory.
cogdb <-args[3] #name of analysis directory.
output.name <-args[4] #name of analysis directory.

coggi <- read_tsv( cogdb, 
                  col_names = c('genename', 'cog', 'gene', 'cosa', 'start',
                                'end', 'x2', 'x3', 'bitscore',
                                'eval', 'function'))

dna <- readDNAStringSet(pathdna)
aa <- readAAStringSet(pathaa)

names(dna) <- str_trim(names(dna))
names(aa) <- str_trim(names(aa))

dna.sel <- dna[coggi$genename]
aa.sel <- aa[coggi$genename]

dna.subseq <- subseq(dna.sel,
                     start = (coggi$start * 3) - 2,
                     end = (coggi$end)*3)
aa.subseq <- subseq(aa.sel, start = coggi$start, end = coggi$end)

whole.dna <- c(dna.subseq, dna[!(names(dna) %in% coggi$genename)])
whole.aa <- c(aa.subseq, aa[!(names(aa) %in% coggi$genename)])

writeXStringSet(whole.dna, str_c(output.name, '.fna'))
writeXStringSet(whole.aa, str_c(output.name, '.faa'))




