# key_biogem_genes
Code to evaluate the seasonality of biogeochemically relevant microbial genes in a coastal ocean microbiome. 

In the `data` folder you can find: 

- The FASTA sequence for each gene variant (sequences_min250.faa.zip).
- The abundance count table for each sequence (abundance_tbl_raw.rds).
- Taxonomy performed with GTDB and Uniref (taxonomy_gtdb.tsv | taxonomy_uniref.tsv). 

In the script folder the code is organized as follows: 

The scripts present the following structure:

```
scripts/
    01_preprocessing/
    02_analysis/
    03_figures/
    03_statistics/
    03_tables/
    utils/
```

In `01_preprocessing`, from the raw data I:
- Annotate the genes of interest.
- Create a taxonomy for each gene variant.
- Generate the tables used in posterior analyses.

In `02_analysis`, intermediate files are calculated to facilitate the obtention of statistics and figures.

In `03_figures, 03_statistics` and `03_tables` the results of the paper are obtained.

In `utils` there are scripts called from the `03_figures` script to avoid redundancies.
