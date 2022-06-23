# Obtaining AnnData files

This Python package supports import of H5AD and 10X data files.

For example, you can visit [https://cellxgene.cziscience.com](https://cellxgene.cziscience.com) and
select from among various choices there and download.

Files used for this example:

```
$ ls -Shlr /mini-corpus/anndata
total 58076592
-rw-r--r--  1 testuser  staff    34K Apr 25 11:13 subset_100_100.h5ad
-rw-r--r--  1 testuser  staff   230K May 11 08:08 pbmc-small.h5ad
-rw-r--r--  1 testuser  staff   4.3M May 10 18:12 10x_pbmc68k_reduced.h5ad
-rw-r--r--@ 1 testuser  staff    27M May 13 22:34 af9d8c03-696c-4997-bde8-8ef00844881b.h5ad
-rw-r--r--@ 1 testuser  staff    30M May 13 23:25 issue-74.h5ad
-rw-r--r--@ 1 testuser  staff    30M May 13 22:45 d4db74ad-a129-4b1a-b9da-1b30db86bbe4-issue-74.h5ad
-rw-r--r--@ 1 testuser  staff    32M May 13 22:17 Puck_200903_10.h5ad
-rw-r--r--@ 1 testuser  staff    36M Apr 25 11:13 local3.h5ad
-rw-r--r--  1 testuser  staff    38M May 11 08:08 pbmc3k_processed.h5ad
-rw-r--r--@ 1 testuser  staff    40M May 13 22:21 brown-adipose-tissue-mouse.h5ad
-rw-r--r--  1 testuser  staff    47M Apr 28 15:15 pbmc3k-krilow.h5ad
-rw-r--r--@ 1 testuser  staff    51M May 13 22:34 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e.h5ad
-rw-r--r--@ 1 testuser  staff    55M May 13 22:33 4056cbab-2a32-4c9e-a55f-c930bc793fb6.h5ad
-rw-r--r--@ 1 testuser  staff    70M Apr 25 22:45 human-kidney-tumors-wilms.h5ad
-rw-r--r--@ 1 testuser  staff    99M May 11 09:57 issue-71.h5ad
-rw-r--r--@ 1 testuser  staff   117M May 13 22:20 adult-mouse-cortical-cell-taxonomy.h5ad
-rw-r--r--@ 1 testuser  staff   119M Apr 25 22:47 longitudinal-profiling-49.h5ad
-rw-r--r--@ 1 testuser  staff   221M Apr 25 22:46 single-cell-transcriptomes.h5ad
-rw-r--r--@ 1 testuser  staff   230M Apr 29 17:13 Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection.h5ad
-rw-r--r--@ 1 testuser  staff   231M Apr 25 11:13 vieira19_Alveoli_and_parenchyma_anonymised.processed.h5ad
-rw-r--r--@ 1 testuser  staff   312M May 11 08:12 not
-rw-r--r--@ 1 testuser  staff   329M Apr 25 11:13 local2.h5ad
-rw-r--r--@ 1 testuser  staff   357M May 10 20:08 issue-69.h5ad
-rw-r--r--@ 1 testuser  staff   376M Apr 25 22:48 acute-covid19-cohort.h5ad
-rw-r--r--@ 1 testuser  staff   686M Apr 25 22:50 autoimmunity-pbmcs.h5ad
-rw-r--r--@ 1 testuser  staff   712M May 13 22:22 developmental-single-cell-atlas-of-the-murine-lung.h5ad
-rw-r--r--@ 1 testuser  staff   2.5G May 13 22:35 tabula-sapiens-stromal.h5ad
-rw-r--r--@ 1 testuser  staff   3.2G May 13 22:30 tabula-sapiens-epithelial.h5ad
-rw-r--r--@ 1 testuser  staff   5.6G Apr 25 23:04 integrated-human-lung-cell-atlas.h5ad
-rw-r--r--@ 1 testuser  staff   5.7G May 13 22:38 tabula-sapiens-immune.h5ad
-rw-r--r--@ 1 testuser  staff   6.6G May 13 22:40 azimuth-meta-analysis.h5ad
```

```
$ du -hs /mini-corpus/anndata
 27G  /mini-corpus/anndata
```
