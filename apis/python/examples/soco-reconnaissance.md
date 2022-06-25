Next, let's do some cross-cutting queries over schemas of all SOMAs in the collection. The goal is
-- in preparation for a collection-level query -- to find out which `obs` columns, and which values
in those columns, are most likely to be promising in terms of yielding results given our
mini-corpus.

## Total cell-counts

The mini-corpus we prepared is 29 SOMAs, 26GB total:

```
$ du -hs /mini-corpus/tiledb-data
 26G  /mini-corpus/tiledb-data

$ ls /mini-corpus/tiledb-data | wc -l
      29
```

This collection includes data on about 2.4 million cells:

```
import tiledbsc
soco = tiledbsc.SOMACollection('/mini-corpus/soco')

>>> sum(len(soma.obs) for soma in soco)
2464363

>>> [len(soma.obs) for soma in soco]
[264824, 4636, 6288, 2223, 59506, 100, 2638, 982538, 385, 67794, 2638, 104148, 44721, 3799, 11574, 1679, 3589, 700, 584884, 16245, 4603, 3726, 4636, 7348, 3589, 40268, 12971, 4232, 80, 82478, 97499, 38024]
```

```
>>> for soma in soco:
...     print(len(soma.obs), soma.name)
...
264824 tabula-sapiens-immune
4636 wilms-tumors-seurat
6288 issue-69
2223 brown-adipose-tissue-mouse
59506 acute-covid19-cohort
100 subset_100_100
2638 pbmc3k_processed
982538 azimuth-meta-analysis
385 local3
67794 developmental-single-cell-atlas-of-the-murine-lung
2638 pbmc3k-krilow
104148 tabula-sapiens-epithelial
44721 Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection
3799 adipocytes-seurat
11574 longitudinal-profiling-49
1679 adult-mouse-cortical-cell-taxonomy
3589 issue-74
700 10x_pbmc68k_reduced
584884 integrated-human-lung-cell-atlas
16245 issue-71
4603 4056cbab-2a32-4c9e-a55f-c930bc793fb6
3726 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
4636 human-kidney-tumors-wilms
7348 local2
3589 d4db74ad-a129-4b1a-b9da-1b30db86bbe4-issue-74
40268 single-cell-transcriptomes
12971 vieira19_Alveoli_and_parenchyma_anonymised.processed
4232 af9d8c03-696c-4997-bde8-8ef00844881b
80 pbmc-small
82478 tabula-sapiens-stromal
97499 autoimmunity-pbmcs
38024 Puck_200903_10
```

## Cell-counts before running a query

Before running a query, we may wish to know how many cells will be involved in the result:

```
>>> [soma.obs.query('cell_type == "B cell"').size for soma in soco if 'cell_type' in soma.obs.keys()]
[514982, 0, 0, 14283, 245240, 391446, 0, 0, 0, 125154, 0, 29060, 0, 0, 311259, 176, 6120, 2480, 0, 0, 0, 26325, 5220, 0, 12750, 0]

>>> sum([soma.obs.query('cell_type == "B cell"').size for soma in soco if 'cell_type' in soma.obs.keys()])
1684495
```

```
>>> [soma.obs.query('cell_type == "leukocyte"').size for soma in soco if 'cell_type' in soma.obs.keys()]
[59436, 5616, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5472, 0, 0, 0, 0, 0, 0, 3753]

>>> sum([soma.obs.query('cell_type == "leukocyte"').size for soma in soco if 'cell_type' in soma.obs.keys()])
74277
```

## Datasets having all three of obs.cell_type, obs.tissue, and obs.feature_name

```
names = sorted([
  soma.name for soma in soco
    if 'cell_type' in soma.obs.keys() and 'tissue' in soma.obs.keys() and 'feature_name' in soma.var.keys()
])
for name in names: print(name)
```

```
0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
4056cbab-2a32-4c9e-a55f-c930bc793fb6
Puck_200903_10
acute-covid19-cohort
adult-mouse-cortical-cell-taxonomy
af9d8c03-696c-4997-bde8-8ef00844881b
autoimmunity-pbmcs
azimuth-meta-analysis
brown-adipose-tissue-mouse
developmental-single-cell-atlas-of-the-murine-lung
human-kidney-tumors-wilms
integrated-human-lung-cell-atlas
local2
local3
longitudinal-profiling-49
single-cell-transcriptomes
tabula-sapiens-epithelial
tabula-sapiens-immune
tabula-sapiens-stromal
```

## Show counts of obs_ids and var_ids across the collection

Using [./collection-counts.py](collection-counts.py) we can answer questions such as _How many cells will
be involved if I do a query?_ Since these pre-counts operate on the smaller `obs` arrays, they run
faster than going ahead and doing full queries (as shown below) on the larger `X` arrays.

```
----------------------------------------------------------------
Per-SOMA values for cell_type

SOMA acute-covid19-cohort
                                 count
name                                  
monocyte                         29878
CD8-positive, alpha-beta T cell   8658
CD4-positive, alpha-beta T cell   6726
B cell                            6131
natural killer cell               3248
alpha-beta T cell                 1659
dendritic cell                    1038
platelet                          1007
plasmablast                        586
plasmacytoid dendritic cell        575

...
```

```
----------------------------------------------------------------
Counts of SOMAs having cell_type

obs_label cell_type
                                                    count
name                                                     
B cell                                                  5
platelet                                                5
dendritic cell                                          4
mature NK T cell                                        4
neutrophil                                              4
plasma cell                                             4
animal cell                                             3
myeloid cell                                            3
natural killer cell                                     2
plasmablast                                             2
plasmacytoid dendritic cell                             2
erythrocyte                                             2
CD4-positive, alpha-beta T cell                         1
CD8-positive, alpha-beta T cell                         1
alpha-beta T cell                                       1
monocyte                                                1
eukaryotic cell                                         1
epithelial cell of nephron                              1
leukocyte                                               1
mesenchymal stem cell                                   1
native cell                                             1
CD14-low, CD16-positive monocyte                        1
CD14-positive monocyte                                  1
CD16-negative, CD56-bright natural killer cell,...      1
CD16-positive, CD56-dim natural killer cell, human      1
CD4-positive, alpha-beta memory T cell                  1
CD8-positive, alpha-beta memory T cell                  1
T cell                                                  1
conventional dendritic cell                             1
gamma-delta T cell                                      1
hematopoietic stem cell                                 1
immature B cell                                         1
memory B cell                                           1
mucosal invariant T cell                                1
naive B cell                                            1
naive thymus-derived CD4-positive, alpha-beta T...      1
naive thymus-derived CD8-positive, alpha-beta T...      1
regulatory T cell                                       1

Counts of SOMAs having tissue

obs_label tissue
        count
name         
blood       5
kidney      1
```

```
----------------------------------------------------------------
Collection-wide counts of values of cell_type

obs_label cell_type
                                                    count
name
monocyte                                            29878
naive thymus-derived CD4-positive, alpha-beta T...  26887
CD14-positive monocyte                              23648
myeloid cell                                        10261
naive B cell                                         8679
CD8-positive, alpha-beta T cell                      8658
B cell                                               8524
mature NK T cell                                     7755
CD16-positive, CD56-dim natural killer cell, human   6948
CD4-positive, alpha-beta T cell                      6726
CD8-positive, alpha-beta memory T cell               6224
erythrocyte                                          3918
CD16-negative, CD56-bright natural killer cell,...   3638
natural killer cell                                  3474
CD4-positive, alpha-beta memory T cell               3276
platelet                                             2926
mesenchymal stem cell                                2811
naive thymus-derived CD8-positive, alpha-beta T...   2387
CD14-low, CD16-positive monocyte                     1923
plasmablast                                          1825
T cell                                               1697
alpha-beta T cell                                    1659
epithelial cell of nephron                           1216
dendritic cell                                       1061
plasma cell                                          1025
animal cell                                           661
plasmacytoid dendritic cell                           650
conventional dendritic cell                           543
native cell                                           465
memory B cell                                         457
regulatory T cell                                     306
neutrophil                                            302
hematopoietic stem cell                               270
immature B cell                                       255
mucosal invariant T cell                              223
gamma-delta T cell                                    216
leukocyte                                             144
eukaryotic cell                                        28

TOTAL count    181544
dtype: int64
Collection-wide counts of values of tissue

obs_label tissue
         count
name
blood   176908
kidney    4636

TOTAL count    181544
dtype: int64
Collection-wide counts of values of cell_type_ontology_term_id

obs_label cell_type_ontology_term_id
            count
name
CL:0000576  29878
CL:0000895  26887
CL:0001054  23648
CL:0000763  10261
CL:0000788   8679
CL:0000625   8658
CL:0000236   8524
CL:0000814   7755
CL:0000939   6948
CL:0000624   6726
CL:0000909   6224
CL:0000232   3918
CL:0000938   3638
CL:0000623   3474
CL:0000897   3276
CL:0000233   2926
CL:0000134   2811
CL:0000900   2387
CL:0002396   1923
CL:0000980   1825
CL:0000084   1697
CL:0000789   1659
CL:1000449   1216
CL:0000451   1061
CL:0000786   1025
CL:0000548    661
CL:0000784    650
CL:0000990    543
CL:0000003    465
CL:0000787    457
CL:0000815    306
CL:0000775    302
CL:0000037    270
CL:0000816    255
CL:0000940    223
CL:0000798    216
CL:0000738    144
CL:0000255     28

TOTAL count    181544
dtype: int64
...
```

## Conclusion

From these we conclude that `obs.cell_type == "B cell"` and `obs.tissue == "blood"`, and
`var.feature_name == "MT-CO3"` (acquired similarly but not shown here) are likeliest to produce the
largest result set, given our local-disk mini-corpus.
