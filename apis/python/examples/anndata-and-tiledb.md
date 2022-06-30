If you're familiar with AnnData, you'll recognize much about the SOMA data model.

:::{.callout-tip}
See also the
[public TileDB Cloud notebook](https://cloud.tiledb.com/notebooks/details/johnkerl-tiledb/d3d7ff44-dc65-4cd9-b574-98312c4cbdbd/preview).
:::

## Single AnnData file

```
>>> import anndata

>>> ann = anndata.read_h5ad('pbmc3k_processed.h5ad')

>>> ann.obs.keys()
Index(['n_genes', 'percent_mito', 'n_counts', 'louvain'], dtype='object')

>>> ann.obs.shape
(2638, 4)

>>> ann.obs
                  n_genes  percent_mito  n_counts          louvain
index
AAACATACAACCAC-1      781      0.030178    2419.0      CD4 T cells
AAACATTGAGCTAC-1     1352      0.037936    4903.0          B cells
AAACATTGATCAGC-1     1131      0.008897    3147.0      CD4 T cells
AAACCGTGCTTCCG-1      960      0.017431    2639.0  CD14+ Monocytes
AAACCGTGTATGCG-1      522      0.012245     980.0         NK cells
...                   ...           ...       ...              ...
TTTCGAACTCTCAT-1     1155      0.021104    3459.0  CD14+ Monocytes
TTTCTACTGAGGCA-1     1227      0.009294    3443.0          B cells
TTTCTACTTCCTCG-1      622      0.021971    1684.0          B cells
TTTGCATGAGAGGC-1      454      0.020548    1022.0          B cells
TTTGCATGCCTCAC-1      724      0.008065    1984.0      CD4 T cells

[2638 rows x 4 columns]

>>> ann.obs.groupby(["louvain"]).size()
louvain
CD4 T cells          1144
CD14+ Monocytes       480
B cells               342
CD8 T cells           316
NK cells              154
FCGR3A+ Monocytes     150
Dendritic cells        37
Megakaryocytes         15
dtype: int64

>>> ann.var
         n_cells
index
TNFRSF4      155
CPSF3L       202
ATAD3C         9
C1orf86      501
RER1         608
...          ...
ICOSLG        34
SUMO3        570
SLC19A1       31
S100B         94
PRMT2        588

[1838 rows x 1 columns]

>>> ann.X
array([[-0.17146951, -0.28081203, -0.04667679, ..., -0.09826884,
        -0.2090951 , -0.5312034 ],
       [-0.21458222, -0.37265295, -0.05480444, ..., -0.26684406,
        -0.31314576, -0.5966544 ],
       [-0.37688747, -0.2950843 , -0.0575275 , ..., -0.15865596,
        -0.17087643,  1.379     ],
       ...,
       [-0.2070895 , -0.250464  , -0.046397  , ..., -0.05114426,
        -0.16106427,  2.0414972 ],
       [-0.19032837, -0.2263336 , -0.04399938, ..., -0.00591773,
        -0.13521303, -0.48211113],
       [-0.33378917, -0.2535875 , -0.05271563, ..., -0.07842438,
        -0.13032717, -0.4713379 ]], dtype=float32)
```

## Single TileDB SOMA

```
>>> import tiledbsc

>>> soma = tiledbsc.SOMA('tiledb://johnkerl-tiledb/pbmc3k_processed')

>>> soma
Name:    pbmc3k_processed
URI:     tiledb://johnkerl-tiledb/pbmc3k_processed
(n_obs, n_var): (2638, 1838)
X:       'data'
obs:     'n_genes', 'percent_mito', 'n_counts', 'louvain'
var:     'n_cells'
obsm:    'X_pca', 'X_tsne', 'X_umap', 'X_draw_graph_fr'
varm:    'PCs'
obsp:    'distances', 'connectivities'
varp:    
raw/X:   'data'
raw/var: 'n_cells'
uns:     draw_graph, louvain, louvain_colors, neighbors, pca

>>> soma.obs.keys()
['n_genes', 'percent_mito', 'n_counts', 'louvain']

>>> soma.obs.shape()
(2638, 4)

>>> soma.obs.df()
                  n_genes  percent_mito  n_counts          louvain
obs_id
AAACATACAACCAC-1      781      0.030178    2419.0      CD4 T cells
AAACATTGAGCTAC-1     1352      0.037936    4903.0          B cells
AAACATTGATCAGC-1     1131      0.008897    3147.0      CD4 T cells
AAACCGTGCTTCCG-1      960      0.017431    2639.0  CD14+ Monocytes
AAACCGTGTATGCG-1      522      0.012245     980.0         NK cells
...                   ...           ...       ...              ...
TTTCGAACTCTCAT-1     1155      0.021104    3459.0  CD14+ Monocytes
TTTCTACTGAGGCA-1     1227      0.009294    3443.0          B cells
TTTCTACTTCCTCG-1      622      0.021971    1684.0          B cells
TTTGCATGAGAGGC-1      454      0.020548    1022.0          B cells
TTTGCATGCCTCAC-1      724      0.008065    1984.0      CD4 T cells

[2638 rows x 4 columns]

>>> soma.obs.df(attrs=["louvain"]).groupby(["louvain"]).size()
louvain
B cells               342
CD14+ Monocytes       480
CD4 T cells          1144
CD8 T cells           316
Dendritic cells        37
FCGR3A+ Monocytes     150
Megakaryocytes         15
NK cells              154
dtype: int64  

>>> soma.var.df()
        n_cells
var_id
AAGAB        93
AAR2         63
AATF        294
ABCB1        29
ABCC10       30
...         ...
ZRANB3       16
ZSWIM6       19
ZUFSP        54
ZWINT        19
ZYX         403

[1838 rows x 1 columns]

>>> soma.X["data"].df()
                            value
obs_id           var_id
AAACATACAACCAC-1 AAGAB  -0.186726
                 AAR2   -0.142754
                 AATF   -0.334397
                 ABCB1  -0.126290
                 ABCC10 -0.139126
...                           ...
TTTGCATGCCTCAC-1 ZRANB3 -0.019673
                 ZSWIM6 -0.069497
                 ZUFSP  -0.167995
                 ZWINT  -0.096189
                 ZYX    -0.264577

[4848644 rows x 1 columns]
```

## See also

The convenience scripts [../tools/desc-ann](../tools/desc-ann) and [../tools/desc-soma](../tools/desc-soma)
display useful information (beyond `h5ls` etc) about dtypes, matrix types (CSR vs `numpy.ndarray`), etc.
