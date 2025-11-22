# Querying a SOMA experiment

## Overview

In this notebook, we’ll take a quick look at how to query a
`SOMAExperiment` using the `SOMAExperimentAxisQuery` class. This allows
for easy selection of data from a `SOMAMeasurement` by filtering on
annotations stored in each axis data frame (i.e., `obs` and `var`).

``` r
library(tiledbsoma)
```

## Example data

Load the bundled `SOMAExperiment` containing a subsetted version of the
10X genomics [PBMC
dataset](https://satijalab.github.io/seurat-object/reference/pbmc_small.html)
provided by SeuratObject. This will return a `SOMAExperiment` object.

``` r
experiment <- load_dataset("soma-exp-pbmc-small")
experiment
#> <SOMAExperiment>
#>   uri: /tmp/RtmphvFyay/soma-exp-pbmc-small
```

## Querying basics

To perform a query we’ll need to initialize a new
`SOMAExperimentAxisQuery` object, specifying the `SOMAExperiment` and
the `SOMAMeasurement` within the experiment we want to query.

We can see that our current experiment contains only a single
measurement: `"RNA"`.

``` r
experiment$ms
#> <SOMACollection>
#>   uri: file:///tmp/RtmphvFyay/soma-exp-pbmc-small/ms
```

To use larger (or smaller) buffer sizes:

``` r
ctx <- SOMATileDBContext$new(c(soma.init_buffer_bytes = as.character(2 * 1024**3)))
experiment <- SOMAExperimentOpen(experiment$uri, tiledbsoma_ctx = ctx)
```

Alternatively, you can have in your environment
`export TILEDB_SOMA_INIT_BUFFER_BYTES=2147483648` before loading the
data.

Now we can construct our query object.

``` r
query <- SOMAExperimentAxisQuery$new(
  experiment = experiment,
  measurement_name = "RNA"
)
```

Once it’s created, we can use the `query` object to inspect, select, and
extract filtered data from the experiment.

For example, we can use `n_obs` and `n_vars` to determine the number of
observations and variables that passed our filtering criteria. Since we
didn’t specify any filtering criteria, these numbers will match the full
size of the experiment.

Number of observations:

``` r
query$n_obs
#> [1] 80
```

Number of variables:

``` r
query$n_vars
#> [1] 230
```

We can also extract any data component from the experiment. Here we’ll
read in the `obs` data frame from the query using `obs()` which returns
an iterator of
[`arrow::Table`](https://arrow.apache.org/docs/r/reference/Table-class.html).
The iterator is useful when the data is too large to load in memory
allowing to stream the data in chunks. This applies to
[`var()`](https://rdrr.io/r/stats/cor.html) as well.

To load the data in memory you can concatenate all chunks of the
iterator as shown below.

``` r
iterator <- query$obs()
obs <- iterator$concat()
obs
#> Table
#> 80 rows x 9 columns
#> $soma_joinid <int64 not null>
#> $orig.ident <dictionary<values=string, indices=int8>>
#> $nCount_RNA <double>
#> $nFeature_RNA <int32>
#> $RNA_snn_res.0.8 <dictionary<values=string, indices=int8>>
#> $letter.idents <dictionary<values=string, indices=int8>>
#> $groups <large_string>
#> $RNA_snn_res.1 <dictionary<values=string, indices=int8>>
#> $obs_id <large_string>
```

As a reminder `arrow:Table` can be easily cast into a `tibble`

``` r
obs$to_data_frame()
#>    soma_joinid    orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8
#> 1            0 SeuratProject         70           47               0
#> 2            1 SeuratProject         85           52               0
#> 3            2 SeuratProject         87           50               1
#> 4            3 SeuratProject        127           56               0
#> 5            4 SeuratProject        173           53               0
#> 6            5 SeuratProject         70           48               0
#> 7            6 SeuratProject         64           36               0
#> 8            7 SeuratProject         72           45               0
#> 9            8 SeuratProject         52           36               0
#> 10           9 SeuratProject        100           41               0
#> 11          10 SeuratProject         62           31               0
#> 12          11 SeuratProject        191           61               0
#> 13          12 SeuratProject        101           41               0
#> 14          13 SeuratProject         51           26               0
#> 15          14 SeuratProject         99           45               0
#> 16          15 SeuratProject        168           44               0
#> 17          16 SeuratProject         67           33               0
#> 18          17 SeuratProject        135           45               0
#> 19          18 SeuratProject         79           43               0
#> 20          19 SeuratProject        109           41               0
#> 21          20 SeuratProject        298           65               1
#> 22          21 SeuratProject        406           74               1
#> 23          22 SeuratProject        213           48               1
#> 24          23 SeuratProject        231           49               1
#> 25          24 SeuratProject        463           77               1
#> 26          25 SeuratProject         87           42               1
#> 27          26 SeuratProject        327           62               1
#> 28          27 SeuratProject        224           50               1
#> 29          28 SeuratProject        361           76               1
#> 30          29 SeuratProject        353           80               1
#> 31          30 SeuratProject        246           59               0
#> 32          31 SeuratProject        115           51               0
#> 33          32 SeuratProject        189           53               0
#> 34          33 SeuratProject        187           61               0
#> 35          34 SeuratProject        156           48               0
#> 36          35 SeuratProject        164           47               0
#> 37          36 SeuratProject        221           67               0
#> 38          37 SeuratProject        151           59               0
#> 39          38 SeuratProject        126           53               0
#> 40          39 SeuratProject        316           65               0
#> 41          40 SeuratProject        156           60               0
#> 42          41 SeuratProject        139           61               0
#> 43          42 SeuratProject        108           44               0
#> 44          43 SeuratProject         41           32               0
#> 45          44 SeuratProject        146           47               0
#> 46          45 SeuratProject        104           40               0
#> 47          46 SeuratProject        126           48               0
#> 48          47 SeuratProject         94           55               0
#> 49          48 SeuratProject        204           52               0
#> 50          49 SeuratProject         99           45               0
#> 51          50 SeuratProject        371           75               1
#> 52          51 SeuratProject        387           83               1
#> 53          52 SeuratProject        139           50               1
#> 54          53 SeuratProject         99           42               1
#> 55          54 SeuratProject        443           77               1
#> 56          55 SeuratProject        417           75               0
#> 57          56 SeuratProject        502           81               1
#> 58          57 SeuratProject        324           76               1
#> 59          58 SeuratProject        292           71               1
#> 60          59 SeuratProject        443           81               0
#> 61          60 SeuratProject        787           88               0
#> 62          61 SeuratProject        612           69               1
#> 63          62 SeuratProject        286           68               0
#> 64          63 SeuratProject        462           86               1
#> 65          64 SeuratProject        872           96               1
#> 66          65 SeuratProject        709           94               1
#> 67          66 SeuratProject        745           84               1
#> 68          67 SeuratProject        328           72               1
#> 69          68 SeuratProject        389           73               1
#> 70          69 SeuratProject        754           83               0
#> 71          70 SeuratProject        212           38               0
#> 72          71 SeuratProject        172           29               0
#> 73          72 SeuratProject        168           37               0
#> 74          73 SeuratProject        210           33               0
#> 75          74 SeuratProject        228           39               0
#> 76          75 SeuratProject        527           47               0
#> 77          76 SeuratProject        202           30               0
#> 78          77 SeuratProject        157           29               0
#> 79          78 SeuratProject        150           30               0
#> 80          79 SeuratProject        233           76               1
#>    letter.idents groups RNA_snn_res.1         obs_id
#> 1              A     g2             0 ATGCCAGAACGACT
#> 2              A     g1             0 CATGGCCTGTGCAT
#> 3              B     g2             0 GAACCTGATGAACC
#> 4              A     g2             0 TGACTGGATTCTCA
#> 5              A     g2             0 AGTCAGACTGCACA
#> 6              A     g1             0 TCTGATACACGTGT
#> 7              A     g1             0 TGGTATCTAAACAG
#> 8              A     g1             0 GCAGCTCTGTTTCT
#> 9              A     g1             0 GATATAACACGCAT
#> 10             A     g1             0 AATGTTGACAGTCA
#> 11             A     g2             2 AGGTCATGAGTGTC
#> 12             A     g1             2 AGAGATGATCTCGC
#> 13             A     g2             2 GGGTAACTCTAGTG
#> 14             A     g2             2 CATGAGACACGGGA
#> 15             A     g2             2 TACGCCACTCCGAA
#> 16             A     g1             2 CTAAACCTGTGCAT
#> 17             A     g2             2 GTAAGCACTCATTC
#> 18             A     g1             2 TTGGTACTGAATCC
#> 19             A     g1             2 CATCATACGGAGCA
#> 20             A     g2             2 TACATCACGCTAAC
#> 21             B     g1             1 TTACCATGAATCGC
#> 22             B     g1             1 ATAGGAGAAACAGA
#> 23             B     g2             1 GCGCACGACTTTAC
#> 24             B     g2             1 ACTCGCACGAAAGT
#> 25             B     g1             1 ATTACCTGCCTTAT
#> 26             B     g2             1 CCCAACTGCAATCG
#> 27             B     g2             1 AAATTCGAATCACG
#> 28             B     g2             1 CCATCCGATTCGCC
#> 29             B     g2             1 TCCACTCTGAGCTT
#> 30             B     g1             1 CATCAGGATGCACA
#> 31             A     g1             0 CTAAACCTCTGACA
#> 32             A     g1             2 GATAGAGAAGGGTG
#> 33             A     g1             0 CTAACGGAACCGAT
#> 34             A     g2             0 AGATATACCCGTAA
#> 35             A     g1             0 TACTCTGAATCGAC
#> 36             A     g1             0 GCGCATCTTGCTCC
#> 37             A     g2             0 GTTGACGATATCGG
#> 38             A     g1             0 ACAGGTACTGGTGT
#> 39             A     g1             0 GGCATATGCTTATC
#> 40             A     g2             0 CATTACACCAACTG
#> 41             A     g1             0 TAGGGACTGAACTC
#> 42             A     g2             2 GCTCCATGAGAAGT
#> 43             A     g2             0 TACAATGATGCTAG
#> 44             A     g2             0 CTTCATGACCGAAT
#> 45             A     g1             2 CTGCCAACAGGAGC
#> 46             A     g2             2 TTGCATTGAGCTAC
#> 47             A     g1             0 AAGCAAGAGCTTAG
#> 48             A     g2             0 CGGCACGAACTCAG
#> 49             A     g1             0 GGTGGAGATTACTC
#> 50             A     g2             0 GGCCGATGTACTCT
#> 51             B     g1             1 CGTAGCCTGTATGC
#> 52             B     g2             1 TGAGCTGAATGCTG
#> 53             B     g2             2 CCTATAACGAGACG
#> 54             B     g2             1 ATAAGTTGGTACGT
#> 55             B     g1             1 AAGCGACTTTGACG
#> 56             A     g1             1 ACCAGTGAATACCG
#> 57             B     g1             1 ATTGCACTTGCTTT
#> 58             B     g1             1 CTAGGTGATGGTTG
#> 59             B     g2             1 GCACTAGACCTTTA
#> 60             A     g1             0 CATGCGCTAGTCAC
#> 61             A     g1             2 TTGAGGACTACGCA
#> 62             B     g1             1 ATACCACTCTAAGC
#> 63             A     g1             2 CATATAGACTAAGC
#> 64             B     g1             1 TTTAGCTGTACTCT
#> 65             B     g1             2 GACATTCTCCACCT
#> 66             B     g2             1 ACGTGATGCCATGA
#> 67             B     g2             1 ATTGTAGATTCCCG
#> 68             B     g1             1 GATAGAGATCACGA
#> 69             B     g1             1 AATGCGTGGACGGA
#> 70             A     g1             2 GCGTAAACACGGTT
#> 71             A     g2             0 ATTCAGCTCATTGG
#> 72             A     g1             0 GGCATATGGGGAGT
#> 73             A     g2             0 ATCATCTGACACCA
#> 74             A     g2             0 GTCATACTTCGCCT
#> 75             A     g1             0 TTACGTACGTTCAG
#> 76             A     g1             0 GAGTTGTGGTAGCT
#> 77             A     g2             0 GACGCTCTCTCTCG
#> 78             A     g1             0 AGTCTTACTTCGGA
#> 79             A     g2             0 GGAACACTTCAGAC
#> 80             B     g1             1 CTTGATTGATCTTC
```

Alternatively, you can use the iterator, which retrieves data in chunks
that are smaller than the `soma.init_buffer_bytes` context field. You
can use the iterator’s method `$read_next()` to load a chunk in memory.

``` r
iterator <- query$obs()
iterator$read_next()
#> Table
#> 80 rows x 9 columns
#> $soma_joinid <int64 not null>
#> $orig.ident <dictionary<values=string, indices=int8>>
#> $nCount_RNA <double>
#> $nFeature_RNA <int32>
#> $RNA_snn_res.0.8 <dictionary<values=string, indices=int8>>
#> $letter.idents <dictionary<values=string, indices=int8>>
#> $groups <large_string>
#> $RNA_snn_res.1 <dictionary<values=string, indices=int8>>
#> $obs_id <large_string>
```

In this example the full `obs` table is relatively small and fits all in
one chunk.

For a bigger `SOMADataFrame` you can check if the iteration has finished
by checking the logical `$read_complete()`.

Here we demonstrate by creating a new iterator.

``` r
iterator <- experiment$obs$read()
iterator$read_complete()
#> [1] FALSE
```

``` r
iterator$read_next()
#> Table
#> 80 rows x 9 columns
#> $soma_joinid <int64 not null>
#> $orig.ident <dictionary<values=string, indices=int8>>
#> $nCount_RNA <double>
#> $nFeature_RNA <int32>
#> $RNA_snn_res.0.8 <dictionary<values=string, indices=int8>>
#> $letter.idents <dictionary<values=string, indices=int8>>
#> $groups <large_string>
#> $RNA_snn_res.1 <dictionary<values=string, indices=int8>>
#> $obs_id <large_string>
iterator$read_complete()
#> [1] TRUE
```

We can also access the expression via `X()`.

Similarly to `obs()` and [`var()`](https://rdrr.io/r/stats/cor.html),
`X()` is intended for iteration, but in this case we have access to two
different iterators, and thus `X()` returns a reader that gives you
access to an iterator for
[`arrow::Table`](https://arrow.apache.org/docs/r/reference/Table-class.html)
and one for `Matrix::sparse_matrix`.

Let’s take a look at the Arrow Table iterator:

``` r
reader <- query$X(layer_name = "counts")
table_irerator <- reader$tables()
table_irerator$read_next()
#> Table
#> 4456 rows x 3 columns
#> $soma_dim_0 <int64 not null>
#> $soma_dim_1 <int64 not null>
#> $soma_data <double not null>
```

As in the `obs` example the data is small enough to fit in one chunk.
For bigger data you can user `iterator$read_complete()` to check the
status of iteration and `iterator$concat()` to concatenate the rest of
the chunks.

The iterator for `Matrix::sparse_matrix` works in the same way. Keep in
mind that the matrix format is `dgTMatrix` as it is the most
memory-efficient and the only format type that can be easily iterated.
And most importantly, the resulting object is a “view” of the full
matrix with the original shape and indexes but only with data
corresponding to the query coordinates or filters (see section below).

``` r
reader <- query$X(layer_name = "counts")
iterator <- reader$sparse_matrix()
str(iterator$read_next())
#> Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
#>   ..@ i       : int [1:4456] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..@ j       : int [1:4456] 1 5 8 11 22 30 33 34 36 38 ...
#>   ..@ Dim     : int [1:2] 80 230
#>   ..@ Dimnames:List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>   ..@ x       : num [1:4456] 1 1 3 1 1 4 1 5 1 1 ...
#>   ..@ factors : list()
```

## Adding filters

Adding filters requires creating a `SOMAAxisQuery` object that allows
you to define coordinates, value filters, or both for an axis.

Here we’ll create a query for `obs` that slices the first 40 rows, and
then filters that subset based on the `nCount_RNA` column.

``` r
obs_query <- SOMAAxisQuery$new(
  coords = list(soma_joinid = 0:39),
  value_filter = "nCount_RNA > 100"
)
```

To apply this filter we’ll pass it to a new `SOMAExperimentAxisQuery`
object.

``` r
query <- SOMAExperimentAxisQuery$new(
  experiment = experiment,
  measurement_name = "RNA",
  obs_query = obs_query
)
```

Let’s see how many observations this query identified.

``` r
query$n_obs
#> [1] 26
```

As before, we can load the `obs` data frame into memory but now it only
includes the filtered observations.

``` r
obs <- query$obs(column_names = c("obs_id", "nCount_RNA"))$concat()
obs$to_data_frame()
#>            obs_id nCount_RNA
#> 1  TGACTGGATTCTCA        127
#> 2  AGTCAGACTGCACA        173
#> 3  AGAGATGATCTCGC        191
#> 4  GGGTAACTCTAGTG        101
#> 5  CTAAACCTGTGCAT        168
#> 6  TTGGTACTGAATCC        135
#> 7  TACATCACGCTAAC        109
#> 8  TTACCATGAATCGC        298
#> 9  ATAGGAGAAACAGA        406
#> 10 GCGCACGACTTTAC        213
#> 11 ACTCGCACGAAAGT        231
#> 12 ATTACCTGCCTTAT        463
#> 13 AAATTCGAATCACG        327
#> 14 CCATCCGATTCGCC        224
#> 15 TCCACTCTGAGCTT        361
#> 16 CATCAGGATGCACA        353
#> 17 CTAAACCTCTGACA        246
#> 18 GATAGAGAAGGGTG        115
#> 19 CTAACGGAACCGAT        189
#> 20 AGATATACCCGTAA        187
#> 21 TACTCTGAATCGAC        156
#> 22 GCGCATCTTGCTCC        164
#> 23 GTTGACGATATCGG        221
#> 24 ACAGGTACTGGTGT        151
#> 25 GGCATATGCTTATC        126
#> 26 CATTACACCAACTG        316
```

As well as the X matrix in two different formats:

[`arrow::Table`](https://arrow.apache.org/docs/r/reference/Table-class.html)

``` r
query$X("counts")$tables()$concat()
#> Table
#> 1485 rows x 3 columns
#> $soma_dim_0 <int64 not null>
#> $soma_dim_1 <int64 not null>
#> $soma_data <double not null>
```

`Matrix::sparse_matrix` in `dgTMatrix` format.

``` r
str(query$X("counts")$sparse_matrix()$concat())
#> Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
#>   ..@ i       : int [1:1485] 3 3 3 3 3 3 3 3 3 3 ...
#>   ..@ j       : int [1:1485] 8 11 22 30 31 32 33 34 37 38 ...
#>   ..@ Dim     : int [1:2] 80 230
#>   ..@ Dimnames:List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>   ..@ x       : num [1:1485] 13 1 6 5 2 1 2 2 1 1 ...
#>   ..@ factors : list()
```

For a re-indexed and re-shaped sparse matrix see the section below.

## Export to an annotated re-indexed sparse matrix

Any component of the queried `SOMAExperiment` can be exported to a
\[sparse matrix\]\[Matrix::sparseMatrix-class\] using the
`to_sparse_matrix()` method.

For example, let’s create a sparse matrix of the filtered expression
data. We’ll create a new query that returns a smaller subset of the data
to make the output easier to read.

``` r
query <- SOMAExperimentAxisQuery$new(
  experiment = experiment,
  measurement_name = "RNA",
  obs_query = SOMAAxisQuery$new(coords = 0:9),
  var_query = SOMAAxisQuery$new(coords = 0:9)
)
```

Then we indicate that we want to access the `"counts"` layer of the
`"X"` collection.

``` r
query$to_sparse_matrix(
  collection = "X",
  layer = "counts"
)
#> 10 x 10 sparse Matrix of class "dgTMatrix"
#>   [[ suppressing 10 column names '0', '1', '2' ... ]]
#>                       
#> 0 . 1 . . . 1 . .  3 .
#> 1 . . . 1 . . . .  7 .
#> 2 . . . . . . . . 11 .
#> 3 . . . . . . . . 13 .
#> 4 . . . 1 . . . .  3 .
#> 5 . . . 1 . . . .  4 .
#> 6 . . . . . . . .  6 .
#> 7 . . . 1 . . . .  4 .
#> 8 . . . . . . . .  2 .
#> 9 . 1 . . . . . . 21 .
```

By default, the dimensions are named using `soma_joinid` values which
are unique to each observation and variable. However, dimension names
can come from any column in the `obs` and `var` arrays that uniquely
identifies each record. For an expression matrix it makes sense to name
the dimensions using cell barcodes and gene names, which are stored in
the `obs_id` and `var_id` columns, respectively.

``` r
query$to_sparse_matrix(
  collection = "X",
  layer = "counts",
  obs_index = "obs_id",
  var_index = "var_id"
)
#> 10 x 10 sparse Matrix of class "dgTMatrix"
#>   [[ suppressing 10 column names 'MS4A1', 'CD79B', 'CD79A' ... ]]
#>                                    
#> ATGCCAGAACGACT . 1 . . . 1 . .  3 .
#> CATGGCCTGTGCAT . . . 1 . . . .  7 .
#> GAACCTGATGAACC . . . . . . . . 11 .
#> TGACTGGATTCTCA . . . . . . . . 13 .
#> AGTCAGACTGCACA . . . 1 . . . .  3 .
#> TCTGATACACGTGT . . . 1 . . . .  4 .
#> TGGTATCTAAACAG . . . . . . . .  6 .
#> GCAGCTCTGTTTCT . . . 1 . . . .  4 .
#> GATATAACACGCAT . . . . . . . .  2 .
#> AATGTTGACAGTCA . 1 . . . . . . 21 .
```

We can use this method for any of the `SOMAExperiment`’s collections.
Let’s access the t-SNE coordinates stored in the `obsm` collection’s
`X_tsne` layer, populating the row names with cell barcodes.

``` r
query$to_sparse_matrix(
  collection = "obsm",
  layer = "X_tsne",
  obs_index = "obs_id"
)
#> 10 x 2 sparse Matrix of class "dgTMatrix"
#>                          0           1
#> ATGCCAGAACGACT   0.8675977  -8.1007483
#> CATGGCCTGTGCAT  -7.3925306  -8.7717451
#> GAACCTGATGAACC -28.2064258   0.2410102
#> TGACTGGATTCTCA  16.3480689 -11.1633255
#> AGTCAGACTGCACA   1.9113998 -11.1929311
#> TCTGATACACGTGT   3.1475998  -9.9369312
#> TGGTATCTAAACAG  17.8526863  -9.8978901
#> GCAGCTCTGTTTCT  -6.4912187  -8.3874434
#> GATATAACACGCAT   1.3305297  -9.6807936
#> AATGTTGACAGTCA  16.9642732  -9.4295980
```

## Export to Seurat

The `query` object also contains methods for loading in results as a
Seurat object (or any of Seurat’s component classes). As with the
`to_sparse_matrix()` method, we can specify the `obs_index` and
`var_index` to use for naming the dimensions of the resulting object.

``` r
query <- SOMAExperimentAxisQuery$new(
  experiment = experiment,
  measurement_name = "RNA"
)

query$to_seurat(
  X_layers = c(counts = "counts", data = "data"),
  obs_index = "obs_id",
  var_index = "var_id"
)
#> Warning: Adding a command log without an assay associated with it
#> An object of class Seurat 
#> 230 features across 80 samples within 1 assay 
#> Active assay: RNA (230 features, 0 variable features)
#>  2 layers present: counts, data
#>  2 dimensional reductions calculated: pca, tsne
```
