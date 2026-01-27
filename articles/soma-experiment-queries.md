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
#>   uri: /tmp/Rtmpdf17at/soma-exp-pbmc-small
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
#>   uri: file:///tmp/Rtmpdf17at/soma-exp-pbmc-small/ms
```

To use larger (or smaller) buffer sizes:

``` r
ctx <- SOMATileDBContext$new(c(soma.init_buffer_bytes = as.character(2 * 1024**3)))
#> Warning: SOMATileDBContext$new() was deprecated in tiledbsoma 2.3.0.
#> ℹ Use `SOMAContext` instead.
#> ℹ The deprecated feature was likely used in the R6 package.
#>   Please report the issue at <https://github.com/r-lib/R6/issues>.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
experiment <- SOMAExperimentOpen(experiment$uri, tiledbsoma_ctx = ctx)
#> Warning: SOMAExperimentOpen(tiledbsoma_ctx) was deprecated in tiledbsoma 2.3.0.
#> ℹ Use `context` instead.
#> ℹ The deprecated feature was likely used in the tiledbsoma package.
#>   Please report the issue at
#>   <https://github.com/single-cell-data/TileDB-SOMA/issues>.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
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
#> # A data frame: 80 × 9
#>    soma_joinid orig.ident  nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents
#>          <int> <fct>            <dbl>        <int> <fct>           <fct>        
#>  1           0 SeuratProj…         70           47 0               A            
#>  2           1 SeuratProj…         85           52 0               A            
#>  3           2 SeuratProj…         87           50 1               B            
#>  4           3 SeuratProj…        127           56 0               A            
#>  5           4 SeuratProj…        173           53 0               A            
#>  6           5 SeuratProj…         70           48 0               A            
#>  7           6 SeuratProj…         64           36 0               A            
#>  8           7 SeuratProj…         72           45 0               A            
#>  9           8 SeuratProj…         52           36 0               A            
#> 10           9 SeuratProj…        100           41 0               A            
#> # ℹ 70 more rows
#> # ℹ 3 more variables: groups <chr>, RNA_snn_res.1 <fct>, obs_id <chr>
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
#> # A data frame: 26 × 2
#>    obs_id         nCount_RNA
#>    <chr>               <dbl>
#>  1 TGACTGGATTCTCA        127
#>  2 AGTCAGACTGCACA        173
#>  3 AGAGATGATCTCGC        191
#>  4 GGGTAACTCTAGTG        101
#>  5 CTAAACCTGTGCAT        168
#>  6 TTGGTACTGAATCC        135
#>  7 TACATCACGCTAAC        109
#>  8 TTACCATGAATCGC        298
#>  9 ATAGGAGAAACAGA        406
#> 10 GCGCACGACTTTAC        213
#> # ℹ 16 more rows
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
