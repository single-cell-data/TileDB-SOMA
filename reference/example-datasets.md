# SOMA Example Datasets

Access example SOMA objects bundled with the tiledbsoma package.  
  
Use `list_datasets()` to list the available datasets and
`load_dataset()` to load a dataset into memory using the appropriate
SOMA class. The `extract_dataset()` function returns the path to the
extracted dataset without loading it into memory.

## Usage

``` r
list_datasets()

extract_dataset(name, dir = tempdir())

load_dataset(name, dir = tempdir(), tiledbsoma_ctx = NULL, context = NULL)
```

## Arguments

- name:

  The name of the dataset.

- dir:

  The directory where the dataset will be extracted to (default:
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

- tiledbsoma_ctx:

  Optional (DEPRECATED) TileDB “Context” object that defaults to `NULL`.

- context:

  Optional `SOMAContext` object used for TileDB operations. If a context
  is not provided, then the default context will be used. Call
  `set_default_context` once before other SOMA operations to configure
  the default context.

## Value

`list_datasets()`: returns a character vector of the available datasets.

`extract_dataset()`: returns the path to the extracted dataset.

`load_dataset()`: returns a SOMA object.

## Details

The SOMA objects are stored as `tar.gz` files in the package's `extdata`
directory. Calling `load_dataset()` extracts the `tar.gz` file to the
specified `dir`, inspects its metadata to determine the appropriate SOMA
class to instantiate, and returns the SOMA object.

## Examples

``` r
list_datasets()
#> [1] "soma-dataframe-pbmc3k-processed-obs" "soma-exp-pbmc-small-pre-1.15"       
#> [3] "soma-exp-pbmc-small"                

dir <- withr::local_tempfile(pattern = "pbmc-small")
dir.create(dir, recursive = TRUE)
dest <- extract_dataset("soma-exp-pbmc-small", dir)
list.files(dest)
#> [1] "__group"            "__meta"             "__tiledb_group.tdb"
#> [4] "ms"                 "obs"                "uns"               
dir <- withr::local_tempfile(pattern = "pbmc_small")
dir.create(dir, recursive = TRUE)
(exp <- load_dataset("soma-exp-pbmc-small", dir))
#> <SOMAExperiment>
#>   uri: /tmp/RtmpbAgXbM/pbmc_small284672471e06/soma-exp-pbmc-small
```
