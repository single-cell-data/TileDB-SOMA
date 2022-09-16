# tiledbsoma

<!-- badges: start -->
[![r-cmd-check](https://github.com/TileDB-Inc/tiledbsoma/workflows/r-cmd-check/badge.svg)](https://github.com/TileDB-Inc/tiledbsoma/actions)
<!-- badges: end -->

An R implementation of the [Stack of Matrices, Annotated][soma-spec] (SOMA) API based on [TileDB](https://tiledb.com). SOMA is an open data model for representing annotated matrices, like those commonly used for single-cell data analysis. With tiledbsoma, users can import from and export to in-memory formats used by popular toolchains like [Seurat][], [Bioconductor][bioc], and even [AnnData][] using the companion [Python package][tiledbsoma-py].

## Installation

You can install the development version of *tiledbsoma* from [GitHub](https://github.com/TileDB-Inc/tiledbsoma) with:

``` r
# install.packages("remotes")
remotes::install_github("tiledb-inc/tiledbsoma")
```

<!-- link -->
[tiledb]: https://tiledb.com
[soma-spec]: https://github.com/single-cell-data/SOMA
[seurat]: https://satijalab.org/seurat/
[bioc]: https://www.bioconductor.org/packages/release/bioc/html/Seurat.html
[bioc-se]: https://www.bioconductor.org/packages/SummarizedExperiment/
[bioc-sce]: https://www.bioconductor.org/packages/SingleCellExperiment/
[anndata]: https://anndata.readthedocs.io
[tiledbsoma-py]: https://github.com/single-cell-data/TileDB-SingleCell
