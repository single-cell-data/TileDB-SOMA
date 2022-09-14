# tiledbsc

<!-- badges: start -->
[![R-CMD-check](https://github.com/TileDB-Inc/tiledbsc/workflows/R-CMD-check/badge.svg)](https://github.com/TileDB-Inc/tiledbsc/actions)
<!-- badges: end -->

An R implementation of the [Stack of Matrices, Annotated][soma-spec] (SOMA) API based on [TileDB](https://tiledb.com). SOMA is an open data model for representing annotated matrices, like those commonly used for single-cell data analysis. With tiledbsc, users can import from and export to in-memory formats used by popular toolchains like [Seurat][], [Bioconductor][bioc], and even [AnnData][] using the companion [Python package][tiledbsc-py].

## Installation

You can install the development version of *tiledbsc* from [GitHub](https://github.com/TileDB-Inc/tiledbsc) with:

``` r
# install.packages("remotes")
remotes::install_github("tiledb-inc/tiledbsc")
```

<!-- link -->
[tiledb]: https://tiledb.com
[soma-spec]: https://github.com/single-cell-data/SOMA
[seurat]: https://satijalab.org/seurat/
[bioc]: https://www.bioconductor.org/packages/release/bioc/html/Seurat.html
[bioc-se]: https://www.bioconductor.org/packages/SummarizedExperiment/
[bioc-sce]: https://www.bioconductor.org/packages/SingleCellExperiment/
[anndata]: https://anndata.readthedocs.io
[tiledbsc-py]: https://github.com/single-cell-data/TileDB-SingleCell
