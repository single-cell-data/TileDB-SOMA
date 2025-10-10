# Developer Tools for the TileDB-SOMA R API

This directory provides a series of scripts to aid developers in maintaining the TileDB-SOMA R API

### Documentation for Developer Scripts

 - [`fetch_releases.R`](#fetch_releases.R)

## `fetch_releases.R`

This script fetches the GitHub releases of TileDB-SOMA and stores the versions and release dates within the R package. This information is used to manage [deprecations](https://github.com/single-cell-data/TileDB-SOMA/blob/main/dev_docs/POLICIES.md) in the TileDB-SOMA R API, as the R API includes both a deprecation stage _and_ a defunct stage. This information is stores in `inst/extdata/releases.dcf` as related to the root of the R package, and is accessible from the R package with

```r
system.file("extdata", "releases.dcf", package = "tiledbsoma")
```

## Dependencies

This script depends on the following R packages

 - [gh](https://cran.r-project.org/package=gh)
 - [rprojroot](https://cran.r-project.org/package=rprojroot)

This script also requires a GitHub token; for more details about setting up a GitHub token for use with the gh package, please see [`?gh::gh_token`](https://gh.r-lib.org/reference/gh_token.html)

### Usage

From the root the R package, run

```shell
Rscript tools/fetch_releases.R
```
