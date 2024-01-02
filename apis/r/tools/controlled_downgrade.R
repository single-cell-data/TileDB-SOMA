#!/usr/bin/env Rscript

## The TileDB-R package builds against TileDB Core, and while we can ensure a _minimum_ version
## we cannot ensure an exact version.  That can lead to cases as in late 2023 / early 2024 where
## TileDB-R 0.23.0 is at CRAN with a default of TileDB Core 2.19.0, yet TileDB-SOMA is still set
## up to build with TileDB Core 2.18.2.  This mismatch is undesirable but a little difficult to
## avoid for example in continuous integration where R will always go to the newest package version.
##
## As an experiment, Paul and I have set up two helper repos. The first ensures we can still build
## the _previous_ TileDB-R version at r-universe (which would otherwise pick newest versions from
## CRAN and/or GitHub).  The second deploys these builds in a particular R repository hosted on
## GitHun from which we can install.  Together, these provide a simple 'time machine' ability to
## the previous binary.
##
## To illustrate, by default we end with these packages (when running on Ubuntu 22.04 as for CI)
##   > show_package_versions()
##   tiledbsoma:    1.6.0
##   tiledb-r:      0.23.0
##   tiledb core:   2.19.0
##   libtiledbsoma: 2.18.2
##   R:             R version 4.3.2 (2023-10-31)
##   OS:            Ubuntu 22.04.1 LTS
##   >
##
## Whereas once we do a controlled downgrade we get
##
##   > library(tiledbsoma)
##   > show_package_versions()
##   tiledbsoma:    1.6.0
##   tiledb-r:      0.22.0
##   tiledb core:   2.18.2
##   libtiledbsoma: 2.18.2
##   R:             R version 4.3.2 (2023-10-31)
##   OS:            Ubuntu 22.04.1 LTS
##   >

if (!grepl("Ubuntu 22.04", utils::osVersion))   	# we currently only support Ubuntu 22.04
    q()

core_from_tiledb <- tiledb::tiledb_version(compact=TRUE)
core_from_soma <- tiledbsoma:::libtiledbsoma_version(compact=TRUE)

if (core_from_tiledb == core_from_soma) 			# nothing to do: versions align
    q()

if (core_from_tiledb > core_from_soma) { 			# if TileDB-R has newer TileDB Core
    if (requireNamespace("bspm", quietly=TRUE))		#   if we have bspm (with r2u, likely in CI)
        bspm::disable() 							#     disable it to permit direct binary install

    url <- "https://eddelbuettel.github.io/tiledb-r-versioned-drat/0.22.0/bin/linux/jammy/4.3/src/contrib"

    install.packages("tiledb", contriburl=url)
}
