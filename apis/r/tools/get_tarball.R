#!/usr/bin/env Rscript

# Please use scripts/update-tiledb-version.py in this repo's base directory to
# update URLs in this file as well as the FindTileDB_EP.cmake file within the
# libtiledbsoma directory.

if ( ! dir.exists("inst/") ) {
    stop("No 'inst/' directory. Exiting.", call. = FALSE)
}

isMac <- Sys.info()["sysname"] == "Darwin"
isLinux <- Sys.info()["sysname"] == "Linux"

if (isMac) {
    arch <- system('uname -m', intern = TRUE)
    url <- "https://github.com/TileDB-Inc/TileDB/releases/download/2.19.1/tiledb-macos-x86_64-2.19.1-29ceb3e7.tar.gz"
} else if (isLinux) {
    url <- "https://github.com/TileDB-Inc/TileDB/releases/download/2.19.1/tiledb-linux-x86_64-2.19.1-29ceb3e7.tar.gz"
} else {
    stop("Unsupported platform for downloading artifacts. Please have TileDB Core installed locally.")
}

tarball <- "tiledb.tar.gz"
if (!file.exists(tarball)) download.file(url, tarball, quiet=TRUE)
if (!dir.exists("inst/tiledb")) untar(tarball, exdir="inst/tiledb")
