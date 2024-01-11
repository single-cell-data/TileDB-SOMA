#!/usr/bin/env Rscript

## version pinning info
tiledb_core_version <- "2.19.0"
# 8-nybble hash for 2.19.0 only. Please see https://github.com/TileDB-Inc/TileDB/pull/4599.
tiledb_core_sha1 <- "fa30a88a"

if ( ! dir.exists("inst/") ) {
    stop("No 'inst/' directory. Exiting.", call. = FALSE)
}

makeUrl <- function(arch, ver=tiledb_core_version, sha1=tiledb_core_sha1) {
    sprintf("https://github.com/TileDB-Inc/TileDB/releases/download/%s/tiledb-%s-%s-%s.tar.gz", ver, arch, ver, sha1)
}

isMac <- Sys.info()["sysname"] == "Darwin"
isLinux <- Sys.info()["sysname"] == "Linux"

if (isMac) {
    arch <- system('uname -m', intern = TRUE)
    url <- makeUrl(paste0("macos-", arch))
} else if (isLinux) {
    url <- makeUrl("linux-x86_64")
} else {
    stop("Unsupported platform for downloading artifacts. Please have TileDB Core installed locally.")
}

tarball <- "tiledb.tar.gz"
if (!file.exists(tarball)) download.file(url, tarball, quiet=TRUE)
if (!dir.exists("inst/tiledb")) untar(tarball, exdir="inst/tiledb")
