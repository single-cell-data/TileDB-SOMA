#!/usr/bin/env Rscript

## version pinning info
+tiledb_core_version <- "2.17.0-rc0"
+tiledb_core_sha1 <- "46b9ca5"

if ( ! dir.exists("inst/") ) {
    stop("No 'inst/' directory. Exiting.", call. = FALSE)
}

makeUrl <- function(arch, ver=tiledb_core_version, sha1=tiledb_core_sha1) {
    sprintf("https://github.com/TileDB-Inc/TileDB/releases/download/%s/tiledb-%s-%s-%s.tar.gz", ver, arch, ver, sha1)
}

isX86 <- Sys.info()["machine"] == "x86_64"
isMac <- Sys.info()["sysname"] == "Darwin"
isLinux <- Sys.info()["sysname"] == "Linux"

if (isMac && isX86) {
    url <- makeUrl("macos-x86_64")
} else if (isMac && !isX86) {
    url <- makeUrl("macos-arm64")
} else if (isLinux) {
    url <- makeUrl("linux-x86_64")
} else {
    stop("Unsupported platform for downloading artifacts. Please have TileDB Core installed locally.")
}

tarball <- "tiledb.tar.gz"
if (!file.exists(tarball)) download.file(url, tarball, quiet=TRUE)
if (!dir.exists("inst/tiledb")) untar(tarball, exdir="inst/tiledb")
