#!/usr/bin/env Rscript

## version pinning info (NB: temporary from test repo url below)
tiledb_soma_version <- "0.0.1.27"
tiledb_soma_sha1 <- "76172ec"
##tiledb_soma_repo <- "single-cell-data/TileDB-SOMA"
tiledb_soma_repo <- "eddelbuettel/tldbsm2"

if ( ! dir.exists("inst/") ) {
    stop("No 'inst/' directory. Exiting.", call. = FALSE)
}

makeUrl <- function(arch, repo=tiledb_soma_repo, ver=tiledb_soma_version, sha1=tiledb_soma_sha1) {
    sprintf("https://github.com/%s/releases/download/%s/libtiledbsoma-%s-%s-%s.tar.gz", repo, ver, arch, ver, sha1)
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

filename <- basename(url)
cat("** Accessing", url, "as", filename, "\n")
if (!file.exists(filename)) {
    cat("** Downloading", url, "\n")
    download.file(url, filename, quiet=TRUE)
}
if (!dir.exists("inst/tiledbsoma")) {
    cat("** Unarchiving", filename, "\n")
    untar(filename, exdir="inst/tiledbsoma")
}
