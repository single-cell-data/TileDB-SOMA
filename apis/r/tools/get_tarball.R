#!/usr/bin/env Rscript

res <- jsonlite::fromJSON("https://api.github.com/repos/tiledb-inc/tiledb/releases/latest")
cat("Seeing", res$name, "\n")

## todo: chipset for macOS to detect arm
isX86 <- Sys.info()["machine"] == "x86_64"
isMac <- Sys.info()['sysname'] == "Darwin"
isLinux <- Sys.info()[["sysname"]] == "Linux"
macosver <- ""

if (isMac) {
    if (isX86) {
        url <- res$assets$browser_download_url[grepl("macos-x86", res$assets$browser_download_url)]
        macosver <- "-mmacosx-version-min=10.14"
    } else {
        url <- res$assets$browser_download_url[grepl("macos-arm64", res$assets$browser_download_url)]
    }
} else if (isLinux) {
    ## use [1] as [2] is the noavx2 one
    url <- res$assets$browser_download_url[grepl("linux", res$assets$browser_download_url)][1]
} else {
    stop("Unsupported platform for downloading artifacts. Please have TileDB Core installed locally.")
}

tarball <- "tiledb.tar.gz"
if (!file.exists(tarball)) download.file(url, tarball, quiet=TRUE)
if (!dir.exists("tiledb")) untar(tarball, exdir="tiledb")
if (!dir.exists("inst/tiledb")) untar(tarball, exdir="inst/tiledb")
archincl <- paste0("-I", system.file("include", package="arch"))

mkvar <- readLines("src/Makevars.in")
mkvar <- gsub("@cxx20_macos@", macosver, mkvar)
mkvar <- gsub("@tiledb_include@", paste("-I../inst/tiledb/include", archincl), mkvar)
mkvar <- gsub("@tiledb_libs@", "-ltiledb -L../inst/tiledb/lib", mkvar)
mkvar <- gsub("@tiledb_rpath@", "-Wl,-rpath,'$$ORIGIN/../tiledb/lib'", mkvar)
writeLines(mkvar, "src/Makevars")
