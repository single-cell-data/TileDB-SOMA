#!/usr/bin/env Rscript

## todo: chipset for macOS to detect arm
isX86 <- Sys.info()["machine"] == "x86_64"
isMac <- Sys.info()['sysname'] == "Darwin"
isLinux <- Sys.info()[["sysname"]] == "Linux"
macosver <- ""

if (isMac) {
    if (isX86) {
        url <- "https://github.com/TileDB-Inc/TileDB/releases/download/2.11.1/tiledb-macos-x86_64-2.11.1-15a1161.tar.gz"
        macosver <- "-mmacosx-version-min=10.14"
    } else {
        url <- "https://github.com/TileDB-Inc/TileDB/releases/download/2.11.1/tiledb-macos-arm64-2.11.1-15a1161.tar.gz"
    }
} else if (isLinux) {
    url <- "https://github.com/TileDB-Inc/TileDB/releases/download/2.11.1/tiledb-linux-x86_64-2.11.1-15a1161.tar.gz"
} else {
    stop("Unsupported platform for downloading artifacts. Please have TileDB Core installed locally.")
}

tarball <- "tiledb.tar.gz"
if (!file.exists(tarball)) download.file(url, tarball, quiet=TRUE)
if (!dir.exists("tiledb")) untar(tarball, exdir="tiledb")
if (!dir.exists("inst/tiledb")) untar(tarball, exdir="inst/tiledb")
archincl <- paste0("-I", system.file("include", package="arch"))

mkvar <- readLines("src/Makevars.in")
mkvar <- gsub("@cxx17_macos@", macosver, mkvar)
mkvar <- gsub("@tiledb_include@", paste("-I../inst/tiledb/include", archincl), mkvar)
mkvar <- gsub("@tiledb_libs@", "-ltiledb -L../inst/tiledb/lib", mkvar)
mkvar <- gsub("@tiledb_rpath@", "-Wl,-rpath,'$$ORIGIN/../tiledb/lib'", mkvar)
writeLines(mkvar, "src/Makevars")
