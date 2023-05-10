#!/usr/bin/env Rscript

isX86 <- Sys.info()["machine"] == "x86_64"
isMac <- Sys.info()['sysname'] == "Darwin"
isLinux <- Sys.info()[["sysname"]] == "Linux"
macosver <- ""

if (isMac) {
    if (isX86) {
        url <- "https://github.com/TileDB-Inc/TileDB/releases/download/2.15.2/tiledb-macos-x86_64-2.15.2-90f30eb.tar.gz"
        macosver <- "-mmacosx-version-min=10.14"
    } else {
        url <- "https://github.com/TileDB-Inc/TileDB/releases/download/2.15.2/tiledb-macos-arm64-2.15.2-90f30eb.tar.gz"
    }
} else if (isLinux) {
    url <- "https://github.com/TileDB-Inc/TileDB/releases/download/2.15.2/tiledb-linux-x86_64-2.15.2-90f30eb.tar.gz"
} else {
    stop("Unsupported platform for downloading artifacts. Please have TileDB Core installed locally.")
}

if (any(Sys.which(c("cmake", "git")) == ""))
    stop("The 'cmake' and 'git' commands are required.")

tarball <- "tiledb.tar.gz"
if (!file.exists(tarball)) download.file(url, tarball, quiet=TRUE)
if (!dir.exists("inst/tiledb")) untar(tarball, exdir="inst/tiledb")

soma_url <- "https://github.com/single-cell-data/TileDB-SOMA/archive/refs/tags/1.2.3.tar.gz"
soma_tarball <- "tiledbsoma.tar.gz"
if (!file.exists(soma_tarball)) download.file(soma_url, soma_tarball, quiet=TRUE)
if (!dir.exists("inst/tiledbsoma")) untar(soma_tarball, exdir="inst/tiledbsoma")
cat("** building libtiledbsoma\n")
system("tools/build_libtiledbsoma.sh")

mkvar <- readLines("src/Makevars.in")
mkvar <- gsub("@cxx17_macos@", macosver, mkvar)
mkvar <- gsub("@tiledb_include@", "-I../inst/tiledb/include -I../inst/tiledbsoma/include", mkvar)
mkvar <- gsub("@tiledb_libs@", "-ltiledb -L../inst/tiledb/lib -ltiledbsoma -L../inst/tiledbsoma/lib", mkvar)
mkvar <- gsub("@tiledb_rpath@", "-Wl,-rpath,'$$ORIGIN/../tiledb/lib' -Wl,-rpath,'$$ORIGIN/../tiledbsoma/lib'", mkvar)
writeLines(mkvar, "src/Makevars")
