#!/usr/bin/env Rscript

if ( ! dir.exists("inst/") ) {
    stop("No 'inst/' directory. Exiting.", call. = FALSE)
}

filename="libtiledbsoma.tar.gz"
if ( ! file.exists(filename) ) {
    stop("No ", filename, ". Exiting.", call. = FALSE)
}

cat("** Using", filename, "\n")
if (!dir.exists("inst/tiledbsoma")) {
    cat("** Unarchiving", filename, "\n")
    untar(filename, exdir="inst/tiledbsoma")
}
