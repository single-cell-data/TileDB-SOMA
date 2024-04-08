## .onLoad is called whether code from the package is used and the packages is 'loaded'. An
## example is calling `tiledbsoma::show_package_versions()`. So this is most elementary check,
## .onAttach is also called when the package is 'attached' via 'library(tiledbsoma)'
## During package build and byte-code compilation and load check, both are called.
.onLoad <- function(libname, pkgname) {
    rpkg_lib_version <- tiledb::tiledb_version(compact=TRUE)
    soma_lib_version <- libtiledbsoma_version(compact=TRUE)
    if (rpkg_lib_version != soma_lib_version) {
        msg <- sprintf("TileDB Core version %s used by TileDB-R package, but TileDB-SOMA uses %s",
                       sQuote(rpkg_lib_version), sQuote(soma_lib_version))
        stop(msg, call. = FALSE)
    }
}

## A _possible_ .onAttach() function -- these are not allowed to use cat() etc but _must_ communicate
## via packageStartupMessage() as this function can be 'muzzled' as desired. See Writing R Extensions.
## .onAttach <- function(libname, pkgname) {
##     if (interactive()) {
##         packageStartupMessage("TileDB-SOMA R package ", packageVersion(pkgname),
##                               " with TileDB Embedded ", format(tiledb::tiledb_version(TRUE)),
##                               " on ", utils::osVersion,
##                               ".\nSee https://github.com/single-cell-data for more information ",
##                               "about the SOMA project.")
##     }
## }
