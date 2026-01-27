#' @import nanoarrow
#' @importFrom Rcpp evalCpp
#' @importFrom methods as is
#' @importFrom Matrix as.matrix
#' @importFrom bit64 as.integer64
#' @importFrom arrow concat_arrays
#' @importFrom tiledb query_layout
#' @importFrom rlang is_scalar_logical is_scalar_character
#' @useDynLib tiledbsoma, .registration=TRUE
#'
NULL

.pkgenv <- new.env(parent = emptyenv())

## .onLoad is called whether code from the package is used and the packages is 'loaded'. An
## example is calling `tiledbsoma::show_package_versions()`. So this is most elementary check,
## .onAttach is also called when the package is 'attached' via 'library(tiledbsoma)'
## During package build and byte-code compilation and load check, both are called.
.onLoad <- function(libname, pkgname) {
  .pkgenv$libname <- libname
  .pkgenv$pkgname <- pkgname
  ## create a slot for somactx in per-package enviroment, do no fill it yet to allow 'lazy load'
  .pkgenv[["somactx"]] <- NULL
}

## An .onAttach() function is not allowed to use cat() etc but _must_ communicate via
## packageStartupMessage() as this function can be 'muzzled' as desired. See Writing R Extensions.
.onAttach <- function(libname, pkgname) {
  rpkg_lib <- get_tiledb_version(compact = FALSE)
  # Check major and minor but not micro: sc-50464
  rpkg_lib_version <- paste(rpkg_lib[["major"]], rpkg_lib[["minor"]], sep = ".")
  soma_lib_version <- libtiledbsoma_version(
    compact = TRUE,
    major_minor_only = TRUE
  )
  if (rpkg_lib_version != soma_lib_version) {
    msg <- sprintf(
      "TileDB Core version %s used by TileDB-R package, but TileDB-SOMA uses %s",
      sQuote(rpkg_lib_version),
      sQuote(soma_lib_version)
    )
    packageStartupMessage(msg)
  }
  if (interactive()) {
    packageStartupMessage(
      "TileDB-SOMA R package ",
      utils::packageVersion(pkgname),
      " with TileDB Embedded ",
      format(get_tiledb_version(TRUE)),
      " on ",
      utils::osVersion,
      ".\nSee https://github.com/single-cell-data for more information ",
      "about the SOMA project."
    )
  }
}

#' Create and cache a SOMA Context Object
#'
#' @param config A named character vector with \dQuote{key} and \dQuote{value}
#' pairs defining the configuration setting.
#'
#' @return An external pointer object containing a shared pointer instance
#' of \code{SOMAContext}.
#'
#' @keywords internal
#'
#' @export
#'
soma_context <- function(config = NULL) {
  .deprecate(when = "2.3.0", what = "soma_context()", details = "Use new `SOMAContext` class instead.")
  context <- SOMAContext$new(config)
  .pkgenv[["somactx"]] <- context
  return(SOMAContext$new(config))
}
