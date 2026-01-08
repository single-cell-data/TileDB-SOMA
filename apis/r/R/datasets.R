#' SOMA Example Datasets
#'
#' Access example SOMA objects bundled with the tiledbsoma package.\cr
#' \cr
#' Use \code{list_datasets()} to list the available datasets and
#' \code{load_dataset()} to load a dataset into memory using the appropriate
#' SOMA class. The \code{extract_dataset()} function returns the path to the
#' extracted dataset without loading it into memory.
#'
#' The SOMA objects are stored as \code{tar.gz} files in the package's
#' \code{extdata} directory. Calling \code{load_dataset()} extracts the
#' \code{tar.gz} file to the specified \code{dir}, inspects its metadata to
#' determine the appropriate SOMA class to instantiate, and returns the
#' SOMA object.
#'
#' @name example-datasets
#' @rdname example-datasets
#'
NULL

#' @rdname example-datasets
#'
#' @return \code{list_datasets()}: returns a character vector of the
#' available datasets.
#'
#' @export
#'
#' @examples
#' list_datasets()
#'
list_datasets <- function() {
  data_dir <- example_data_dir()
  files <- dir(data_dir, pattern = "tar\\.gz$")
  return(tools::file_path_sans_ext(basename(files), compression = TRUE))
}

#' @rdname example-datasets
#'
#' @param name The name of the dataset.
#' @param dir The directory where the dataset will be extracted to
#' (default: \code{\link[base]{tempdir}()}).
#'
#' @return \code{extract_dataset()}: returns the path to the extracted dataset.
#'
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' dir <- withr::local_tempfile(pattern = "pbmc-small")
#' dir.create(dir, recursive = TRUE)
#' dest <- extract_dataset("soma-exp-pbmc-small", dir)
#' list.files(dest)
#'
extract_dataset <- function(name, dir = tempdir()) {
  data_dir <- example_data_dir()
  tarfiles <- list_datasets()

  stopifnot(
    "The specified directory does not exist" = dir.exists(dir),
    "Provide the name of a single dataset" = is_scalar_character(name),
    assert_subset(name, tarfiles, type = "dataset")
  )

  # Extract tar.gz file to dir
  tarfile <- dir(data_dir, pattern = name, full.names = TRUE)
  tarfile <- file.path(data_dir, sprintf("%s.tar.gz", name))
  stopifnot("The specified dataset does not exist" = file.exists(tarfile))

  dataset_uri <- file.path(dir, name)
  utils::untar(tarfile, exdir = dataset_uri)
  dataset_uri
}

#' @rdname example-datasets
#'
#' @param tiledbsoma_ctx Optional (DEPRECATED) TileDB \dQuote{Context} object
#' that defaults to \code{NULL}.
#' @param context Optional \code{SOMAContext} object used for TileDB operations.
#' If a context is not provided, then the default context will be used. Call
#' \code{set_default_context} once before other SOMA operations to configure
#' the default context.
#'
#' @return \code{load_dataset()}: returns a SOMA object.
#'
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' dir <- withr::local_tempfile(pattern = "pbmc_small")
#' dir.create(dir, recursive = TRUE)
#' (exp <- load_dataset("soma-exp-pbmc-small", dir))
#'
#' \dontshow{
#' exp$close()
#' }
#'
load_dataset <- function(name, dir = tempdir(), tiledbsoma_ctx = NULL, context = NULL) {
  dataset_uri <- extract_dataset(name, dir)

  # Inspect the object's metadata
  context <- get_soma_context(context, tiledbsoma_ctx, what="load_dataset(tiledbsoma_ctx)")
  metadata <- get_all_metadata(
    dataset_uri,
    is_array = switch(
      get_tiledb_object_type(dataset_uri, ctxxp = context$handle),
      ARRAY = TRUE,
      GROUP = FALSE,
      stop("The dataset is not a TileDB Array or Group", call. = FALSE)
    ),
    ctxxp = context$handle
  )
  return(switch(
    metadata$soma_object_type %||% "",
    SOMAExperiment = SOMAExperimentOpen(
      dataset_uri,
      tiledbsoma_ctx = tiledbsoma_ctx,
      context = context
    ),
    SOMADataFrame = SOMADataFrameOpen(
      dataset_uri,
      tiledbsoma_ctx = tiledbsoma_ctx,
      context = context
    ),
    stop("The dataset is an unsupported SOMA object", call. = FALSE)
  ))
}

example_data_dir <- function() {
  system.file("extdata", package = "tiledbsoma", mustWork = TRUE)
}
