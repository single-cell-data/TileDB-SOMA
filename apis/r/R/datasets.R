#' SOMA Example Datasets
#'
#' Access example SOMA objects bundled with the tiledbsoma package. Use
#' `list_datasets()` to list the available datasets and `load_dataset()` to load
#'  a dataset into memory using the appropriate SOMA class. The
#'  `extract_dataset()` method returns the path to the extracted dataset without
#'  loading it into memory.
#'
#' @details
#' The SOMA objects are stored as `tar.gz` files in the package's `raw-data`
#' directory. Calling `load_dataset()` extracts the `tar.gz` file to the
#' specified `dir`, inspects its metadata to determine the appropriate SOMA
#' class to instantiate, and returns the SOMA object.
#'
#' @examples
#' soma_pbmc_small <- load_dataset("soma-exp-pbmc-small")
#'
#' @name example-datasets
NULL

#' @rdname example-datasets
#' @return A character vector of the available datasets.
#' @importFrom tools file_path_sans_ext
#' @export
list_datasets <- function() {
  data_dir <- system.file("raw-data", package = "tiledbsoma", mustWork = TRUE)
  stopifnot("The dataset directory does not exist" = dir.exists(data_dir))
  files <- dir(data_dir, pattern = "tar.gz")
  tools::file_path_sans_ext(basename(files), compression = TRUE)
}

#' @rdname example-datasets
#' @param name The name of the dataset.
#' @param dir The directory where the dataset will be extracted to (default:
#' `tempdir()`).
#' @return The path to the extracted dataset.
#' @export
extract_dataset <- function(name, dir = tempdir()) {
  data_dir <- system.file("raw-data", package = "tiledbsoma")
  tarfiles <- list_datasets()

  stopifnot(
    "The specified directory does not exist" = dir.exists(dir),
    "Provide the name of a single dataset" = is_scalar_character(name),
    assert_subset(name, tarfiles, type = "dataset")
  )

  # Extract tar.gz file to dir
  tarfile <- dir(data_dir, pattern = name, full.names = TRUE)
  stopifnot("The specified dataset does not exist" = file.exists(tarfile))

  dataset_uri <- file.path(dir, name)
  untar(tarfile, exdir = dataset_uri)
  dataset_uri
}

#' @rdname example-datasets
#' @inheritParams extract_dataset
#' @return An SOMA object.
#' @export
load_dataset <- function(name, dir = tempdir()) {
  dataset_uri <- extract_dataset(name, dir)

  # Inspect the object's metadata
  object <- switch(
    tiledb::tiledb_object_type(dataset_uri),
    "ARRAY" = TileDBArray$new(dataset_uri, internal_use_only = "allowed_use"),
    "GROUP" = TileDBGroup$new(dataset_uri, internal_use_only = "allowed_use"),
    stop("The dataset is not a TileDB Array or Group", call. = FALSE)
  )

  # Instantiate the proper SOMA object
  switch(
    object$get_metadata("soma_object_type"),
    "SOMAExperiment" = SOMAExperimentOpen(dataset_uri),
    "SOMADataFrame" = SOMADataFrameOpen(dataset_uri),
    stop("The dataset is an unsupported SOMA object", call. = FALSE)
  )
}
