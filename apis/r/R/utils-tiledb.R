#' Check if TileDB context contains a tag
#' @param key The name of the option
#' @return `TRUE` if the tag is set, `FALSE` otherwise
#' @noRd
tiledb_ctx_has_key <- function(key) {
  cfg <- tiledb::config(tiledb::tiledb_get_context())
  unname(!is.na(cfg[key]))
}

#' Retrieve the value of a TileDB context option
#' @param key The name of the option
#' @param default The default value to return if the option is not set
#' @return The value of the option, or `default` if the option is not set
#' @noRd
tiledb_ctx_get_key <- function(key, default = NULL) {
  if (!tiledb_ctx_has_key(key)) return(default)
  unname(tiledb::config(tiledb::tiledb_get_context())[key])
}

#' Set the value of a TileDB context option
#' @param key The name of the option
#' @param value The value to set
#' @noRd
tiledb_ctx_set_key <- function(key, value) {
  cfg <- tiledb::config(tiledb::tiledb_get_context())
  cfg[key] <- value
  invisible(tiledb::tiledb_ctx(config = cfg))
}

#' Toggle legacy validity mode in the global context
#'
#' Legacy validity mode is enabled if:
#' - the array does not contain the metadata tag `SOMA_LEGACY_VALIDITY_KEY`
#' - the array contains the metadata tag `"soma_legacy_validity"` and it is set
#'   to `"true"`
#'
#' Legacy validity mode is disabled if:
#' - `"r.legacy_validity_mode"`` is `"true"` in the global context and the array
#'   contains the metadata tag `SOMA_LEGACY_VALIDITY_KEY` and it is set to
#'   `"false"`
#'
#' Note this has to open and close the array to read the metadata.
#'
#' @param arr A [`tiledb::tiledb_array`]
#' @importFrom utils packageVersion
#' @noRd
toggle_tiledb_legacy_mode_if_needed <- function(arr, verbose = FALSE) {
  stopifnot(inherits(arr, "tiledb_array"))
  if (utils::packageVersion("tiledb") < "0.18.0.3") return()
  if (verbose) message(sprintf("Checking legacy validity mode for array: '%s'", arr@uri))

  # Get the metadata "soma_legacy_validity" tag
  tiledb::tiledb_array_open(arr, "READ")
  legacy_md_value <- tiledb::tiledb_get_metadata(arr, SOMA_LEGACY_VALIDITY_KEY)
  tiledb::tiledb_array_close(arr)

  # Has the user set the global context to use legacy mode?
  legacy_ctx_value <- tiledb_ctx_get_key(TILEDB_LEGACY_KEY)

  # Should we enable legacy mode?
  if (is.null(legacy_ctx_value)) { # Validity mode is unset in the global context
    if (is.null(legacy_md_value) || legacy_md_value == "true") {
      tiledb_ctx_set_key(TILEDB_LEGACY_KEY, "true")
      if (verbose) message("Enabled legacy validity mode")
    }
  } else { # Validity mode was previously set in the global context
    if (is.null(legacy_md_value)) {
      if (legacy_ctx_value == "false") {
        warning(
          sprintf("This array does not contain the '%s' metadata but legacy mode is disabled", SOMA_LEGACY_VALIDITY_KEY),
          call. = FALSE
        )
      }
    } else {
      if (legacy_ctx_value == "true" && legacy_md_value == "false") {
        stop("Legacy mode is enabled but this array was created without it", call. = FALSE)
      } else if (legacy_ctx_value == "false" && legacy_md_value == "true") {
        stop("Legacy mode is disabled but this array was created with it", call. = FALSE)
      }
    }
  }
}
