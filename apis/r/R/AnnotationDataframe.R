#' Single-cell Annotation Data Frame
#'
#' Base class for data frames with rows aligned to the observations or features
#' of a [`SOMA`]. Used to store a heterogeneous collection of
#' annotations/measurements.
#' @export

AnnotationDataframe <- R6::R6Class(
  classname = "AnnotationDataframe",
  inherit = AnnotationArray,

  public = list(
    #' @field verbose Print status messages
    verbose = TRUE,

    #' @description Retrieves the values from the array's dimension.
    ids = function() {
      arr <- self$object

      # Reset query layout to avoid (sc18770)
      tiledb::query_layout(arr) <- character()
      tiledb::attrs(arr) <- NA_character_
      arr[][[1]]
    },

    #' @description Ingest annotation data
    #' @param x a [`data.frame`]
    #' @param index_col Name to use for the TileDB array's dimension that will
    #' contain the data.frame row names.
    from_dataframe = function(x, index_col) {
      stopifnot(
        "'x' must be a data.frame" = is.data.frame(x),
        "'x' must have character row names" = has_character_rownames(x),
        "'index_col' must be a scalar character" = is_scalar_character(index_col)
      )

      # convert rownames to a column
      x[[index_col]] <- rownames(x)
      if (!self$exists()) {
        # TODO: Replace with configurable SOMAOptions class
        capacity <- switch(basename(self$uri),
          obs = 256L,
          var = 2048L,
          10000L
        )
        private$create_empty_array(x, index_col, capacity = capacity)
      } else {
        if (self$verbose) {
          message(
            sprintf("Updating existing %s at '%s'", self$class(), self$uri)
          )
        }
      }
      # Check for legacy validity mode metadata tag
      toggle_tiledb_legacy_mode_if_needed(self$object, self$verbose)
      private$ingest_data(x)
    },

    #' @description Retrieve the annotation data from TileDB
    #' @param attrs A character vector of the attribute names to retrieve. By
    #' default, all attributes are retrieved.
    #' @return A [`data.frame`] with row names containing values from the index
    #'    dimension
    to_dataframe = function(attrs = NULL) {
      if (self$verbose) {
        message(
          sprintf("Reading %s into memory from '%s'", self$class(), self$uri)
        )
      }

      # Check for legacy validity mode metadata tag
      toggle_tiledb_legacy_mode_if_needed(self$object, self$verbose)

      arr <- self$object
      tiledb::attrs(arr) <- attrs %||% character()
      tiledb::return_as(arr) <- "data.frame"
      df <- arr[]

      dimname <- self$dimnames()
      data.frame(
        df[setdiff(colnames(df), dimname)],
        row.names = df[[dimname]]
      )
    }
  )
)
