#' Single-cell Annotation Matrix
#'
#' Base class for matrix-like data with rows aligned to the observations or
#' features of a [`SOMA`].
#'
#' @param batch_mode logical, if `TRUE`, batch query mode is enabled, which
#' provides the ability to detect partial query results and resubmit until
#' all results are retrieved.
#' @export

AnnotationMatrix <- R6::R6Class(
  classname = "AnnotationMatrix",
  inherit = AnnotationArray,

  public = list(
    #' @field verbose Print status messages
    verbose = TRUE,

    #' @description Ingest annotation matrix
    #' @param x any `matrix`-like object coercible to a
    #' [`TsparseMatrix`][`Matrix::TsparseMatrix-class`] with string dimensions.
    #' @param index_col Name to use for the TileDB array's dimension that will
    #' contain the matrix row names.
    from_matrix = function(x, index_col) {
      stopifnot(
        "Must define 'index_col' to provide a dimension name" = !missing(index_col),
        "'index_col' must be a scalar character" = is_scalar_character(index_col)
      )
      private$validate_matrix(x)

      # convert to a data frame containing the index column
      x <- as.data.frame(x)
      x[[index_col]] <- rownames(x)

      if (!self$exists()) {
        private$create_empty_array(x, index_col)
      } else {
        if (self$verbose) {
          message(
            sprintf("Updating existing %s at '%s'", self$class(), self$uri)
          )
        }
      }
      private$ingest_data(x)
    },

    #' @description Retrieve the annotation data from TileDB
    #' @param attrs A character vector of the attribute names to retrieve. By
    #' default, all attributes are retrieved.
    #' @return A [`matrix`]
    to_matrix = function(attrs = NULL, batch_mode = FALSE) {

      df <- private$read_data(
        attrs = attrs,
        batch_mode = batch_mode,
        return_as = "data.frame"
      )

      index_col <- self$dimnames()
      attr_cols <- setdiff(colnames(df), index_col)

      mat <- as.matrix(df[attr_cols])
      rownames(mat) <- df[[index_col]]

      return(mat)
    }
  )
)
