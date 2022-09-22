#' Single-cell Annotation Pairwise Matrix
#'
#' Base class for matrix-like data storing pairwise relationships with rows and
#' columns aligned to either the observations or features of `X`
#' ([`AssayMatrix`]).
#'
#' @export

AnnotationPairwiseMatrix <- R6::R6Class(
  classname = "AnnotationPairwiseMatrix",
  inherit = AnnotationArray,

  public = list(

    #' @description Ingest annotation matrix
    #' @param x any `matrix`-like object coercible to a
    #' [`TsparseMatrix`][`Matrix::TsparseMatrix-class`] with string dimensions.
    #' @param index_cols Names to use for the TileDB array's dimensions that
    #' will contain the matrix row/column names.
    #' @param value_col Name to use for the TileDB array's attribute that will
    #' contain the matrix values.
    from_matrix = function(x, index_cols, value_col = "value") {
      stopifnot(
        "Must provide 'index_cols' to name the index columns" = !missing(index_cols),
        "'index_cols' must be a character vector of length 2" = length(index_cols) == 2,
        "'value_col' must be a scalar character" = is_scalar_character(value_col)
      )
      private$validate_matrix(x)

      x <- matrix_to_coo(
        x = x,
        index_cols = index_cols,
        value_col = value_col
      )

      if (!self$exists()) {
        private$create_empty_array(x, index_cols)
      } else {
        if (self$verbose) {
          message(
            sprintf("Updating existing %s at '%s'", self$class(), self$uri)
          )
        }
      }
      private$ingest_data(x)
    },

    #' @description Read annotation data from TileDB into a matrix
    #' @return A [`matrix`]
    to_matrix = function() {
      if (self$verbose) {
        message(
          sprintf("Reading %s into matrix from '%s'", self$class(), self$uri)
        )
      }
      self$tiledb_array(return_as = "matrix")[]
    },

    #' @description Read annotation data from TileDB into a data frame
    #' @param attrs Specify one or more attributes to retrieve. If `NULL`,
    #' all attributes are retrieved.
    #' @return A [`data.frame`]
    to_dataframe = function(attrs = NULL) {
      if (self$verbose) {
        message(
          sprintf("Reading %s into dataframe from '%s'", self$class(), self$uri)
        )
      }
      attrs <- attrs %||% character()
      self$tiledb_array(attrs = attrs, return_as = "data.frame")[]
    },

    #' @description Read annotation data from TileDB into a sparse matrix
    #' @return A [`Matrix::dgTMatrix-class`].
    to_sparse_matrix = function() {
      dataframe_to_dgtmatrix(
        self$to_dataframe(),
        index_cols = self$dimnames()
      )[[1]]
    },

    #' @description Read annotation data from TileDB into Seurat Graph
    #' @return A [`SeuratObject::Graph-class`]
    #' @importFrom SeuratObject DefaultAssay
    to_seurat_graph = function() {
      assay <- self$get_metadata(key = "assay_used")
      object <- SeuratObject::as.Graph(self$to_sparse_matrix())
      SeuratObject::DefaultAssay(object) <- assay
      object
    }
  )
)
