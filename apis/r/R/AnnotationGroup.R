#' Abstract Annotation Group
#'
#' @description
#' This is an abstract class to augment the base TiledbGroup class with fields
#' necessary for the child Annotation- classes. Currently it only adds the
#' `dimension_name` field.
AnnotationGroup <- R6::R6Class(
  classname = "AnnotationGroup",
  inherit = TileDBGroup,

  public = list(
    #' @field dimension_name Optional name of the dimension shared by all
    #' arrays within the group (typically `obs_id` or `var_id`).
    dimension_name = NULL,

    #' @description Create a new `TileDBGroup`-based Annotation class.
    #' @param uri URI for the TileDB group.
    #' @param dimension_name Optional name of the dimension shared by all
    #' of the group's member arrays.
    #' @param verbose Print status messages
    #' @param config optional configuration
    #' @param ctx optional tiledb context
    initialize = function(uri, dimension_name = NULL, verbose = TRUE, config = NULL, ctx = NULL) {
      super$initialize(uri, verbose, config, ctx)
      self$dimension_name <- dimension_name
      self
    },

    #' @description Set dimension values to slice from the array members.
    #' @param dims a named list of character vectors. Each must correspond to a
    #' dimension shared by all array members.
    #' @param attr_filter a TileDB query condition for attribute filtering
    #' pushdown.
    set_query = function(dims = NULL, attr_filter = NULL) {
      # stopifnot(
          # TODO: Utilize AnnotationGroup's dimension_name field for validation
        # "All 'dims' element names must match an array dimension" =
          # all(names(dims) %in% self$dimension_name)
      # )

      # capture unevaluated expression as a character vector
      attr_filter <- deparse(substitute(attr_filter))

      # Assuming all members are TileDBArrays
      for (member in names(self$members)) {
        self$members[[member]]$set_query(dims, attr_filter)
      }
    },

    #' @description Reset the group member queries.
    #' @param dims Clear the defined dimension ranges?
    #' @param attr_filter Clear the defined attribute filters?
    #' @return NULL
    reset_query = function(dims = TRUE, attr_filter = TRUE) {
      for (member in names(self$members)) {
        self$members[[member]]$reset_query(dims, attr_filter)
      }
    }
  )
)
