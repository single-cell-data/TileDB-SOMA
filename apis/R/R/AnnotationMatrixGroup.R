#' Single-cell Annotation Matrix Group
#'
#' @description
#' Class for representing a TileDB group containing one or more
#' [`AnnotationMatrix`] arrays that share a common dimension name.
#' @export
AnnotationMatrixGroup <- R6::R6Class(
  classname = "AnnotationMatrixGroup",
  inherit = AnnotationGroup,

  public = list(

    #' @description Add a new [`AnnotationMatrix`] array to the group.
    #' @param data a [`matrix`] of annotation data to ingest. The `matrix` rows
    #' must be aligned to the [`SOMA`] dimension indicated by the group's
    #' `dimension_name`.
    #' @param name Name of the new variable annotation matrix.
    #' @param metadata Named list of metadata to add.
    add_annotation_matrix = function(data, name, metadata = NULL) {
      if (missing(name)) {
        stop("Must specify a `name` for the new AnnotationMatrix")
      }
      if (missing(data)) {
        stop("Must provide a `matrix` to ingest into the new AnnotationMatrix")
      }
      # TODO: Verify that the matrix is aligned to the group's dimension

      # create the new array
      array_uri <- file_path(self$uri, name)
      array <- AnnotationMatrix$new(
        uri = array_uri,
        verbose = self$verbose
      )

      array$from_matrix(data, self$dimension_name)
      if (!is.null(metadata)) array$add_metadata(metadata)
      if (is.null(self$members[[name]])) self$add_member(array, name)
    }
  ),

  private = list(
    instantiate_members = function() {
      lapply(self$list_member_uris(), AnnotationMatrix$new, verbose = self$verbose)
    }
  )
)
