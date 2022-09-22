#' Single-cell Annotation Matrix Group
#'
#' @description
#' Class for representing a TileDB group containing one or more
#' [`AnnotationPairwiseMatrix`] arrays that share a common dimension name.
#' @export
AnnotationPairwiseMatrixGroup <- R6::R6Class(
  classname = "AnnotationPairwiseMatrixGroup",
  inherit = AnnotationGroup,

  public = list(

    #' @description Add a new [`AnnotationPairwiseMatrix`] array to the group.
    #' @param data a [`matrix`] of annotation data to ingest. The `matrix` rows
    #' must be aligned to the dimension indicated by the group's
    #' `dimension_name`.
    #' @param name Name of the new pairwise annotation matrix.
    #' @param metadata Named list of metadata to add.
    add_matrix = function(data, name, metadata = NULL) {
      if (missing(name)) {
        stop("Must specify a `name` for the new AnnotationPairwiseMatrix")
      }
      if (missing(data)) {
        stop("Must provide a `matrix` to ingest into the new AnnotationPairwiseMatrix")
      }

      # TODO: Verify that the matrix is aligned to the group's dimension
      # create the new array
      array_uri <- file_path(self$uri, name)
      array <- AnnotationPairwiseMatrix$new(
        uri = array_uri,
        verbose = self$verbose
      )

      index_cols <- paste(self$dimension_name, c("i", "j"), sep = "_")
      array$from_matrix(data, index_cols)
      if (!is.null(metadata)) array$add_metadata(metadata)
      if(is.null(self$members[[name]])) self$add_member(array, name)
    },

    #' @description Convert a [`SeuratObject::Graph`] object to
    #' [`AnnotationPairwiseMatrix`].
    #'
    #' @details
    #' ## On-Disk Format
    #'
    #' Arrays are named `graph_<technique>`.
    #'
    #' ## Metadata
    #'
    #' - `assay_used`: Name of the assay used to generate the graph.
    #' - `graph_technique`: Name of the technique used to generate the graph.
    #' used.
    #'
    #' @param object A [`SeuratObject::Graph`] object.
    #' @param technique Name of the technique used to generate the graph
    #' (typically, `nn` or `snn`).
    add_seurat_graph = function(object, technique) {
      stopifnot(
        "Must provide a Seurat 'Graph' object" = inherits(object, "Graph"),
        "'technique' must be a scalar character" = is_scalar_character(technique)
      )

      prefix <- "graph_"
      assay <- SeuratObject::DefaultAssay(object)
      array_name <- paste0(prefix, technique)

      self$add_matrix(
        data = as(object, "dgTMatrix"),
        name = array_name,
        metadata = list(
          assay_used = assay,
          graph_technique = technique
        )
      )
    }
  ),

  private = list(
    instantiate_members = function() {
      lapply(self$list_member_uris(), AnnotationPairwiseMatrix$new, verbose = self$verbose)
    }
  )
)
