#' SOMA Collection
#'
#' @description
#' Class for representing a `SOMACollection`, which may contain of one or more
#' [`SOMA`]s.
#' @importFrom SeuratObject CreateSeuratObject Reductions Idents
#' @export
SOMACollection <- R6::R6Class(
  classname = "SOMACollection",
  inherit = TileDBGroup,

  public = list(
    #' @field uns Named list of unstructured objects.
    uns = list(),

    #' @description Create a new `SOMACollection`. The existing TileDB group is
    #'   opened at the specified array `uri` if one is present, otherwise a new
    #'   array group is created. The `members` field is populated with
    #'   `SOMA` objects for each URI passed explicitly to `soma_uris`, as
    #'   well `SOMA` objects discovered within the `SOMACollection` object's
    #'   TileDB group.
    #'
    #' @param uri URI of the TileDB group
    #' @param verbose Print status messages
    #' @param config optional configuration
    #' @param ctx optional tiledb context
    initialize = function(uri, verbose = TRUE, config = NULL, ctx = NULL) {
      super$initialize(uri, verbose, config, ctx)

      # For compatibility with SCDatasets created with <=0.1.2 we look for a
      # misc directory first and treat it as uns
      if ("misc" %in% names(self$members)) {
        warning("Found deprecated 'misc' directory in SOMACollection.")
        self$uns <- self$get_member("misc")
      } else {
        if ("uns" %in% names(self$members)) {
          self$uns <- self$get_member("uns")
        } else {
          self$uns <- TileDBGroup$new(
            uri = file_path(self$uri, "uns"),
            verbose = self$verbose
          )
          self$add_member(self$uns, name = "uns")
        }
      }

      # Special handling of Seurat commands array
      if ("commands" %in% names(self$uns$members)) {
        self$uns$members$commands <- CommandsArray$new(
          uri = self$uns$members$commands$uri,
          verbose = self$verbose
        )
      }

      self
    },

    #' @description Set query parameters to slice by dimension values or filter
    #' by attribute values that are applied to all SOMAs within the collection.
    #'
    #' @details
    #' See `SOMA$set_query()` for more information about querying mechanics.
    #'
    #' @param obs_ids,var_ids character vector containing observation- or
    #' variable-identifiers.
    #' @param obs_attr_filter,var_attr_filter a TileDB query condition for
    #' attribute filtering pushdown.
    set_query = function(
      obs_ids = NULL,
      var_ids = NULL,
      obs_attr_filter = NULL,
      var_attr_filter = NULL
    ) {
      stopifnot(
        "SOMACollection must contain a SOMA to query" =
          length(self$somas) > 0,
        "'obs_ids' must be a character vector" =
          is.null(obs_ids) || is.character(obs_ids),
        "'var_ids' must be a character vector" =
          is.null(var_ids) || is.character(var_ids)
      )

      # capture unevaluated expressions as character vectors
      obs_attr_filter <- deparse(substitute(obs_attr_filter))
      var_attr_filter <- deparse(substitute(var_attr_filter))

      for (soma in names(self$somas)) {
        self$members[[soma]]$set_query(
          obs_ids = obs_ids,
          var_ids = var_ids,
          obs_attr_filter = obs_attr_filter,
          var_attr_filter = var_attr_filter
        )
      }
    },

    #' @description Reset query dimensions and attribute filters.
    #' @return NULL
    reset_query = function() {
      stopifnot("No SOMAs to reset" = length(self$somas) > 0)
      for (member in names(self$somas)) {
        self$members[[member]]$reset_query()
      }
    },

    #' @description Convert a Seurat object to a TileDB-backed `SOMACollection`.
    #'
    #' ## Assays
    #' Each `[SeuratObject::Assay`] is converted to a [`SOMA`] and written to
    #' a nested TileDB group with a URI of `./soma_<assay>` where `<assay>`
    #' is the name of the Seurat assay.
    #'
    #' ## Identities
    #'
    #' Cell identities in the [`SeuratObject::Seurat`] are maintained by
    #' creating an `active_ident` attribute in `obs` that stores the factor
    #' levels as a character vector.
    #
    #' ## Dimensionality Reductions
    #'
    #' Dimensionality reduction results are stored as `obsm` and `varm` arrays
    #' within an `SOMA`. The [`SeuratObject::DimReduc`] object's `key` slot is
    #' used to determine which `SOMA` to store the results in. The array names
    #' are `(obsm|varm)_dimreduction_<name>`, where `<name>` is the name of the
    #' dimensionality reduction method (e.g., `"pca"`).
    #'
    #' @param object A [`SeuratObject::Seurat`] object.
    from_seurat = function(object) {
      stopifnot(inherits(object, "Seurat"))

      idents <- SeuratObject::Idents(object)
      if (nlevels(idents) > 1L) {
        object <- SeuratObject::AddMetaData(
          object = object,
          metadata = as.character(idents),
          col.name = "active_ident"
        )
      }

      assays <- SeuratObject::Assays(object)
      for (assay in assays) {
        if (is.null(self$members[[assay]])) {
          assay_uri <- file_path(self$uri, paste0("soma_", assay))
          soma <- SOMA$new(assay_uri, verbose = self$verbose, config = self$config, ctx = self$context)
          self$add_member(soma, name = assay)
        } else {
          soma <- self$members[[assay]]
        }
        assay_object <- object[[assay]]
        soma$from_seurat_assay(assay_object, obs = object[[]])
      }

      reductions <- SeuratObject::Reductions(object)
      if (!is_empty(reductions)) {
        for (reduction in reductions) {
          reduction_object <- SeuratObject::Reductions(object, slot = reduction)
          assay <- SeuratObject::DefaultAssay(reduction_object)
          self$members[[assay]]$add_seurat_dimreduction(
            object = reduction_object,
            technique = reduction
          )
        }
      }

      graphs <- SeuratObject::Graphs(object)
      if (!is_empty(graphs)) {
        for (graph in graphs) {
          graph_object <- SeuratObject::Graphs(object, slot = graph)
          assay <- SeuratObject::DefaultAssay(graph_object)
          technique <- sub(paste0(assay, "_"), "", graph, fixed = TRUE)
          self$members[[assay]]$obsp$add_seurat_graph(
            object = graph_object,
            technique = technique
          )
        }
      }

      commandNames <- SeuratObject::Command(object)
      if (!is_empty(commandNames)) {
        namedListOfCommands <- lapply(commandNames, SeuratObject::Command,  object=object)
        names(namedListOfCommands) <- commandNames

        commandsArray <- CommandsArray$new(
          uri = file_path(self$uns$uri, "commands"),
          verbose = self$verbose
        )
        commandsArray$from_named_list_of_commands(namedListOfCommands)
        if (is.null(self$uns$members$commands)) {
          self$uns$add_member(commandsArray, name = "commands")
        }
      }

      if (self$verbose) {
        msg <- sprintf("Finished converting Seurat object to %s", self$class())
        message(msg)
      }
    },

    #' @description Convert to a [SeuratObject::Seurat] object.
    #' @param project [`SeuratObject::Project`] name for the `Seurat` object
    #' @param batch_mode logical, if `TRUE`, batch query mode is enabled for
    #' retrieving `X` layers. See
    #' [`AssayMatrix$to_dataframe()`][`AssayMatrix`] for more information.
    to_seurat = function(project = "SeuratProject", batch_mode = FALSE) {
      stopifnot(is_scalar_character(project))

      assays <- lapply(
        X = self$somas,
        FUN = function(x) x$to_seurat_assay(batch_mode = batch_mode)
      )
      nassays <- length(assays)

      # cell-level obs metadata is stored in each soma, so for now we
      # just take the first soma's obs metadata
      obs_df <- self$somas[[1]]$obs$to_dataframe()

      # retain cell identities before restoring cell-level metadata
  idents <- obs_df$active_ident
      if (!is.null(idents)) {
        idents <- setNames(idents, rownames(obs_df))
        obs_df$active_ident <- NULL
      }

      object <- SeuratObject::CreateSeuratObject(
        counts = assays[[1]],
        project = project,
        meta.data = obs_df
      )

      if (!is.null(idents)) {
        SeuratObject::Idents(object) <- idents[SeuratObject::Cells(object)]
      }

      if (nassays > 1) {
        for (i in seq(2, nassays)) {
          assay <- names(assays)[i]
          object[[assay]] <- assays[[assay]]
        }
      }

      # dimreductions
      # Retrieve list of all techniques used in any soma's obsm/varm
      # dimensionality reduction arrays. The association between assay and
      # dimreduction is maintained by the DimReduc's `assay.used` slot.
      dimreductions <- lapply(
        self$somas,
        function(x) x$get_seurat_dimreductions_list()
      )
      object@reductions <- Reduce(base::c, dimreductions)

      # graphs
      graph_arrays <- lapply(self$somas,
        function(x) x$get_annotation_pairwise_matrix_arrays(prefix = "graph_")
      )
      if (!is_empty(graph_arrays)) {
        graph_arrays <- unlist(graph_arrays)
        graphs <- lapply(graph_arrays, function(x) x$to_seurat_graph())
        # TODO: Bit of a hack to recreate the graph names
        names(graphs) <- sub("\\.(obs|var)p\\.graph", "", names(graphs))
        object@graphs <- graphs
      }

      # command history
      if ("commands" %in% names(self$uns$members)) {
        commands_array <- self$uns$get_member("commands")
        object@commands <- commands_array$to_named_list_of_commands()
      }

      return(object)
    },

    #' @description List the [`SOMA`] URIs in the collection.
    #' @return A vector of URIs for each [`SOMA`] in the collection.
    soma_uris = function() {
      vapply_char(self$somas, function(x) x$uri)
    }
  ),

  active = list(
    #' @field somas Retrieve [`SOMA`] members.
    somas = function(value) {
      if (!missing(value)) {
        stop("somas is read-only, use 'add_member()' to add a new SOMA")
      }
      Filter(function(x) inherits(x, "SOMA"), self$members)
    }
  ),

  private = list(

    instantiate_members = function() {

      # with the exception of 'uns' all members should be SOMA objects
      # TODO: Use group metadata to indicate each member's class
      member_uris <- self$list_member_uris()

      # TODO: Remove misc check when SCDatasets/SCGroups/misc are defunct
      if ("misc" %in% names(member_uris)) {
        member_uris <- rename(member_uris, c(uns = "misc"))
      }

      uns_uri <- member_uris[names(member_uris) == "uns"]
      soma_uris <- member_uris[names(member_uris) != "uns"]
      names(soma_uris) <- sub("^(scgroup|soma)_", "", names(soma_uris))

      # TODO: Remove this switch when SCDatasets/SCGroups/misc are defunct
      if (self$class() == "SCDataset") {
        somas <- suppressWarnings(lapply(
          X = soma_uris,
          FUN = SCGroup$new,
          verbose = self$verbose,
          config = self$config,
          ctx = self$context
        ))
      } else {
        somas <- lapply(
          X = soma_uris,
          FUN = SOMA$new,
          verbose = self$verbose,
          config = self$config,
          ctx = self$context
        )
      }

      c(
        somas,
        lapply(uns_uri, TileDBGroup$new, verbose = self$verbose, config = self$config, ctx = self$context)
      )
    }
  )
)

#' Single-cell Dataset
#'
#' @description
#' Class for representing the now-deprecated SCDataset object, which has been
#' renamed to [`SOMACollection`].
#' @export
SCDataset <- R6::R6Class(
  classname = "SCDataset",
  inherit = SOMACollection,

  public = list(
    #' @description Create a new SCDataset object.
    #'
    #' @param uri URI of the TileDB group
    #' @param verbose Print status messages
    #' @param config optional configuration
    #' @param ctx optional tiledb context
    initialize = function(uri, verbose = TRUE, config = NULL, ctx = NULL) {
      .Deprecated(
        new = "SOMACollection",
        old = "SCDataset",
        package = "tiledbsc"
      )
      super$initialize(uri, verbose, config, ctx)
    },

    #' @description List the [`SOMA`] (formerly `SCGroup`) URIs in the
    #' collection.
    #' @return A vector of URIs for each [`SOMA`] in the collection.
    scgroup_uris = function() {
      .Deprecated(
        new = "soma_uris",
        old = "scgroup_uri",
        package = "tiledbsc"
      )
      self$soma_uris
    }
  ),

  active = list(
    #' @field scgroups Retrieve the [`SOMA`] (formerly `SCGroup`) members.
    scgroups = function(value) {
      if (!missing(value)) {
        stop("scgroups is read-only, use 'add_member()' to add a new SOMA")
      }
      .Deprecated(new = "somas", old = "scgroups", package = "tiledbsc")
      self$somas
    },

    #' @field misc An alias for `uns`.
    misc = function(value) {
      if (!missing(value)) stop("misc is read-only")
      .Deprecated(new = "uns", old = "misc", package = "tiledbsc")
      self$uns
    }
  )
)
