#' SOMA Experiment
#'
#' @description `SOMAExperiment` is a specialized [`SOMACollection`],
#' representing one or more modes of measurement across a single collection of
#' cells (aka a "multimodal dataset") with pre-defined fields: `obs` and `ms`
#' (see _Active Bindings_ below for details). (lifecycle: experimental)
#'
#' @export
SOMAExperiment <- R6::R6Class(
  classname = "SOMAExperiment",
  inherit = SOMACollectionBase,

  public = list(
    #' @description Loads the experiment as a \code{\link[SeuratObject]{Seurat}}
    #' object
    #'
    #' @param X_layers A named list of named character vectors describing the
    #' measurements to load and the layers within those measurements to read in;
    #' for example: \preformatted{
    #' list(
    #'   RNA = c(counts = "counts", data = "logcounts"),
    #'   ADT = c(counts = "counts")
    #' )
    #' }
    #' @template param-obs-index
    #' @param var_index A named character of column names in \code{var} for
    #' each measurement to use as feature names; for example: \preformatted{
    #' c(RNA = "gene_name", ADT = "protein_name")
    #' }
    #' Uses \code{paste0("feature", var_joinids())} by default
    #' @template param-obs-column-names
    #' @param var_column_names A named list of character vectors describing the
    #' columns in \code{var} for each measurement to load in as feature-level
    #' metadata; for example: \preformatted{
    #' list(
    #'   RNA = c("vst.mean", "vst.variable"),
    #'   ADT = c("ensembl_id")
    #' )
    #' }
    #' By default, loads in entire feature-level metadata for all measurements
    #' described in \code{X_layers}
    #' @param obsm_layers A named list of character vectors describing the
    #' arrays in \code{obsm} for each measurement to load in as
    #' dimensional reductions; for example: \preformatted{
    #' list(
    #'   RNA = c("pca", "umap"),
    #'   ADT = c("adtpca", "adtumap")
    #' )
    #' }
    #' By default, loads in all dimensional reductions for all measurements
    #' described in \code{X_layers}
    #' @param varm_layers A named list of named character vectors describing the
    #' arrays in \code{varm} to load in as feature loadings and which array in
    #' \code{obsm} they correspond to; for example: \preformatted{
    #' list(
    #'   RNA = c(pca = "PCs"),
    #'   ADT = c(adtpca = "ADTPCs")
    #' )
    #' }
    #' By default, will try to determine \code{varm_layers} from
    #' \code{obsm_layers} and load in all loadings for all dimensional
    #' reductions for all measurements described in \code{X_layers}
    #' @param obsp_layers A named list of character vectors describing the
    #' arrays in \code{obsp} for each measurement to load in as
    #' nearest-neighbor graphs; for example: \preformatted{
    #' list(
    #'   RNA = c("RNA_nn", "RNA_snn"),
    #'   ADT = c("ADT_nn")
    #' )
    #' }
    #' By default, loads in all nearest-neighbor graphs for all measurements for
    #' all measurements described in \code{X_layers}
    #'
    #' @return A \code{\link[SeuratObject]{Seurat}} object
    #'
    to_seurat = function(
      X_layers,
      obs_index = NULL,
      var_index = NULL,
      obs_column_names = NULL,
      var_column_names = NULL,
      obsm_layers = NULL,
      varm_layers = NULL,
      obsp_layers = NULL
    ) {
      .check_seurat_installed()
      stopifnot(
        "'X_layers' must be a named list" = is_named_list(
          X_layers,
          allow_empty = FALSE
        ),
        "'obs_index' must be a single character value" = is.null(obs_index) ||
          is_scalar_character(obs_index),
        "'var_index' must be a named character vector" = is_character_or_null(var_index),
        "'var_column_names' must be a named list" = is.null(var_column_names) ||
          is_named_list(var_column_names, allow_empty = FALSE),
        "'obsm_layers' must be a named list" = is.null(obsm_layers) ||
          is_scalar_logical(obsm_layers) ||
          is_named_list(obsm_layers, allow_empty = FALSE),
        "'varm_layers' must be a named list" = is.null(varm_layers) ||
          is_scalar_logical(varm_layers) ||
          is_named_list(varm_layers, allow_empty = FALSE),
        "'obsp_layers' must be a named list" = is.null(obsp_layers) ||
          is_scalar_logical(obsp_layers) ||
          is_named_list(obsp_layers, allow_empty = FALSE)
      )
      # Check `X_layers`
      if (!all(names(X_layers) %in% self$ms$names())) {
        msg <- paste(
          "The following measurements could not be found in this experiment:",
          string_collapse(setdiff(x = names(X_layers), y = self$ms$names()))
        )
        stop(paste(strwrap(msg), collapse = '\n'), call. = FALSE)
      }
      layer_check <- vapply_lgl(
        X = X_layers,
        FUN = function(x) {
          return(is.character(x) && is_named(x, allow_empty = FALSE))
        }
      )
      if (!all(layer_check)) {
        stop("All entries in 'X_layers' must be named lists", call. = FALSE)
      }
      layers <- names(X_layers)
      nlayers <- length(X_layers)
      null_list <- stats::setNames(
        object = vector(mode = 'list', length = nlayers),
        nm = layers
      )
      # Check `obs_index`
      if (is_scalar_character(obs_index)) {
        obs_index <- match.arg(obs_index, choices = self$obs$attrnames())
      }
      # Check `var_index`
      var_index <- var_index %||% null_list
      if (length(var_index) == 1L) {
        var_index <- stats::setNames(
          object = rep_len(x = var_index, length.out = nlayers),
          nm = layers
        )
      }
      stopifnot(
        "There must be one 'var_index' for every X layer" = length(var_index) == nlayers,
        "'var_index' must be named" = is_named(var_index, allow_empty = FALSE),
        "'var_index' must have the same names as 'X_layers'" = all(names(var_index) %in% layers)
      )
      # Check `var_column_names`
      var_column_names <- var_column_names %||% null_list
      stopifnot(
        "'var_column_names' must have the same names as 'X_layers'" = all(names(var_column_names) %in% layers)
      )
      # Check `obsm_layers`
      sublayer_check <- function(x, named = FALSE) {
        checks <- c(
          is.null(x),
          is_scalar_logical(x),
          if (isFALSE(named)) {
            is.character(x)
          } else {
            is.character(x) && is_named(x, allow_empty = FALSE)
          }
        )
        return(Reduce(f = `||`, x = checks))
      }
      obsm_layers <- obsm_layers %||% null_list
      if (is_scalar_logical(obsm_layers)) {
        obsm_layers <- stats::setNames(
          object = rep_len(x = obsm_layers, length.out = nlayers),
          nm = layers
        )
      }
      stopifnot(
        "'obsm_layers' must have the same names as 'X_layers'" = all(names(obsm_layers) %in% layers),
        "Every entry in 'obsm_layers' must be a character vector" = all(vapply_lgl(obsm_layers, sublayer_check))
      )
      # Check `varm_layers`
      varm_layers <- varm_layers %||% null_list
      if (is_scalar_logical(varm_layers)) {
        varm_layers <- stats::setNames(
          object = rep_len(x = varm_layers, length.out = nlayers),
          nm = layers
        )
      }
      stopifnot(
        "'varm_layers' must have the same names as 'X_layers'" = all(names(varm_layers) %in% layers),
        "Every entry in 'varm_layers' must be a named character vector" = all(vapply_lgl(
          X = varm_layers,
          FUN = sublayer_check,
          named = TRUE
        ))
      )
      # Check `obsp_layers`
      obsp_layers <- obsp_layers %||% null_list
      if (is_scalar_logical(obsp_layers)) {
        obsp_layers <- stats::setNames(
          object = rep_len(x = obsp_layers, length.out = nlayers),
          nm = layers
        )
      }
      stopifnot(
        "'obsp_layers' must have the same names as 'X_layers'" = all(names(obsp_layers) %in% layers),
        "Every entry in 'obsp_layers' must be a character vector" = all(vapply_lgl(obsp_layers, sublayer_check))
      )
      # Load in the first assay as the default assay
      active <- names(X_layers)[1L]
      query <- SOMAExperimentAxisQuery$new(
        experiment = self,
        measurement_name = active
      )
      object <- query$to_seurat(
        X_layers = X_layers[[active]],
        obs_index = obs_index,
        var_index = var_index[[active]],
        obs_column_names = obs_column_names,
        var_column_names = var_column_names[[active]],
        obsm_layers = obsm_layers[[active]],
        varm_layers =  varm_layers[[active]],
        obsp_layers = obsp_layers[[active]]
      )
      # Add alternate assays
      for (assay in setdiff(x = layers, y = active)) {
        query <- SOMAExperimentAxisQuery$new(
          experiment = self,
          measurement_name = assay
        )
        obj <- tryCatch(
          expr = query$to_seurat_assay(
            X_layers = X_layers[[assay]],
            obs_index = obs_index,
            var_index = var_index[[assay]],
            var_column_names = var_column_names[[assay]]
          ),
          error = function(e) {
            warning(conditionMessage(e), call. = FALSE, immediate. = TRUE)
            return(NULL)
          }
        )
        if (is.null(obj)) {
          next
        }
        object[[assay]] <- obj
        # Add reductions
        embeddings <- obsm_layers[[assay]]
        skip_reducs <- isFALSE(obsm_layers) || rlang::is_na(obsm_layers)
        obsm <- tryCatch(expr = self$ms$get(assay)$get('obsm'), error = null)
        if (is.null(obsm)) {
          if (!skip_reducs) {
            warning(
              'Dimensional reductions were requested for assay',
              sQuote(assay),
              ', but no reductions found',
              call. = FALSE,
              immediate. = TRUE
            )
          }
          skip_reducs <- TRUE
        }
        if (!skip_reducs) {
          if (isTRUE(embeddings)) {
            embeddings <- NULL
          }
          loadings <- varm_layers[['loadings']]
          if (isTRUE(loadings)) {
            loadings <- NULL
          }
          reductions <- .get_seurat_reductions(
            query = query,
            obsm_layers = embeddings,
            varm_layers = loadings,
            obs_index = obs_index,
            var_index = var_index[[assay]]
          )
          if (length(reductions)) {
            for (reduc in names(reductions)) {
              object[[reduc]] <- reductions
            }
          }
        }
        # Add graphs
        graphs <- obsp_layers[[assay]]
        obsp <- tryCatch(expr = self$ms$get(assay)$get('obsp'), error = null)
        if (is.null(obsp)) {
          if (!(isFALSE(graphs) || rlang::is_na(graphs))) {
            ''
          }
        }
      }
      return(object)
    }
  ),

  active = list(
    #' @field obs a [`SOMADataFrame`] containing primary annotations on the
    #' observation axis. The contents of the `soma_joinid` column define the
    #' observation index domain, `obs_id`. All observations for the
    #' `SOMAExperiment` must be defined in this dataframe.
    obs = function(value) {
      private$get_or_set_soma_field(value, "obs", "SOMADataFrame")
    },

    #' @field ms a [`SOMACollection`] of named [`SOMAMeasurement`]s.
    ms = function(value) {
      private$get_or_set_soma_field(value, "ms", "SOMACollection")
    }
  )
)
