#' @importFrom rlang is_na
#' @importFrom methods as new validObject
#'
NULL

#' `SOMAExperiment` Axis Query
#' @description Perform an axis-based query against a [`SOMAExperiment`].
#'
#' `SOMAExperimentAxisQuery` allows easy selection and extraction of data from a
#' single [`SOMAMeasurement`] in a [`SOMAExperiment`], by `obs`/`var` (axis)
#' coordinates and/or value filter. The primary use for this class is slicing
#' [`SOMAExperiment`] `X` layers by `obs` or `var` value and/or coordinates.
#' (lifecycle: experimental)
#'
#' ## X Layer Support
#'
#' Slicing on [`SOMASparseNDArray`] `X` matrices is supported;
#' slicing on [`SOMADenseNDArray`] is not supported at this time.
#'
#' ## Result Size
#' `SOMAExperimentAxisQuery` query class assumes it can store the full result of
#' both axis dataframe queries in memory, and only provides incremental access
#' to the underlying X NDArray. Accessors such as `n_obs` and `n_vars` codify
#' this in the class.
#'
#' @importFrom arrow concat_arrays
#' @export
SOMAExperimentAxisQuery <- R6::R6Class(
  classname = "SOMAExperimentAxisQuery",

  public = list(
    #' @description Create a new `SOMAExperimentAxisQuery` object.
    #' @param experiment A [`SOMAExperiment`] object.
    #' @param measurement_name The name of the measurement to query.
    #' @param obs_query,var_query An [`SOMAAxisQuery`] object for the obs/var
    #' axis.
    initialize = function(
      experiment,
      measurement_name,
      obs_query = NULL,
      var_query = NULL
    ) {

      stopifnot(
        "experiment must be a SOMAExperiment" =
          inherits(experiment, "SOMAExperiment"),
        "SOMAExperiment does not exist" = experiment$exists(),
        "Must specify a single measurement to query" =
          is_scalar_character(measurement_name),
        "Measurement does not exist in the experiment" =
          measurement_name %in% experiment$ms$names(),
        is.null(obs_query) || inherits(obs_query, "SOMAAxisQuery"),
        is.null(var_query) || inherits(var_query, "SOMAAxisQuery")
      )

      private$.experiment <- experiment
      private$.measurement_name <- measurement_name
      private$.obs_query <- obs_query %||% SOMAAxisQuery$new()
      private$.var_query <- var_query %||% SOMAAxisQuery$new()
      private$.joinids <- JoinIDCache$new(self)
      private$.indexer <- SOMAAxisIndexer$new(self)
    },

    #' @description Retrieve obs [`arrow::Table`]
    #' @param column_names A character vector of column names to retrieve
    obs = function(column_names = NULL) {
      obs_query <- self$obs_query
      self$obs_df$read(
        coords = recursively_make_integer64(obs_query$coords),
        value_filter = obs_query$value_filter,
        column_names = column_names
      )$concat()
    },

    #' @description Retrieve var [`arrow::Table`]
    #' @param column_names A character vector of column names to retrieve
    var = function(column_names = NULL) {
      var_query <- self$var_query
      self$var_df$read(
        coords = recursively_make_integer64(var_query$coords),
        value_filter = var_query$value_filter,
        column_names = column_names
      )$concat()
    },

    #' @description Retrieve `soma_joinids` as an [`arrow::Array`] for `obs`.
    obs_joinids = function() {
      arrow::concat_arrays(private$.joinids$obs())
    },

    #' @description Retrieve `soma_joinids` as an [`arrow::Array`] for `var`.
    var_joinids = function() {
      arrow::concat_arrays(private$.joinids$var())
    },

    #' @description Retrieves an `X` layer as a link{SOMASparseNDArrayRead}
    #' @param layer_name The name of the layer to retrieve.
    X = function(layer_name) {
      stopifnot(
        "Must specify an X layer name" = !missing(layer_name),
        "Must specify a single X layer name" = is_scalar_character(layer_name),
        assert_subset(layer_name, self$ms$X$names(), "layer")
      )

      x_layer <- self$ms$X$get(layer_name)
      stopifnot(
        "X layer must be a SOMASparseNDArray" =
          inherits(x_layer, "SOMASparseNDArray")
      )

      # TODO: Stop converting to vectors when SOMAArrayReader supports arrow arrays
      x_layer$read(coords = list(
        self$obs_joinids()$as_vector(),
        self$var_joinids()$as_vector()
      ))$tables()$concat()
    },

    #' @description Reads the entire query result as a list of
    #' [`arrow::Table`]s. This is a low-level routine intended to be used by
    #' loaders for other in-core formats, such as `Seurat`, which can be created
    #' from the resulting Tables.
    #'
    #' @param X_layers The name(s) of the `X` layer(s) to read and return.
    #' @param obs_column_names,var_column_names Specify which column names in
    #' `var` and `obs` dataframes to read and return.
    read = function(
      X_layers = NULL, obs_column_names = NULL, var_column_names = NULL) {
      stopifnot(
        "'X_layers' must be a character vector" =
          is.null(X_layers) || is.character(X_layers),
        assert_subset(X_layers, self$ms$X$names(), "layer"),
        "'obs_column_names' must be a character vector" =
          is.null(obs_column_names) || is.character(obs_column_names),
        assert_subset(obs_column_names, self$obs_df$colnames(), "column"),
        "'var_column_names' must be a character vector" =
          is.null(var_column_names) || is.character(var_column_names),
        assert_subset(var_column_names, self$var_df$colnames(), "column")
      )

      x_collection <- self$ms$X
      X_layers <- X_layers %||% x_collection$names()

      # Named list of SOMASparseNDArrays
      x_arrays <- Map(
        f = function(layer_name) {
          x_layer <- x_collection$get(layer_name)
          stopifnot(
            "X layer must be a SOMASparseNDArray" =
              inherits(x_layer, "SOMASparseNDArray")
          )
          x_layer
        },
        layer_name = X_layers
      )

      # TODO: parallelize with futures
      obs_ft <- self$obs(obs_column_names)
      var_ft <- self$var(var_column_names)

      x_matrices <- lapply(x_arrays, function(x_array) {
          x_array$read(coords = list(
            self$obs_joinids()$as_vector(),
            self$var_joinids()$as_vector()
          ))$tables()$concat()
        }
      )

      SOMAAxisQueryResult$new(
        obs = obs_ft, var = var_ft, X_layers = x_matrices
      )
    },
    #' @description Loads the query as a \code{\link[SeuratObject]{Seurat}} object
    #'
    #' @template param-x-layers-v3
    #' @template param-obs-index
    #' @template param-var-index
    #' @param obs_column_names Names of columns in \code{obs} to add as
    #' cell-level meta data; by default, loads all columns
    #' @template param-var-column-names
    #' @param obsm_layers Names of arrays in \code{obsm} to load in as the
    #' cell embeddings; pass \code{FALSE} to suppress loading in any
    #' dimensional reductions; by default, loads all dimensional
    #' reduction information
    #' @param varm_layers Named vector of arrays in \code{varm} to load in as
    #' the feature loadings; names must be names of array in \code{obsm} (eg.
    #' \code{varm_layers = c(X_pca = 'PCs')}); will try to determine
    #' \code{varm_layers} from \code{obsm_layers}
    #' @param obsp_layers Names of arrays in \code{obsp} to load in as
    #' \code{\link[SeuratObject]{Graph}s}; by default, loads all graphs
    #'
    #' @return A \code{\link[SeuratObject]{Seurat}} object
    #'
    to_seurat = function(
      X_layers = c(counts = 'counts', data = 'logcounts'),
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
        "'obs_index' must be a single character value" = is.null(obs_index) ||
          (is_scalar_character(obs_index) && !is.na(obs_index)),
        "'obs_column_names' must be a character vector" = is.null(obs_column_names) ||
          is.character(obs_column_names) ||
          is_scalar_logical(obs_column_names),
        "'obsm_layers' must be a character vector" = is.null(obsm_layers) ||
          is.character(obsm_layers) ||
          is_scalar_logical(obsm_layers),
        "'varm_layers' must be a named character vector" = is.null(varm_layers) ||
          (is.character(varm_layers) && is_named(varm_layers, allow_empty = FALSE)) ||
          is_scalar_logical(varm_layers),
        "'obsp_layers' must be a character vector" = is.null(obsp_layers) ||
          is.character(obsp_layers) ||
          is_scalar_logical(obsp_layers)
      )
      tryCatch(
        expr = self$obs_df,
        error = function(...) {
          stop("No 'obs' found", call. = FALSE)
        }
      )
      # Load in the cells
      cells <- if (is.null(obs_index)) {
        paste0('cell', self$obs_joinids()$as_vector())
      } else {
        obs_index <- match.arg(
          arg = obs_index,
          choices = self$obs_df$attrnames()
        )
        self$obs(obs_index)$GetColumnByName(obs_index)$as_vector()
      }
      # Load in the assay
      assay <- self$to_seurat_assay(
        X_layers = X_layers,
        obs_index = obs_index,
        var_index = var_index,
        var_column_names = var_column_names
      )
      object <- SeuratObject::CreateSeuratObject(
        counts = assay,
        assay = private$.measurement_name
      )
      # Load in cell-level meta data
      if (isTRUE(obs_column_names)) {
        obs_column_names <- NULL
      }
      obs_column_names <- obs_column_names %||% setdiff(
        x = self$obs_df$attrnames(),
        y = obs_index
      )
      if (!(isFALSE(obs_column_names) || rlang::is_na(obs_column_names))) {
        obs <- as.data.frame(
          x = self$obs(obs_column_names)$to_data_frame()
        )
        row.names(obs) <- cells
        object[[names(obs)]] <- obs
      }
      # Load in reductions
      ms_embed <- tryCatch(expr = self$ms$obsm$names(), error = null)
      skip_reducs <- isFALSE(obsm_layers) || rlang::is_na(obsm_layers)
      if (is.null(ms_embed)) {
        if (!(skip_reducs || is.null(obsm_layers))) {
          warning("No reductions found", call. = FALSE, immediate. = TRUE)
        }
        skip_reducs <- TRUE
      }
      if (!skip_reducs) {
        names(ms_embed) <- .anndata_to_seurat_reduc(ms_embed)
        if (isTRUE(obsm_layers)) {
          obsm_layers <- NULL
        }
        obsm_layers <- obsm_layers %||% ms_embed
        # Match loadings to embeddings
        ms_load <- tryCatch(expr = self$ms$varm$names(), error = null)
        if (isTRUE(varm_layers)) {
          varm_layers <- NULL
        } else if (rlang::is_na(varm_layers)) {
          varm_layers <- FALSE
        }
        if (is.null(ms_load)) {
          if (!(isFALSE(varm_layers) || is.null(varm_layers))) {
            warning("No loadings found", call. = FALSE, immediate. = TRUE)
          }
          varm_layers <- FALSE
        }
        if (!isFALSE(varm_layers)) {
          names(ms_load) <- ms_embed[.anndata_to_seurat_reduc(ms_load, 'loadings')]
          varm_layers <- varm_layers %||% ms_load
          reduc_misisng <- setdiff(x = names(varm_layers), y = names(ms_load))
          if (length(reduc_misisng) == length(varm_layers)) {
            warning(
              "None of the reductions specified in 'varm_layers' can be found",
              call. = FALSE,
              immediate. = TRUE
            )
            varm_layers <- FALSE
          } else if (length(reduc_misisng)) {
            warning(
              paste(
                strwrap(paste(
                  "The reductions for the following loadings cannot be found in 'varm':",
                  sQuote(varm_layers[reduc_misisng]),
                  collapse = ', '
                )),
                collapse = '\n'
              ),
              call. = FALSE,
              immediate. = TRUE
            )
            varm_layers <- varm_layers[!names(varm_layers) %in% reduc_misisng]
          }
        }
        # Read in reductions and add to `object`
        for (embed in obsm_layers) {
          if (embed %in% names(ms_embed)) {
            embed <- ms_embed[embed]
          }
          rname <- .anndata_to_seurat_reduc(embed)
          reduc <- withCallingHandlers(
            expr = tryCatch(
              expr = self$to_seurat_reduction(
                obsm_layer = embed,
                varm_layer = ifelse(
                  embed %in% names(varm_layers),
                  yes = varm_layers[embed],
                  no = FALSE
                ),
                obs_index = obs_index,
                var_index = var_index
              ),
              error = err_to_warn
            ),
            noArrayWarning = function(w) {
              invokeRestart("muffleWarning")
            }
          )
          if (!inherits(reduc, 'DimReduc')) {
            next
          }
          object[[rname]] <- reduc
        }
      }
      # Load in graphs
      ms_graphs <- tryCatch(expr = self$ms$obsp$names(), error = null)
      skip_graphs <- isFALSE(obsp_layers) || rlang::is_na(obsp_layers)
      if (is.null(ms_graphs)) {
        if (!(skip_graphs || is.null(obsp_layers))) {
          warning("No graphs found in 'obsp'", call. = FALSE, immediate. = TRUE)
        }
        skip_graphs <- TRUE
      }
      if (!skip_graphs) {
        if (isTRUE(obsp_layers)) {
          obsp_layers <- NULL
        }
        obsp_layers <- obsp_layers %||% ms_graphs
        for (grph in obsp_layers) {
          mat <- withCallingHandlers(
            expr = tryCatch(
              expr = self$to_seurat_graph(obsp_layer = grph, obs_index = obs_index),
              error = err_to_warn
            ),
            noArrayWarning = function(w) {
              invokeRestart("muffleWarning")
            }
          )
          if (!inherits(mat, 'Graph')) {
            next
          }
          object[[grph]] <- mat
        }
      }
      # Validate and return
      validObject(object)
      return(object)
    },
    #' @description Loads the query as a Seurat \code{\link[SeuratObject]{Assay}}
    #'
    #' @return An \code{\link[SeuratObject]{Assay}} object
    #'
    to_seurat_assay = function(
      X_layers = c(counts = 'counts', data = 'logcounts'),
      obs_index = NULL,
      var_index = NULL,
      var_column_names = NULL
    ) {
      version <- 'v3'
      .check_seurat_installed()
      stopifnot(
        "'X_layers' must be a named character vector" = is.character(X_layers) &&
          is_named(X_layers, allow_empty = FALSE),
        "'version' must be a single character value" = is_scalar_character(version),
        "'obs_index' must be a single character value" = is.null(obs_index) ||
          (is_scalar_character(obs_index) && !is.na(obs_index)),
        "'var_index' must be a single character value" = is.null(var_index) ||
          (is_scalar_character(var_index) && !is.na(var_index)),
        "'var_column_names' must be a character vector" = is.null(var_column_names) ||
          is.character(var_column_names) ||
          is_scalar_logical(var_column_names)
      )
      match.arg(version, choices = 'v3')
      features <- if (is.null(var_index)) {
        paste0('feature', self$var_joinids()$as_vector())
      } else {
        var_index <- match.arg(
          arg = var_index,
          choices = self$var_df$attrnames()
        )
        self$var(var_index)$GetColumnByName(var_index)$as_vector()
      }
      cells <- if (is.null(obs_index)) {
        paste0('cell', self$obs_joinids()$as_vector())
      } else {
        obs_index <- match.arg(
          arg = obs_index,
          choices = self$obs_df$attrnames()
        )
        self$obs(obs_index)$GetColumnByName(obs_index)$as_vector()
      }
      # Check the layers
      assert_subset(x = X_layers, y = self$ms$X$names(), type = 'X_layer')
      # Read in the assay
      obj <- switch(
        EXPR = version,
        v3 = {
          assert_subset(
            x = names(X_layers),
            y = c('counts', 'data', 'scale.data'),
            type = 'Seurat slot'
          )
          private$.to_seurat_assay_v3(
            counts = tryCatch(expr = X_layers[['counts']], error = null),
            data = tryCatch(expr = X_layers[['data']], error = null),
            scale_data = tryCatch(expr = X_layers[['scale.data']], error = null),
            cells = cells,
            features = features
          )
        }
      )
      # Set the key
      SeuratObject::Key(obj) <- SeuratObject::Key(
        object = tolower(private$.measurement_name),
        quiet = TRUE
      )
      # Add feature-level meta data
      if (isTRUE(var_column_names)) {
        var_column_names <- NULL
      }
      var_column_names <- var_column_names %||% setdiff(
        x = self$var_df$attrnames(),
        y = var_index
      )
      if (!(isFALSE(var_column_names) || rlang::is_na(var_column_names))) {
        var <- as.data.frame(self$var(var_column_names)$to_data_frame())
        row.names(var) <- features
        obj[[names(var)]] <- var
      }
      validObject(obj)
      return(obj)
    },
    #' @description Loads the query as a Seurat
    #' \link[SeuratObject:DimReduc]{dimensional reduction}
    #'
    #' @param obsm_layer Name of array in \code{obsm} to load as the
    #' cell embeddings
    #' @param varm_layer Name of the array in \code{varm} to load as the
    #' feature loadings; by default, will try to determine \code{varm_layer}
    #' from \code{obsm_layer}
    #'
    #' @return A \code{\link[SeuratObject]{DimReduc}} object
    #'
    to_seurat_reduction = function(
      obsm_layer,
      varm_layer = NULL,
      obs_index = NULL,
      var_index = NULL
    ) {
      .check_seurat_installed()
      stopifnot(
        "'obsm_layer' must be a single character value" = is_scalar_character(obsm_layer),
        "'varm_layer' must be a single character value" = is.null(varm_layer) ||
          is_scalar_character(varm_layer) ||
          is_scalar_logical(varm_layer),
        "one of 'obsm_layer' or 'varm_layer' must be provided" =
          (is_scalar_character(obsm_layer) || is_scalar_logical(obsm_layer)) ||
          (is_scalar_character(varm_layer) || is_scalar_logical(varm_layer)),
        "'obs_index' must be a single character value" = is.null(obs_index) ||
          (is_scalar_character(obs_index) && !is.na(obs_index)),
        "'var_index' must be a single character value" = is.null(var_index) ||
          (is_scalar_character(var_index) && !is.na(var_index))
      )
      # Check embeddings/loadings
      ms_embed <- tryCatch(expr = self$ms$obsm$names(), error = null)
      ms_load <- tryCatch(expr = self$ms$varm$names(), error = null)
      if (is.null(ms_embed) && is.null(ms_load)) {
        warning(warningCondition(
          "No reductions present",
          class = c("noObsmWarning", "noArrayWarning"),
          call = NULL
        ))
        return(NULL)
      }
      if (is.null(ms_embed)) {
        stop("No embeddings in obsm present", call. = FALSE)
      }
      names(ms_embed) <- .anndata_to_seurat_reduc(ms_embed)
      if (is.null(ms_load) && !is.null(varm_layer)) {
        warning(warningCondition(
          "No loadings present in 'varm'",
          class = c("noVarmWarning", "noArrayWarning"),
          call = NULL
        ))
        varm_layer <- NULL
      } else {
        names(ms_load) <- .anndata_to_seurat_reduc(ms_load, 'loadings')
      }
      # Check provided names
      assert_subset(
        x = obsm_layer,
        y = c(ms_embed, names(ms_embed)),
        type = 'cell embedding'
      )
      if (is_scalar_character(varm_layer)) {
        assert_subset(
          x = varm_layer,
          y = c(ms_load, names(ms_load)),
          'feature loading'
        )
      }
      # Find Seurat name
      seurat <- c(
        embeddings = unname(.anndata_to_seurat_reduc(obsm_layer)),
        loadings = tryCatch(
          expr = unname(.anndata_to_seurat_reduc(varm_layer, 'loadings')),
          error = null
        )
      )
      if (length(seurat) == 2L && !identical(seurat[['embeddings']], y = seurat[['loadings']])) {
        stop(
          paste(
            strwrap(paste0(
              "The embeddings requested (",
              sQuote(obsm_layer),
              ") do not match the loadings requested (",
              sQuote(varm_layer),
              "); using the embeddings to create a Seurat name (",
              sQuote(seurat[['embeddings']]),
              ")"
            )),
            collapse = '\n'
          ),
          call. = FALSE,
          immediate. = TRUE
        )
      }
      seurat <- seurat[['embeddings']]
      # Create a Seurat key
      key <- SeuratObject::Key(
        object = switch(
          EXPR = seurat,
          pca = 'PC',
          ica = 'IC',
          tsne = 'tSNE',
          toupper(seurat)
        ),
        quiet = TRUE
      )
      # Read in cell embeddings
      # Translate Seurat name to AnnData name
      if (obsm_layer %in% names(ms_embed)) {
        obsm_layer <- ms_embed[obsm_layer]
      }
      # Get cell names
      cells <- if (is.null(obs_index)) {
        paste0('cell', self$obs_joinids()$as_vector())
      } else {
        obs_index <- match.arg(
          arg = obs_index,
          choices = self$obs_df$attrnames()
        )
        self$obs(obs_index)$GetColumnByName(obs_index)$as_vector()
      }
      embed <- self$ms$obsm$get(obsm_layer)
      coords <- list(
        cells = self$obs_joinids()$as_vector(),
        dims = seq_len(as.integer(embed$shape()[2L])) - 1L
      )
      embed_mat <- if (inherits(embed, 'SOMASparseNDArray')) {
        this_mat <- embed$read()$sparse_matrix(zero_based=TRUE)$concat()
        this_mat <- this_mat$take(coords$cells, coords$dims)
        this_mat <- this_mat$get_one_based_matrix()
        this_mat <- as(this_mat, "CsparseMatrix")
        as.matrix(this_mat)
      } else if (inherits(embed, 'SOMADenseNDArray')) {
        warning(
          paste(
            strwrap(paste(
              "Embeddings for",
              sQuote(obsm_layer),
              "are encoded as dense instead of sparse; all arrays should be saved as",
              sQuote('SOMASparseNDArrays')
            )),
            collapse = '\n'
          ),
          call. = FALSE
        )
        embed$read_dense_matrix(unname(coords))
      } else {
        stop("Unknown SOMA Array type: ", class(embed)[1L], call. = FALSE)
      }
      # Set matrix names
      rownames(embed_mat) <- cells
      colnames(embed_mat) <- paste0(key, seq_len(ncol(embed_mat)))
      # Autoset loadings if needed
      if (is.null(varm_layer) || isTRUE(varm_layer)) {
        if (seurat %in% c(names(ms_load), ms_load)) {
          varm_layer <- seurat
        }
      }
      # Read in feature loadings
      if (is_scalar_character(varm_layer)) {
        # Translate Seurat name to AnnData name
        if (varm_layer %in% names(ms_load)) {
          varm_layer <- ms_load[varm_layer]
        }
        # Get feature names
        features <- if (is.null(var_index)) {
          paste0('feature', self$var_joinids()$as_vector())
        } else {
          var_index <- match.arg(
            arg = var_index,
            choices = self$var_df$attrnames()
          )
          self$var(var_index)$GetColumnByName(var_index)$as_vector()
        }
        loads <- self$ms$varm$get(varm_layer)
        coords <- list(
          features = self$var_joinids()$as_vector(),
          dims = seq_len(as.integer(loads$shape()[2L])) - 1L
        )
        load_mat <- if (inherits(loads, 'SOMASparseNDArray')) {
          this_mat <- loads$read()$sparse_matrix(zero_based=TRUE)$concat()
          this_mat <- embed$read_sparse_matrix_zero_based()
          this_mat <- this_mat$take(coords$features, coords$dims)
          this_mat <- this_mat$get_one_based_matrix()
          this_mat <- as(this_mat, "CsparseMatrix")
          as.matrix(this_mat)
        } else if (inherits(loads, 'SOMADenseNDArray')) {
          warning(
            paste(
              strwrap(paste(
                "Loadings for",
                sQuote(varm_layer),
                "are encoded as dense instead of sparse; all arrays should be saved as",
                sQuote('SOMASparseNDArrays')
              )),
              collapse = '\n'
            ),
            call. = FALSE,
            immediate. = TRUE
          )
          loads$read_dense_matrix(unname(coords))
        } else {
          stop("Unknown SOMA Array type: ", class(loads)[1L], call. = FALSE)
        }
        # Set matrix names
        rownames(load_mat) <- features
        colnames(load_mat) <- paste0(key, seq_len(ncol(load_mat)))
        if (!is.null(embed_mat) && ncol(load_mat) != ncol(embed_mat)) {
          stop("The loadings do not match the embeddings", call. = FALSE)
        }
      } else {
        load_mat <- NULL
      }
      # Create the DimReduc
      return(SeuratObject::CreateDimReducObject(
        embeddings = embed_mat,
        loadings = load_mat %||% methods::new('matrix'),
        assay = private$.measurement_name,
        global = !seurat %in% c('pca', 'ica'),
        key = key
      ))
    },
    #' @description Loads the query as a Seurat \link[SeuratObject:Graph]{graph}
    #'
    #' @param obsp_layer Name of array in \code{obsp} to load as the graph
    #'
    #' @return A \code{\link[SeuratObject]{Graph}} object
    #'
    to_seurat_graph = function(obsp_layer, obs_index = NULL) {
      .check_seurat_installed()
      stopifnot(
        "'obsp_layer' must be a single character value" = is_scalar_character(obsp_layer),
        "'obs_index' must be a single character value" = is.null(obs_index) ||
          (is_scalar_character(obs_index) && !is.na(obs_index))
      )
      # Check graph name
      ms_graph <- tryCatch(expr = self$ms$obsp$names(), error = null)
      if (is.null(ms_graph)) {
        warning(
          warningCondition(
            "No graphs present",
            class = c("noObspWarning", "noArrayWarning")
          ),
          call. = FALSE,
          immediate. = TRUE
        )
        return(NULL)
      }
      # Check provided graph name
      obsp_layer <- match.arg(arg = obsp_layer, choices = ms_graph)
      mat <- self$ms$obsp$get(obsp_layer)$read()$sparse_matrix(zero_based=FALSE)$concat()$get_one_based_matrix()
      mat <- as(mat, "CsparseMatrix")
      idx <- self$obs_joinids()$as_vector() + 1L
      mat <- mat[idx, idx]
      mat <- as(mat, 'Graph')
      cells <- if (is.null(obs_index)) {
        paste0('cell', self$obs_joinids()$as_vector())
      } else {
        obs_index <- match.arg(
          arg = obs_index,
          choices = self$obs_df$attrnames()
        )
        self$obs(obs_index)$GetColumnByName(obs_index)$as_vector()
      }
      dimnames(mat) <- list(cells, cells)
      SeuratObject::DefaultAssay(mat) <- private$.measurement_name
      validObject(mat)
      return(mat)
    }
  ),

  active = list(

    #' @field experiment The parent [`SOMAExperiment`] object.
    experiment = function(value) {
      if (!missing(value)) read_only_error("experiment")
      private$.experiment
    },

    #' @field indexer The [`SOMAAxisIndexer`] object.
    indexer = function(value) {
      if (!missing(value)) read_only_error("indexer")
      private$.indexer
    },

    #' @field obs_query The `obs` [`SOMAAxisQuery`] object.
    obs_query = function(value) {
      if (!missing(value)) read_only_error("obs_query")
      private$.obs_query
    },

    #' @field var_query The `var` [`SOMAAxisQuery`] object.
    var_query = function(value) {
      if (!missing(value)) read_only_error("var_query")
      private$.var_query
    },

    #' @field n_obs The number of `obs` axis query results.
    n_obs = function(value) {
      if (!missing(value)) read_only_error("n_obs")
      length(self$obs_joinids())
    },

    #' @field n_vars The number of `var` axis query results.
    n_vars = function(value) {
      if (!missing(value)) read_only_error("n_vars")
      length(self$var_joinids())
    },

    #' @field obs_df The `obs` [`SOMADataFrame`] object.
    obs_df = function(value) {
      if (!missing(value)) read_only_error("obs_df")
      self$experiment$obs
    },

    #' @field var_df The `var` [`SOMADataFrame`] object for the specified
    #' `measurement_name`.
    var_df = function(value) {
      if (!missing(value)) read_only_error("var_df")
      self$ms$var
    },

    #' @field ms The [`SOMAMeasurement`] object for the specified
    #' `measurement_name`.
    ms = function(value) {
      if (!missing(value)) read_only_error("ms")
      self$experiment$ms$get(private$.measurement_name)
    }
  ),

  private = list(
    .experiment = NULL,
    .measurement_name = NULL,
    .obs_query = NULL,
    .var_query = NULL,
    .joinids = NULL,
    .indexer = NULL,
    .as_matrix = function(table, repr = 'C', transpose = FALSE) {
      stopifnot(
        "'table' must be an Arrow table" = inherits(table, 'Table'),
        "'repr' must be a single character value" = is_scalar_character(repr),
        "'transpose' must be a single logical value" = is_scalar_logical(transpose),
        "'table' must have column names 'soma_dim_0', 'soma_dim_1', and 'soma_data'" =
          all(c('soma_dim_0', 'soma_dim_1', 'soma_data') %in% table$ColumnNames())
      )
      repr <- match.arg(arg = repr, choices = c('C', 'R', 'T', 'D'))
      obs <- table$GetColumnByName('soma_dim_0')$as_vector()
      var <- table$GetColumnByName('soma_dim_1')$as_vector()
      mat <- Matrix::sparseMatrix(
        i = self$indexer$by_obs(obs)$as_vector() + 1L,
        j = self$indexer$by_var(var)$as_vector() + 1L,
        x = table$GetColumnByName('soma_data')$as_vector(),
        dims = c(self$n_obs, self$n_vars),
        repr = switch(EXPR = repr, D = 'T', repr)
      )
      if (isTRUE(transpose)) {
        mat <- Matrix::t(mat)
      }
      if (repr == 'D') {
        mat <- as.matrix(mat)
      }
      return(mat)
    },
    .to_seurat_assay_v3 = function(
      counts,
      data,
      scale_data = NULL,
      cells = NULL,
      features = NULL
    ) {
      .check_seurat_installed()
      stopifnot(
        "'data' must be a single character value" = is.null(data) ||
          is_scalar_character(data),
        "'counts' must be a single character value" = is.null(counts) ||
          is_scalar_character(counts),
        "one of 'counts' or 'data' must be provided" = is_scalar_character(counts) ||
          is_scalar_character(data),
        "'scale_data' must be a single character value" = is.null(scale_data) ||
          is_scalar_character(scale_data),
        "'cells' must be a character vector" = is.character(cells),
        "'features' must be a character vector" = is.character(features)
      )
      if (!length(x = cells) == self$n_obs) {
        stop("'cells' must have a length of ", self$n_obs, call. = FALSE)
      }
      if (!length(x = features) == self$n_vars) {
        stop("'features' must have a length of ", self$n_vars, call. = FALSE)
      }
      dnames <- list(features, cells)
      # Read in `data` slot
      if (is_scalar_character(data)) {
        dmat <- private$.as_matrix(
          table = self$X(data),
          repr = 'C',
          transpose = TRUE
        )
        dimnames(dmat) <- dnames
        obj <- SeuratObject::CreateAssayObject(data = dmat)
      }
      # Add the `counts` slot
      if (is_scalar_character(counts)) {
        cmat <- private$.as_matrix(
          table = self$X(counts),
          repr = 'C',
          transpose = TRUE
        )
        dimnames(cmat) <- dnames
        obj <- if (is_scalar_character(data)) {
          SeuratObject::SetAssayData(
            object = obj,
            slot = 'counts',
            new.data = cmat
          )
        } else {
          SeuratObject::CreateAssayObject(counts = cmat)
        }
      }
      # Add the `scale.data` slot
      if (is_scalar_character(scale_data)) {
        smat <- private$.as_matrix(
          table = self$X(scale_data),
          repr = 'D',
          transpose = TRUE
        )
        dimnames(smat) <- dnames
        obj <- SeuratObject::SetAssayData(
          object = obj,
          slot = 'scale.data',
          new.data = smat
        )
      }
      # Return the assay
      validObject(obj)
      return(obj)
    }
  )
)

JoinIDCache <- R6::R6Class(
  classname = "JoinIDCache",
  public = list(
    #' @field query The [`SOMAExperimentAxisQuery`] object to build indices for.
    query = NULL,

    initialize = function(query) {
      stopifnot(inherits(query, "SOMAExperimentAxisQuery"))
      self$query <- query
    },

    is_cached = function(axis) {
      stopifnot(axis %in% c("obs", "var"))
      !is.null(switch(axis,
        obs = !is.null(private$cached_obs),
        var = !is.null(private$cached_var)
      ))
    },

    preload = function(pool) {
      if (!is.null(private$cached_obs) && !is.null(private$cached_var)) {
        return(invisible(NULL))
      }

      # TODO: Use futures to parallelize preloading of obs and var joinids
      self$obs()
      self$var()
    },

    obs = function() {
      if (is.null(private$cached_obs)) {
        spdl::info("[JoinIDCache] Loading obs joinids")
        private$cached_obs <- private$load_joinids(
          df = self$query$obs_df,
          axis_query = self$query$obs_query
        )
      }
      private$cached_obs
    },

    set_obs = function(val) {
      private$cached_obs <- val
    },

    var = function() {
      if (is.null(private$cached_var)) {
        spdl::info("[JoinIDCache] Loading var joinids")
        private$cached_var <- private$load_joinids(
          df = self$query$var_df,
          axis_query = self$query$var_query
        )
      }
      private$cached_var
    },

    set_var = function(val) {
      private$cache_var <- val
    }
  ),

  private = list(
    cached_obs = NULL,
    cached_var = NULL,

    # Load joinids from the dataframe corresponding to the axis query
    # @return [`arrow::ChunkedArray`] of joinids
    load_joinids = function(df, axis_query) {
      stopifnot(
        inherits(df, "SOMADataFrame"),
        inherits(axis_query, "SOMAAxisQuery")
      )
      tbl <- df$read(
        coords = axis_query$coords,
        value_filter = axis_query$value_filter,
        column_names = "soma_joinid",
      )$concat()
      tbl$soma_joinid
    }
  )
)
