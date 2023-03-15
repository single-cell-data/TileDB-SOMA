#' @importFrom rlang is_na
#' @importFrom methods new
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
        coords = obs_query$coords,
        value_filter = obs_query$value_filter,
        column_names = column_names
      )
    },

    #' @description Retrieve var [`arrow::Table`]
    #' @param column_names A character vector of column names to retrieve
    var = function(column_names = NULL) {
      var_query <- self$var_query
      self$var_df$read(
        coords = var_query$coords,
        value_filter = var_query$value_filter,
        column_names = column_names
      )
    },

    #' @description Retrieve `soma_joinids` as an [`arrow::Array`] for `obs`.
    obs_joinids = function() {
      arrow::concat_arrays(private$.joinids$obs())
    },

    #' @description Retrieve `soma_joinids` as an [`arrow::Array`] for `var`.
    var_joinids = function() {
      arrow::concat_arrays(private$.joinids$var())
    },

    #' @description Retrieves an `X` layer as an [`arrow::Table`].
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

      # TODO: Stop converting to vectors when SOMAReader supports arrow arrays
      x_layer$read_arrow_table(coords = list(
        self$obs_joinids()$as_vector(),
        self$var_joinids()$as_vector()
      ))
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
          x_array$read_arrow_table(coords = list(
            self$obs_joinids()$as_vector(),
            self$var_joinids()$as_vector()
          ))
        }
      )

      SOMAAxisQueryResult$new(
        obs = obs_ft, var = var_ft, X_layers = x_matrices
      )
    },
    #' @description ...
    #'
    #' @param X_layers ...
    #'
    #' @return ...
    #'
    to_seurat = function(X_layers = c(counts = 'counts', data = 'logcounts')) {
      .check_seurat_installed()
      .NotYetImplemented()
    },
    #' @description Loads the query as a Seurat \code{\link[SeuratObject]{Assay}}
    #'
    #' @param X_layers A named character of X layers to add to the Seurat assay;
    #' names should be one of:
    #' \itemize{
    #'  \item \dQuote{\code{counts}} to add the layer as \code{counts}
    #'  \item \dQuote{\code{data}} to add the layer as \code{data}
    #'  \item \dQuote{\code{scale.data}} to add the layer as \code{scale.data}
    #' }
    #' At least one of \dQuote{\code{counts}} or \dQuote{\code{data}} is required
    #' @param cells_index Name of column in \code{obs} to add as cell names
    #' @param features_index Name of column in \code{var} to add as feature names
    #' @param var_column_names Names of columns in \code{var} to add as
    #' feature-level meta data
    #'
    #' @return An \code{\link[SeuratObject]{Assay}} object
    #'
    to_seurat_assay = function(
      X_layers = c(counts = 'counts', data = 'logcounts'),
      cells_index = NULL,
      features_index = NULL,
      var_column_names = NULL
    ) {
      version <- 'v3'
      .check_seurat_installed()
      stopifnot(
        "'X_layers' must be a named character vector" = is.character(X_layers) &&
          is_named(X_layers, allow_empty = FALSE),
        "'version' must be a single character value" = is_scalar_character(version),
        "'cells_index' must be a single character value" = is.null(cells_index) ||
          (is_scalar_character(cells_index) && !is.na(cells_index)),
        "'features_index' must be a single character value" = is.null(features_index) ||
          (is_scalar_character(features_index) && !is.na(features_index)),
        "'var_column_names' must be a character vector" = is.null(var_column_names) ||
          is.character(var_column_names) ||
          (is.logical(var_column_names) && length(var_column_names) == 1L)
      )
      match.arg(version, choices = 'v3')
      features <- if (is.null(features_index)) {
        paste0('feature', self$var_joinids())
      } else {
        features_index <- match.arg(
          arg = features_index,
          choices = self$var_df$attrnames()
        )
        self$var(features_index)$GetColumnByName(features_index)$as_vector()
      }
      cells <- if (is.null(cells_index)) {
        paste0('cell', self$obs_joinids())
      } else {
        cells_index <- match.arg(
          arg = cells_index,
          choices = self$obs_df$attrnames()
        )
        self$obs(cells_index)$GetColumnByName(cells_index)$as_vector()
      }
      # Check the layers
      layers <- X_layers[X_layers %in% self$ms$X$names()]
      if (!length(layers)) {
        stop("None of the requested X_layers can be found", call. = FALSE)
      } else if (length(layers) != length(X_layers)) {
        warning(
          paste(
            strwrap(paste(
              "The following layers cannot be found in this ExperimentQuery:",
              paste(sQuote(setdiff(x = X_layers, y = layers)), collapse = ', ')
            )),
            collapse = '\n'
          ),
          call. = FALSE,
          immediate. = TRUE
        )
      }
      # Read in the assay
      obj <- switch(
        EXPR = version,
        v3 = {
          if (!all(names(x = X_layers) %in% c('counts', 'data', 'scale.data'))) {
            stop(
              "The names of 'X_layers' must one or more of 'counts', 'data', and 'scale.data'",
              call. = FALSE
            )
          }
          private$.to_seurat_assay_v3(
            counts = tryCatch(expr = layers[['counts']], error = \(...) NULL),
            data = tryCatch(expr = layers[['data']], error = \(...) NULL),
            scale_data = tryCatch(expr = layers[['scale.data']], error = \(...) NULL),
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
        y = features_index
      )
      if (!(isFALSE(var_column_names) || rlang::is_na(var_column_names))) {
        var <- as.data.frame(self$var(var_column_names)$to_data_frame())
        row.names(var) <- features
        obj[[names(var)]] <- var
      }
      return(obj)
    },
    #' @description Loads the query as a Seurat
    #' \link[SeuratObject:DimReduc]{dimensional reduction}
    #'
    #' @param embeddings Name of array in \code{obsm} to load as the
    #' cell embeddings
    #' @param loadings Name of the array in \code{varm} to load as the  feature
    #' loadings; will try to determine \code{loadings} from \code{embeddings}
    #' @param cells_index Name of column in \code{obs} to add as cell names
    #' @param features_index Name of column in \code{var} to add as feature names
    #'
    #' @return A \code{\link[SeuratObject]{DimReduc}} object
    #'
    to_seurat_reduction = function(
      embeddings,
      loadings = NULL,
      cells_index = NULL,
      features_index = NULL
    ) {
      .check_seurat_installed()
      stopifnot(
        "'embeddings' must be a single character value" = is_scalar_character(embeddings),
        "'loadings' must be a single character value" = is.null(loadings) ||
          is_scalar_character(loadings) ||
          is_scalar_logical(loadings),
        "one of 'embeddings' or 'loadings' must be provided" =
          (is_scalar_character(embeddings) || is_scalar_logical(embeddings)) ||
          (is_scalar_character(loadings) || is_scalar_logical(loadings)),
        "'cells_index' must be a single character value" = is.null(cells_index) ||
          (is_scalar_character(cells_index) && !is.na(cells_index)),
        "'features_index' must be a single character value" = is.null(features_index) ||
          (is_scalar_character(features_index) && !is.na(features_index))
      )
      # Check embeddings/loadings
      ms_embed <- tryCatch(
        expr = self$ms$obsm$names(),
        error = \(...) NULL
      )
      ms_load <- tryCatch(
        expr = self$ms$varm$names(),
        error = \(...) NULL
      )
      if (is.null(ms_embed) && is.null(ms_load)) {
        warning("No reductions present", call. = FALSE)
        return(NULL)
      }
      if (is.null(ms_embed)) {
        stop("No embeddings present", call. = FALSE)
      }
      names(ms_embed) <- .anndata_to_seurat_reduc(ms_embed)
      if (is.null(ms_load) && !is.null(loadings)) {
        msg <- "No loadings present"
        if (is.null(embeddings)) {
          stop(msg, call. = FALSE)
        }
        warning(msg, call. = FALSE, immediate. = TRUE)
        loadings <- NULL
      } else {
        names(ms_load) <- .anndata_to_seurat_reduc(ms_load, 'loadings')
      }
      # Check provided names
      if (!embeddings %in% c(ms_embed, names(ms_embed))) {
        stop("Cannot find embeddings ", sQuote(embeddings), call. = FALSE)
      }
      if (is_scalar_character(loadings) && !loadings %in% c(ms_load, names(ms_load))) {
        stop("Cannot find loadings ", sQuote(loadings), call. = FALSE)
      }
      # Find Seurat name
      seurat <- c(
        .anndata_to_seurat_reduc(embeddings),
        tryCatch(
          expr = .anndata_to_seurat_reduc(loadings, 'loadings'),
          error = \(...) NULL
        )
      )
      if (length(seurat) == 2L && !identical(x = seurat[1L], y = seurat[2L])) {
        msg <- paste0(
          "The embeddings requested (",
          sQuote(embeddings),
          ") do not match the loadings requested (",
          sQuote(loadings),
          "); using the embeddings to create a Seurat name (",
          sQuote(seurat[1L]),
          ")"
        )
        warning(
          paste(strwrap(msg), collapse = '\n'),
          call. = FALSE,
          immediate. = TRUE
        )
      }
      seurat <- seurat[1L]
      # Create a Seurat key
      key <- SeuratObject::Key(
        object = switch(EXPR = seurat, pca = 'PC', tsne = 'tSNE', toupper(seurat)),
        quiet = TRUE
      )
      # Read in cell embeddings
      # Translate Seurat name to AnnData name
      if (embeddings %in% names(ms_embed)) {
        embeddings <- ms_embed[embeddings]
      }
      # Get cell names
      cells <- if (is.null(cells_index)) {
        paste0('cell', self$obs_joinids())
      } else {
        cells_index <- match.arg(
          arg = cells_index,
          choices = self$obs_df$attrnames()
        )
        self$obs(cells_index)$GetColumnByName(cells_index)$as_vector()
      }
      embed <- self$ms$obsm$get(embeddings)
      embed_mat <- embed$read_dense_matrix(list(
        self$obs_joinids()$as_vector(),
        seq_len(as.integer(embed$shape()[2L])) - 1L
      ))
      # Set matrix names
      rownames(embed_mat) <- cells
      colnames(embed_mat) <- paste0(key, seq_len(ncol(embed_mat)))
      # Autoset loadings if needed
      if (is.null(loadings) || isTRUE(loadings)) {
        if (seurat %in% c(names(ms_load), ms_load)) {
          loadings <- seurat
        }
      }
      # Read in feature loadings
      if (is_scalar_character(loadings)) {
        # Translate Seurat name to AnnData name
        if (loadings %in% names(ms_load)) {
          loadings <- ms_load[loadings]
        }
        # Get feature names
        features <- if (is.null(features_index)) {
          paste0('feature', self$var_joinids())
        } else {
          features_index <- match.arg(
            arg = features_index,
            choices = self$var_df$attrnames()
          )
          self$var(features_index)$GetColumnByName(features_index)$as_vector()
        }
        loads <- self$ms$varm$get(loadings)
        load_mat <- loads$read_dense_matrix(list(
          self$var_joinids()$as_vector(),
          seq_len(as.integer(loads$shape()[2L])) - 1L
        ))
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
        global = seurat != 'pca',
        key = key
      ))
    },
    #' @description ...
    #'
    #' @return ...
    #'
    to_seurat_graph = function() {
      .NotYetImplemented()
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
      as_matrix <- function(lyr, repr = 'C') {
        repr <- match.arg(arg = repr, choices = c('C', 'R', 'T', 'D'))
        obs <- self$X(lyr)$GetColumnByName('soma_dim_0')$as_vector()
        var <- self$X(lyr)$GetColumnByName('soma_dim_1')$as_vector()
        mat <- Matrix::sparseMatrix(
          i = self$indexer$by_obs(obs)$as_vector() + 1L,
          j = self$indexer$by_var(var)$as_vector() + 1L,
          x = self$X(lyr)$GetColumnByName('soma_data')$as_vector(),
          repr = switch(EXPR = repr, D = 'T', repr)
        )
        mat <- Matrix::t(mat)
        if (repr == 'D') {
          mat <- as.matrix(mat)
        }
        return(mat)
      }
      if (!length(x = cells) == self$n_obs) {
        stop("'cells' must have a length of ", self$n_obs, call. = FALSE)
      }
      if (!length(x = features) == self$n_vars) {
        stop("'features' must have a length of ", self$n_vars, call. = FALSE)
      }
      dnames <- list(features, cells)
      # Read in `data` slot
      if (is_scalar_character(data)) {
        dmat <- as_matrix(lyr = data, repr = 'C')
        dimnames(dmat) <- dnames
        obj <- SeuratObject::CreateAssayObject(data = dmat)
      }
      # Add the `counts` slot
      if (is_scalar_character(counts)) {
        cmat <- as_matrix(lyr = counts, repr = 'C')
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
        smat <- as_matrix(lyr = scale_data, repr = 'D')
        dimnames(smat) <- dnames
        obj <- SeuratObject::SetAssayData(
          object = obj,
          slot = 'scale.data',
          new.data = smat
        )
      }
      # Return the assay
      methods::validObject(obj)
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
      )
      tbl$soma_joinid
    }
  )
)
