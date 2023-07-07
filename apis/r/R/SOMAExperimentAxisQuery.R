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

    #' @description Retrieve obs \link{TableReadIter}
    #' @param column_names A character vector of column names to retrieve
    obs = function(column_names = NULL) {
      obs_query <- self$obs_query
      self$obs_df$read(
        coords = recursively_make_integer64(obs_query$coords),
        value_filter = obs_query$value_filter,
        column_names = column_names
      )
    },

    #' @description Retrieve var [`arrow::Table`]
    #' @param column_names A character vector of column names to retrieve
    var = function(column_names = NULL) {
      var_query <- self$var_query
      self$var_df$read(
        coords = recursively_make_integer64(var_query$coords),
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
      obs_ft <- self$obs(obs_column_names)$concat()
      var_ft <- self$var(var_column_names)$concat()

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

    #' @description Retrieve a collection layer as a sparse matrix with named
    #' dimensions.
    #'
    #' Load any layer from the `X`, `obsm`, `varm`, `obsp`, or `varp`
    #' collections as a [sparse matrix][Matrix::sparseMatrix-class].
    #'
    #' By default the matrix dimensions are named using the `soma_joinid` values
    #' in the specified layer's dimensions (e.g., `soma_dim_0`). However,
    #' dimensions can be named using values from any `obs` or `var` column that
    #' uniquely identifies each record by specifying the `obs_index` and
    #' `var_index` arguments.
    #'
    #' For layers in \code{obsm} or \code{varm}, the column axis (the axis not
    #' indexed by \dQuote{\code{obs}} or \dQuote{\code{var}}) is set to the
    #' range of values present in \dQuote{\code{soma_dim_1}}; this ensures
    #' that gaps in this axis are preserved (eg. when a query for
    #' \dQuote{\code{obs}} that results in selecting entries that are all zero
    #' for a given PC)
    #'
    #' @param collection The [`SOMACollection`] containing the layer of
    #' interest, either: `"X"`, `"obsm"`, `"varm"`, `"obsp"`, or `"varp"`.
    #' @param layer_name Name of the layer to retrieve from the `collection`.
    #' @param obs_index,var_index Name of the column in `obs` or `var`
    #' (`var_index`) containing values that should be used as dimension labels
    #' in the resulting matrix. Whether the values are used as row or column
    #' labels depends on the selected `collection`:
    #'
    #' | Collection | `obs_index`          | `var_index`          |
    #' |-----------:|:---------------------|:---------------------|
    #' | `X`        | row names            | column names         |
    #' | `obsm`     | row names            | ignored              |
    #' | `varm`     | ignored              | row names            |
    #' | `obsp`     | row and column names | ignored              |
    #' | `varp`     | ignored              | row and column names |
    #' @return A [`Matrix::sparseMatrix-class`]
    to_sparse_matrix = function(
      collection, layer_name, obs_index = NULL, var_index = NULL
    ) {
      stopifnot(
        assert_subset(
          x = collection,
          y = c("X", "obsm", "varm", "obsp", "varp"),
          type = "collection"
        ),

        "Must specify a single layer name" = is_scalar_character(layer_name),
        assert_subset(layer_name, self$ms[[collection]]$names(), "layer"),

        "Must specify a single obs index" =
          is.null(obs_index) || is_scalar_character(obs_index),
        assert_subset(obs_index, self$obs_df$colnames(), "column"),

        "Must specify a single var index" =
          is.null(var_index) || is_scalar_character(var_index),
        assert_subset(var_index, self$var_df$colnames(), "column")
      )

      # Retrieve and validate obs/var indices
      obs_labels <- var_labels <- NULL
      if (!is.null(obs_index)) {
        if (collection %in% c("varm", "varp")) {
          spdl::warn("The obs_index is ignored for {} collections", collection)
        } else {
          obs_labels <- self$obs(column_names = obs_index)$concat()[[1]]$as_vector()
        }
        stopifnot(
          "All obs_index values must be unique" = anyDuplicated(obs_labels) == 0
        )
      }

      if (!is.null(var_index)) {
        if (collection %in% c("obsm", "obsp")) {
          spdl::warn("The var_index is ignored for {} collections", collection)
        } else {
          var_labels <- self$var(column_names = var_index)$concat()[[1]]$as_vector()
        }
        stopifnot(
          "All var_index values must be unique" = anyDuplicated(var_labels) == 0
        )
      }

      # Construct coordinates
      coords <- switch(collection,
        X = list(
          soma_dim_0 = self$obs_joinids(),
          soma_dim_1 = self$var_joinids()
        ),
        obsm = list(
          soma_dim_0 = self$obs_joinids()
        ),
        varm = list(
          soma_dim_0 = self$var_joinids()
        ),
        obsp = list(
          soma_dim_0 = self$obs_joinids(),
          soma_dim_1 = self$obs_joinids()
        ),
        varp = list(
          soma_dim_0 = self$var_joinids(),
          soma_dim_1 = self$var_joinids()
        )
      )

      # TODO: Coords must be vectors until read() supports arrow arrays
      coords <- lapply(coords, function(x) x$as_vector())

      # Retrieve coo arrow table with query result
      layer <- self$ms[[collection]]$get(layer_name)
      if (inherits(layer, "SOMADenseNDArray")) {
        tbl <- layer$read_arrow_table(coords = coords)
      } else {
        tbl <- layer$read(coords = coords)$tables()$concat()
      }


      # Reindex the coordinates
      # Constructing a matrix with the joinids produces a matrix with
      # the same shape as the original array, which is not we want. To create
      # a matrix containing only values in the query result we need to
      # reindex the coordinates.
      mat_coords <- switch(collection,
        X = list(
          i = self$indexer$by_obs(tbl$soma_dim_0),
          j = self$indexer$by_var(tbl$soma_dim_1)
        ),
        obsm = list(
          i = self$indexer$by_obs(tbl$soma_dim_0),
          j = tbl$soma_dim_1
        ),
        varm = list(
          i = self$indexer$by_var(tbl$soma_dim_0),
          j = tbl$soma_dim_1
        ),
        obsp = list(
          i = self$indexer$by_obs(tbl$soma_dim_0),
          j = self$indexer$by_obs(tbl$soma_dim_1)
        ),
        varp = list(
          i = self$indexer$by_var(tbl$soma_dim_0),
          j = self$indexer$by_var(tbl$soma_dim_1)
        )
      )

      # Construct the dimension names
      dim_names <- switch(collection,
        X = list(obs_labels, var_labels),
        obsm = {
          soma_dim_1 <- range(tbl$soma_dim_1$as_vector())
          list(obs_labels, seq(min(soma_dim_1), max(soma_dim_1)))
        },
        varm = {
          soma_dim_1 <- range(tbl$soma_dim_1$as_vector())
          list(var_labels, seq(min(soma_dim_1), max(soma_dim_1)))
        },
        obsp = list(obs_labels, obs_labels),
        varp = list(var_labels, var_labels)
      )

      # Use joinids if the dimension names are empty
      dim_names <- Map("%||%", dim_names, coords)

      Matrix::sparseMatrix(
        i = mat_coords$i$as_vector(),
        j = mat_coords$j$as_vector(),
        x = tbl$soma_data$as_vector(),
        index1 = FALSE,
        dims = vapply_int(dim_names, length),
        dimnames = dim_names,
        repr = "T"
      )
    },

    #' @description Loads the query as a \code{\link[SeuratObject]{Seurat}} object
    #'
    #' @param X_layers \Sexpr[results=rd]{tiledbsoma:::rd_outgest_xlayers()}
    #' @param obs_index \Sexpr[results=rd]{tiledbsoma:::rd_outgest_index()}
    #' @param var_index \Sexpr[results=rd]{tiledbsoma:::rd_outgest_index(axis = 'var')}
    #' @param obs_column_names \Sexpr[results=rd]{tiledbsoma:::rd_outgest_metadata_names()}
    #' @param var_column_names \Sexpr[results=rd]{tiledbsoma:::rd_outgest_metadata_names(axis = 'var')}
    #' @param obsm_layers \Sexpr[results=rd]{tiledbsoma:::rd_outgest_mlayers()}
    #' @param varm_layers \Sexpr[results=rd]{tiledbsoma:::rd_outgest_mlayers(axis = 'varm')}
    #' @param obsp_layers \Sexpr[results=rd]{tiledbsoma:::rd_outgest_players()}
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
        self$obs(obs_index)$concat()$GetColumnByName(obs_index)$as_vector()
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
          x = self$obs(obs_column_names)$concat()$to_data_frame()
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
    #' @param X_layers \Sexpr[results=rd]{tiledbsoma:::rd_outgest_xlayers()}
    #' @param obs_index \Sexpr[results=rd]{tiledbsoma:::rd_outgest_index()}
    #' @param var_index \Sexpr[results=rd]{tiledbsoma:::rd_outgest_index(axis = 'var')}
    #' @param var_column_names \Sexpr[results=rd]{tiledbsoma:::rd_outgest_metadata_names(axis = 'var')}
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
        self$var(var_index)$concat()$GetColumnByName(var_index)$as_vector()
      }
      cells <- if (is.null(obs_index)) {
        paste0('cell', self$obs_joinids()$as_vector())
      } else {
        obs_index <- match.arg(
          arg = obs_index,
          choices = self$obs_df$attrnames()
        )
        self$obs(obs_index)$concat()$GetColumnByName(obs_index)$as_vector()
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
        var <- as.data.frame(self$var(var_column_names)$concat()$to_data_frame())
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
    #' @param obs_index \Sexpr[results=rd]{tiledbsoma:::rd_outgest_index()}
    #' @param var_index \Sexpr[results=rd]{tiledbsoma:::rd_outgest_index(axis = 'var')}
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

      spdl::info("Reading obsm layer '{}' into memory", obsm_layer)
      warn_if_dense("obsm", self$ms$obsm$get(obsm_layer))
      embed_mat <- self$to_sparse_matrix(
        collection = "obsm",
        layer_name = obsm_layer,
        obs_index = obs_index
      )

      # Set matrix names
      if (is.null(obs_index)) {
        rownames(embed_mat) <- paste0("cell", rownames(embed_mat))
      }
      colnames(embed_mat) <- paste0(key, seq_len(ncol(embed_mat)))
      spdl::debug("Converting '{}' dgTMatrix to matrix", obsm_layer)
      embed_mat <- as.matrix(embed_mat)

      # Autoset loadings if needed
      if (is.null(varm_layer) || isTRUE(varm_layer)) {
        if (seurat %in% c(names(ms_load), ms_load)) {
          varm_layer <- seurat
        }
      }

      # Read in feature loadings (if present)
      load_mat <- NULL
      if (is_scalar_character(varm_layer)) {
        # Translate Seurat name to AnnData name
        if (varm_layer %in% names(ms_load)) {
          varm_layer <- ms_load[varm_layer]
        }

        spdl::info("Reading varm layer '{}' into memory", varm_layer)
        warn_if_dense("varm", self$ms$varm$get(varm_layer))
        load_mat <- self$to_sparse_matrix(
          collection = "varm",
          layer_name = varm_layer,
          var_index = var_index
        )

        # Set matrix names
        if (is.null(var_index)) {
          rownames(load_mat) <- paste0('feature', rownames(load_mat))
        }
        colnames(load_mat) <- paste0(key, seq_len(ncol(load_mat)))
        spdl::debug("Converting '{}' dgTMatrix to matrix", varm_layer)
        load_mat <- as.matrix(load_mat)

        if (!is.null(embed_mat) && ncol(load_mat) != ncol(embed_mat)) {
          stop("The loadings do not match the embeddings", call. = FALSE)
        }
      }

      # Create the DimReduc
      SeuratObject::CreateDimReducObject(
        embeddings = embed_mat,
        loadings = load_mat %||% methods::new('matrix'),
        assay = private$.measurement_name,
        global = !seurat %in% c('pca', 'ica'),
        key = key
      )
    },
    #' @description Loads the query as a Seurat \link[SeuratObject:Graph]{graph}
    #'
    #' @param obsp_layer Name of array in \code{obsp} to load as the graph
    #' @param obs_index \Sexpr[results=rd]{tiledbsoma:::rd_outgest_index()}
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

      # Retrieve the named TsparseMatrix
      mat <- self$to_sparse_matrix(
        collection = "obsp",
        layer_name = obsp_layer,
        obs_index = obs_index
      )

      # Convert to Seurat graph by way of a CsparseMatrix
      mat <- as(as(mat, "CsparseMatrix"), "Graph")

      if (is.null(obs_index)) {
        dimnames(mat) <- lapply(dimnames(mat), function(x) paste0("cell", x))

      }

      SeuratObject::DefaultAssay(mat) <- private$.measurement_name
      validObject(mat)
      return(mat)
    },
    #' @description Loads the query as a
    #' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
    #'
    #' @param X_layers \Sexpr[results=rd]{tiledbsoma:::rd_outgest_xlayers('sce')}
    #' @param obs_index \Sexpr[results=rd]{tiledbsoma:::rd_outgest_index('sce')}
    #' @param var_index \Sexpr[results=rd]{tiledbsoma:::rd_outgest_index('sce', 'var')}
    #' @param obsm_layers \Sexpr[results=rd]{tiledbsoma:::rd_outgest_mlayers('sce')}
    #' @param obs_column_names \Sexpr[results=rd]{tiledbsoma:::rd_outgest_metadata_names('sce')}
    #' @param var_column_names \Sexpr[results=rd]{tiledbsoma:::rd_outgest_metadata_names('sce', 'var')}
    #' @param obsp_layers \Sexpr[results=rd]{tiledbsoma:::rd_outgest_players('sce')}
    #' @param varp_layers \Sexpr[results=rd]{tiledbsoma:::rd_outgest_players('sce', 'varp')}
    #'
    #' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
    #'
    to_single_cell_experiment = function(
      X_layers = NULL,
      obs_index = NULL,
      var_index = NULL,
      obs_column_names = NULL,
      var_column_names = NULL,
      obsm_layers = NULL,
      # Omission of `varm_layers` parameter is purposeful as
      # SCE objects do not support `varm_layers`
      obsp_layers = NULL,
      varp_layers = NULL
    ) {
      .check_sce_installed()
      stopifnot(
        "'X_layers' must be a character vector" = is_character_or_null(X_layers),
        "'obs_index' must be a single character value" = is.null(obs_index) ||
          (is_scalar_character(obs_index) && !is.na(obs_index)),
        "'var_index' must be a single character value" = is.null(var_index) ||
          (is_scalar_character(var_index) && !is.na(var_index)),
        "'obs_column_names' must be a character vector" = is.null(obs_column_names) ||
          is.character(obs_column_names) ||
          is_scalar_logical(obs_column_names),
        "'var_column_names' must be a character vector" = is.null(var_column_names) ||
          is.character(var_column_names) ||
          is_scalar_logical(var_column_names),
        "'obsm_layers' must be a character vector" = is.null(obsm_layers) ||
          is.character(obsm_layers) ||
          is_scalar_logical(obsm_layers),
        "'obsp_layers' must be a character vector" = is.null(obsp_layers) ||
          is.character(obsp_layers) ||
          is_scalar_logical(obsp_layers),
        "'varp_layers' must be a character vector" = is.null(varp_layers) ||
          is.character(varp_layers) ||
          is_scalar_logical(varp_layers)
      )
      # Load in colData
      obs <- private$.load_df('obs', index = obs_index, column_names = obs_column_names)
      # Load in rowData
      var <- private$.load_df('var', index = var_index, column_names = var_column_names)
      # Check the layers
      X_layers <- pad_names(X_layers %||% self$ms$X$names())
      assert_subset(x = X_layers, y = self$ms$X$names(), type = 'X_layer')
      # Read in the layers
      layers <- lapply(
        X = X_layers,
        FUN = function(layer, var_ids, obs_ids) {
          mat <- Matrix::t(self$to_sparse_matrix(
            collection = 'X',
            layer_name = layer
          ))
          dimnames(mat) <- list(var_ids, obs_ids)
          return(as(object = mat, Class = 'CsparseMatrix'))
        },
        var_ids = row.names(var),
        obs_ids = row.names(obs)
      )
      names(layers) <- names(X_layers)
      # Load in reduced dimensions
      reduced_dims <- private$.load_sce_reduced_dims(
        obsm_layers = obsm_layers,
        obs_ids = row.names(obs)
      )
      # Load in the colPairs
      col_pairs <- private$.load_sce_col_pairs(
        obsp_layers = obsp_layers,
        obs_ids = row.names(obs)
      )
      # Load in the rowPairs
      row_pairs <- private$.load_sce_row_pairs(
        varp_layers = varp_layers,
        var_ids = row.names(var)
      )
      # Create the SingleCellExperiment object
      sce <- SingleCellExperiment::SingleCellExperiment(
        layers,
        reducedDims = reduced_dims,
        rowPairs = row_pairs,
        colPairs = col_pairs,
        mainExpName = private$.measurement_name
      )
      if (ncol(var)) {
        SummarizedExperiment::rowData(sce) <- as(var, 'DataFrame')
      }
      if (ncol(obs)) {
        SummarizedExperiment::colData(sce) <- as(obs, 'DataFrame')
      }
      # Validate and return
      methods::validObject(sce)
      return(sce)
    },
    #' @description Convenience method for \code{$to_single_cell_experiment()}
    #'
    #' @param X_layers \Sexpr[results=rd]{tiledbsoma:::rd_outgest_xlayers('sce')}
    #' @param obs_index \Sexpr[results=rd]{tiledbsoma:::rd_outgest_index('sce')}
    #' @param var_index \Sexpr[results=rd]{tiledbsoma:::rd_outgest_index('sce', 'var')}
    #' @param obsm_layers \Sexpr[results=rd]{tiledbsoma:::rd_outgest_mlayers('sce')}
    #' @param obs_column_names \Sexpr[results=rd]{tiledbsoma:::rd_outgest_metadata_names('sce')}
    #' @param var_column_names \Sexpr[results=rd]{tiledbsoma:::rd_outgest_metadata_names('sce', 'var')}
    #' @param obsp_layers \Sexpr[results=rd]{tiledbsoma:::rd_outgest_players('sce')}
    #' @param varp_layers \Sexpr[results=rd]{tiledbsoma:::rd_outgest_players('sce', 'varp')}
    #'
    #' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
    #'
    to_sce = function(
      X_layers = NULL,
      obs_index = NULL,
      var_index = NULL,
      obs_column_names = NULL,
      var_column_names = NULL,
      obsm_layers = NULL,
      obsp_layers = NULL,
      varp_layers = NULL
    ) {
      return(self$to_single_cell_experiment(
        X_layers = X_layers,
        obs_index = obs_index,
        var_index = var_index,
        obs_column_names = obs_column_names,
        var_column_names = var_column_names,
        obsm_layers = obsm_layers,
        obsp_layers = obsp_layers,
        varp_layers = varp_layers
      ))
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
    # Helper methods to load aspects of a measurement for the ecosystems
    # @description Load the `obs` or `var` data frames for usage with the ecosystems
    # @param df_name Name of SOMADataFrame to load; choose from `"obs"` or `"var"`
    # @param index Name of attribute in `df_name` to use as the row names for
    # the resulting data frame. By default, uses `paste0(df_name, df$soma_join_ids)`;
    # If provided, and `column_names` is `NULL`, `index` will be **excluded** from
    # the resulting data frame
    # @param column name The names of attributes in `df_name` to load; choose from:
    # - `NULL` (default) or `TRUE`: loads all attributes **except** `index`
    # - `FALSE` or `NA`: return a data frame the number of rows as present
    # in `df_name` and zero columns
    # - a character vector of names of attributes to load in
    .load_df = function(df_name = c('obs', 'var'), index = NULL, column_names = NULL) {
      stopifnot(
        is.character(df_name),
        is.null(index) || is_scalar_character(index),
        is.null(column_names) || is.character(column_names) || is_scalar_logical(column_names)
      )
      df_name <- match.arg(arg = df_name)
      switch(
        EXPR = df_name,
        obs = {
          soma_df <- self$obs_df
          soma_reader <- self$obs
          soma_joinids <- self$obs_joinids()
        },
        var = {
          soma_df <- self$var_df
          soma_reader <- self$var
          soma_joinids <- self$var_joinids()
        }
      )
      attr_names <- soma_df$attrnames()
      if (is.character(index)) {
        index <- match.arg(arg = index, choices = attr_names)
      }
      if (isTRUE(column_names)) {
        column_names <- NULL
      }
      column_names <- column_names %||% setdiff(attr_names, index)
      attrs <- if (isFALSE(column_names) || rlang::is_na(column_names)) {
        index
      } else {
        union(index, column_names)
      }
      df <- if (is.null(attrs)) {
        NULL
      } else {
        soma_reader(attrs)$concat()$to_data_frame()
      }
      ids <- if (is.null(index)) {
        paste0(df_name, soma_joinids$as_vector())
      } else {
        df[[index]]
      }
      if (is.null(df)) {
        df <- as.data.frame(matrix(nrow = length(ids), ncol = 0L))
      }
      row.names(df) <- ids
      df <- if (isFALSE(column_names) || rlang::is_na(column_names)) {
        df[, character(), drop = FALSE]
      } else {
        df[, column_names, drop = FALSE]
      }
      return(df)
    },
    .load_m_axis = function(layer, m_axis = c('obsm', 'varm'), type = "Embeddings") {
      stopifnot(
        is_scalar_character(layer),
        is.character(m_axis),
        is_scalar_character(type)
      )
      m_axis <- match.arg(arg = m_axis)
      switch(
        EXPR = m_axis,
        obsm = {
          soma_collection <- self$ms$obsm
          soma_joinids <- self$obs_joinids()
        },
        varm = {
          soma_collection <- self$ms$varm
          soma_joinids <- self$var_joinids()
        }
      )
      soma_axis <- soma_collection$get(layer)
      warn_if_dense(m_axis, soma_axis)
      mat <- self$to_sparse_matrix(collection = m_axis, layer_name = layer)
      spdl::debug("Converting '{}' dgTMatrix to matrix", layer)
      return(as.matrix(mat))
    },
    .load_p_axis = function(layer, p_axis = c('obsp', 'varp'), repr = c('C', 'T', 'R', 'D')) {
      stopifnot(
        is_scalar_character(layer),
        is.character(p_axis),
        is.character(repr)
      )
      p_axis <- match.arg(p_axis)
      repr <- match.arg(repr)
      mat <- self$to_sparse_matrix(collection = p_axis, layer_name = layer)
      return(switch(
        EXPR = repr,
        C = {
          spdl::debug("Converting '{}' TsparseMatrix to CsparseMatrix", layer)
          as(mat, 'CsparseMatrix')
        },
        R = {
          spdl::debug("Converting '{}' TsparseMatrix to RsparseMatrix", layer)
          as(mat, 'RsparseMatrix')
        },
        D = {
          spdl::debug("Converting '{}' TsparseMatrix to matrix", layer)
          as.matrix(mat)
        },
        mat
      ))
    },
    # Helper methods for loading SCE components
    .load_sce_reduced_dims = function(obsm_layers, obs_ids) {
      stopifnot(
        is.null(obsm_layers) || is.character(obsm_layers) || is_scalar_logical(obsm_layers),
        is.character(obs_ids)
      )
      if (isTRUE(obsm_layers)) {
        obsm_layers <- NULL
      }
      ms_obsm <- tryCatch(expr = self$ms$obsm$names(), error = null)
      skip_reduced_dims <- isFALSE(obsm_layers) || rlang::is_na(obsm_layers)
      if (is.null(ms_obsm)) {
        if (!skip_reduced_dims && !is.null(obsm_layers)) {
          warning(
            "Reduced dimensions requested but none were found, skipping",
            call. = FALSE,
            immediate. = TRUE
          )
        }
        skip_reduced_dims <- TRUE
      }
      if (skip_reduced_dims) {
        return(list())
      }
      obsm_layers <- obsm_layers %||% stats::setNames(
        object = ms_obsm,
        nm = .anndata_to_sce_reduced_dim(ms_obsm)
      )
      obsm_layers <- pad_names(obsm_layers)
      assert_subset(x = obsm_layers, y = ms_obsm, type = 'cell embedding')
      reduced_dims <- lapply(
        X = seq_along(obsm_layers),
        FUN = function(i) {
          layer <- obsm_layers[i]
          rd <- names(obsm_layers)[i]
          mat <- private$.load_m_axis(layer = layer, type = "Reduced dimensions")
          dimnames(mat) <- list(
            obs_ids,
            paste0(
              switch(EXPR = rd, PCA = 'PC', ICA = 'IC', TSNE = 'tSNE', rd),
              seq_len(ncol(mat))
            )
          )
          return(mat)
        }
      )
      return(stats::setNames(reduced_dims, nm = names(obsm_layers)))
    },
    .load_sce_col_pairs = function(obsp_layers, obs_ids) {
      stopifnot(
        is.null(obsp_layers) || is.character(obsp_layers) || is_scalar_logical(obsp_layers),
        is.character(obs_ids)
      )
      if (isTRUE(obsp_layers)) {
        obsp_layers <- NULL
      }
      ms_obsp <- tryCatch(expr = self$ms$obsp$names(), error = null)
      skip_col_pairs <- isFALSE(obsp_layers) || rlang::is_na(obsp_layers)
      if (is.null(ms_obsp)) {
        if (!skip_col_pairs && !is.null(obsp_layers)) {
          warning(
            "colPairs requested but none were found, skipping",
            call. = FALSE,
            immediate. = TRUE
          )
        }
        skip_col_pairs <- TRUE
      }
      if (skip_col_pairs) {
        return(list())
      }
      obsp_layers <- obsp_layers %||% ms_obsp
      obsp_layers <- pad_names(obsp_layers)
      assert_subset(x = obsp_layers, y = ms_obsp, type = 'nearest neighbor graph')
      col_pairs <- lapply(
        X = obsp_layers,
        FUN = function(layer) {
          mat <- private$.load_p_axis(layer)
          dimnames(mat) <- list(obs_ids, obs_ids)
          return(.mat_to_hits(mat))
        }
      )
      return(stats::setNames(col_pairs, names(obsp_layers)))
    },
    .load_sce_row_pairs = function(varp_layers, var_ids) {
      stopifnot(
        is.null(varp_layers) || is.character(varp_layers) || is_scalar_logical(varp_layers),
        is.character(var_ids)
      )
      if (isTRUE(varp_layers)) {
        varp_layers <- NULL
      }
      ms_varp <- tryCatch(expr = self$ms$varp$names(), error = null)
      skip_row_pairs <- isFALSE(varp_layers) || rlang::is_na(varp_layers)
      if (is.null(ms_varp)) {
        if (!skip_row_pairs && !is.null(varp_layers)) {
          warning(
            "rowPairs requested but none were found, skipping",
            call. = FALSE,
            immediate. = TRUE
          )
        }
        skip_row_pairs <- TRUE
      }
      if (skip_row_pairs) {
        return(list())
      }
      varp_layers <- varp_layers %||% ms_varp
      varp_layers <- pad_names(varp_layers)
      assert_subset(x = varp_layers, y = ms_varp, type = 'feature network')
      row_pairs <- lapply(
        X = varp_layers,
        FUN = function(layer) {
          mat <- private$.load_p_axis(layer, p_axis = 'varp')
          dimnames(mat) <- list(var_ids, var_ids)
          return(.mat_to_hits(mat))
        }
      )
      return(stats::setNames(row_pairs, names(varp_layers)))
    },
    # Helper methods for loading Seurat assays
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
        # TODO: potentially replace with public$to_sparse_matrix()
        dmat <- private$.as_matrix(
          table = self$X(data)$tables()$concat(),
          repr = 'C',
          transpose = TRUE
        )
        dimnames(dmat) <- dnames
        obj <- SeuratObject::CreateAssayObject(data = dmat)
      }
      # Add the `counts` slot
      if (is_scalar_character(counts)) {
        cmat <- private$.as_matrix(
          table = self$X(counts)$tables()$concat(),
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
          table = self$X(scale_data)$tables()$concat(),
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
