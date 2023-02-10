#' SOMAExperiment Axis Query
#' @description  Axis-based query against a [`SOMAExperiment`].
#'
#' `SOMAExperimentAxisQuery` allows easy selection and extraction of data from a
#' single [`SOMAMeasurement`] in a [`SOMAExperiment`], by `obs`/`var` (axis)
#' coordinates and/or value filter. The primary use for this class is slicing
#' [`SOMAExperiment`] `X` layers by `obs` or `var` value and/or coordinates.
#' [lifecycle: experimental]
#'
#' ## X Layer Support
#'
#' Slicing on [`SOMASparseNDArray`] `X` matrices is supported;
#' [`SOMADenseNDArray`] is not supported at this time.
#'
#' ## Result Size
#' `SOMAExperimentAxisQuery` query class assumes it can store the full result of
#' both axis dataframe queries in memory, and only provides incremental access
#' to the underlying X NDArray. Accessors such as `n_obs` and `n_vars` codify
#' this in the class.
#'
#' @importFrom arrow concat_arrays
#' @export
ExperimentAxisQuery <- R6::R6Class(
  classname = "ExperimentAxisQuery",

  public = list(
    #' @description Create a new `SOMAExperimentAxisQuery` object.
    #' @param experiment A [`SOMAExperiment`] object.
    #' @param measurement_name The name of the measurement to query.
    #' @param obs_query,var_query An [`AxisQuery`] object for the obs/var axis.
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
        is.null(obs_query) || inherits(obs_query, "AxisQuery"),
        is.null(var_query) || inherits(var_query, "AxisQuery")
      )

      private$.experiment <- experiment
      private$.measurement_name <- measurement_name
      private$.obs_query <- obs_query %||% AxisQuery$new()
      private$.var_query <- var_query %||% AxisQuery$new()
      private$.joinids <- JoinIDCache$new(self)
      private$.indexer <- AxisIndexer$new(self)
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
        "Must specify a layer name" = !missing(layer_name),
        "Must specify a single layer name" = is_scalar_character(layer_name),
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
    #' @param X_layers The name of the `X` layer(s) to read and return.
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

      AxisQueryResult$new(
        obs = obs_ft, var = var_ft, X_layers = x_matrices
      )
    }
  ),

  active = list(

    #' @field experiment The parent [`SOMAExperiment`] object.
    experiment = function(value) {
      if (!missing(value)) read_only_error("experiment")
      private$.experiment
    },

    #' @field indexer The [`SOMAIndexer`] object.
    indexer = function(value) {
      if (!missing(value)) read_only_error("indexer")
      private$.indexer
    },

    #' @field obs_query The `obs` [`AxisQuery`] object.
    obs_query = function(value) {
      if (!missing(value)) read_only_error("obs_query")
      private$.obs_query
    },

    #' @field var_query The `var` [`AxisQuery`] object.
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
    .indexer = NULL
  )
)

JoinIDCache <- R6::R6Class(
  classname = "JoinIDCache",
  public = list(
    #' @field query The [`ExperimentAxisQuery`] object to build indices for.
    query = NULL,

    initialize = function(query) {
      stopifnot(inherits(query, "ExperimentAxisQuery"))
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
        inherits(axis_query, "AxisQuery")
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

#' @description Return types for [`ExperimentAxisQuery`] `read()` method.
AxisQueryResult <- R6Class(
  classname = "AxisQueryResult",
  public = list(

    #' @description Create a new `AxisQueryResult` object.
    #' @param obs,var [`arrow::Table`] containing `obs` or `var` query slice.
    #' @param X_layers named list of [`arrow::Table`]s, one for each `X` layer.
    initialize = function(obs, var, X_layers) {
      stopifnot(
        is_arrow_table(obs),
        is_arrow_table(var),
        is_named_list(X_layers),
        all(vapply_lgl(X_layers, is_arrow_table))
      )

      private$.obs <- obs
      private$.var <- var
      private$.X_layers <- X_layers
    }
  ),

  active = list(
    #' @field obs [`arrow::Table`] containing `obs` query slice.
    obs = function(value) {
      if (!missing(value)) read_only_error("obs")
      private$.obs
    },

    #' @field var [`arrow::Table`] containing `var` query slice.
    #' `measurement_name`.
    var = function(value) {
      if (!missing(value)) read_only_error("var")
      private$.var
    },

    #' @field X_layers named list of [`arrow::Table`]s for each `X` layer.
    X_layers = function(value) {
      if (!missing(value)) read_only_error("ms")
      private$.X_layers
    }
  ),

  private = list(
    .obs = NULL,
    .var = NULL,
    .X_layers = NULL
  )
)



#' @description Index obs/var soma_joinids for a given query result.
#'
#' Retrieve the index of the given `obs` or `var` coordinates in the query
#' result. Coordinates outside of the query result will return
#' [`arrow::null()`].
AxisIndexer <- R6::R6Class("AxisIndexer",
  public = list(

    #' @description Create a new `AxisIndexer` object.
    #' @param query The [`ExperimentAxisQuery`] object to build indices for.
    initialize = function(query) {
      stopifnot(inherits(query, "ExperimentAxisQuery"))
      private$.query <- query
    },

    #' @description Get the index of the given `obs` coordinates.
    #' @param coords vector or [`arrow:Array`] of numeric coordinates.
    by_obs = function(coords) {
      arrow::match_arrow(
        x = private$.validate_coords(coords),
        table = private$.obs_index()
      )
    },

    #' @description Get the index of the given `var` coordinates.
    #' @param coords vector or [`arrow:Array`] of numeric coordinates.
    by_var = function(coords) {
      arrow::match_arrow(
        x = private$.validate_coords(coords),
        table = private$.var_index()
      )
    }
  ),

  private = list(
    .cached_obs = NULL,
    .cached_var = NULL,
    .query = NULL,

    # Retrieve index for the obs axis
    .obs_index = function() {
      if (is.null(private$.cached_obs)) {
        private$.cached_obs <- private$.query$obs_joinids()
      }
      private$.cached_obs
    },

    # Retrieve index for the var axis
    .var_index = function() {
      if (is.null(private$.cached_var)) {
        private$.cached_var <- private$.query$var_joinids()
      }
      private$.cached_var
    },

    .validate_coords = function(coords) {
      stopifnot(
        "'coords' must be a numeric vector or arrow Array" =
          is.numeric(coords) || is_arrow_array(coords)
      )
      coords
    }
  )
)
