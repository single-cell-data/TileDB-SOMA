#' @description Single-axis dataframe query with coordinates and a value filter.
#' [lifecycle: experimental]
#' Per dimension, the AxisQuery can have value of:
#'
#' * None - all data
#' * Coordinates - a set of coordinates on the axis dataframe index, expressed
#'   in any type or format supported by ``DataFrame.read()``.
#' * A SOMA ``value_filter`` across columns in the axis dataframe, expressed as
#'   string
#' * Or, a combination of coordinates and value filter.

AxisQuery <- R6Class(
  classname = "AxisQuery",
  public = list(
    value_filter = NULL,
    coords = NULL,

    initialize = function(value_filter = NULL, coords = NULL) {
      self$coords <- validate_read_coords(coords, dimnames = NULL)
      self$value_filter <- validate_read_value_filter(value_filter)
    }
  )
)

#' Validate read coordinates
#' Ensures that coords is a named list
#' @param coords A list of coordinates
#' @param dimnames vector of array dimension names
#' @noRd
validate_read_coords <- function(coords, dimnames = NULL) {
  stopifnot(
    "'coords' must be a named list" =
        is.null(coords) || is_named_list(coords)
  )

  if (!is.null(dimnames)) {
    stopifnot(
      "'dimnames' must be a character vector" = is.character(dimnames),
      "names of 'coords' must correspond to dimension names" =
        all(names(coords) %in% dimnames)
    )
  }
  coords
}

#' Validate read/query value filter
#' @noRd
validate_read_value_filter <- function(value_filter) {
  stopifnot(
    "'value_filter' must be a scalar character" =
      is.null(value_filter) || is_scalar_character(value_filter)
    )
  value_filter
}


#' @importFrom arrow concat_arrays
ExperimentAxisQuery <- R6::R6Class(
  classname = "ExperimentAxisQuery",
  public = list(

    experiment = NULL,
    measurement_name = NULL,
    matrix_axis_query = NULL,
    joinids = NULL,
    indexer = NULL,
    threadpool = NULL,

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

      self$experiment <- experiment
      self$measurement_name <- measurement_name
      self$matrix_axis_query <- list(
        obs = obs_query %||% AxisQuery$new(),
        var = var_query %||% AxisQuery$new()
      )
      self$joinids <- JoinIDCache$new(self)
      self$indexer <- AxisIndexer$new(self)
      self$threadpool <- NULL
    },

    #' @description Retrieve obs [`arrow::Table`]
    #' @param column_names A character vector of column names to retrieve
    obs = function(column_names = NULL) {
      obs_query <- self$matrix_axis_query$obs
      self$obs_df$read(
        coords = obs_query$coords,
        value_filter = obs_query$value_filter,
        column_names = column_names
      )
    },

    #' @description Retrieve var [`arrow::Table`]
    #' @param column_names A character vector of column names to retrieve
    var = function(column_names = NULL) {
      var_query <- self$matrix_axis_query$var
      self$var_df$read(
        coords = var_query$coords,
        value_filter = var_query$value_filter,
        column_names = column_names
      )
    },

    #' @description Retrieve obs `soma_joinids`` as an [`arrow::Array`]
    obs_joinids = function() {
      arrow::concat_arrays(self$joinids$obs())
    },

    #' @description Retrieve obs `soma_joinids`` as an [`arrow::Array`]
    var_joinids = function() {
      arrow::concat_arrays(self$joinids$var())
    },

    #' @description Retrieves an `X` layer as an [`arrow::Table`].
    #' @param layer_name The name of the layer to retrieve.
    X = function(layer_name) {
      stopifnot(
        "Must specify a layer name" = !missing(layer_name),
        assert_subset(layer_name, self$ms$names(), "layer")
      )

      self$ms

      x_layer <- self$ms$get(layer_name)
      stopifnot("X Layers must be SparseNDArrays" = inherits(x_layer, "data.SparseNDArray"))
      browser()
      # x_layer$
      return(x_layer)
    }
  ),

  active = list(
    #' @field The number of `obs` axis query results.
    n_obs = function() {
      length(self$obs_joinids())
    },

    #' @field The number of `var` axis query results.
    n_vars = function() {
      length(self$var_joinids())
    },

    #' @field The `obs` [`SOMADataFrame`] object.
    obs_df = function(value) {
      if (!missing(value)) read_only_error("obs_df")
      self$experiment$obs
    },

    #' @field The `var` [`SOMADataFrame`] object for the specified
    #' `measurement_name`.
    var_df = function(value) {
      if (!missing(value)) read_only_error("var_df")
      self$ms$var
    },

    #' @field The [`SOMAMeasurement`] object for the specified
    #' `measurement_name`.
    ms = function(value) {
      if (!missing(value)) read_only_error("ms")
      self$experiment$ms$get(self$measurement_name)
    }
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
      obs_ft <- pool$submit(function() self$obs)
      var_ft <- pool$submit(function() self$var)
      obs_ft$result()
      var_ft$result()
    },

    obs = function() {
      if (is.null(private$cached_obs)) {
        private$cached_obs <- private$load_joinids(
          df = self$query$obs_df,
          axis_query = self$query$matrix_axis_query$obs
        )
      }
      private$cached_obs
    },

    set_obs = function(val) {
      private$cached_obs <- val
    },

    var = function() {
      if (is.null(private$cached_var)) {
        private$cached_var <- private$load_joinids(
          df = self$query$var_df,
          axis_query = self$query$matrix_axis_query$var
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



#' @description Given a query, provide index-building services for obs/var
#' axis.
AxisIndexer <- R6::R6Class("AxisIndexer",
  public = list(
    #' @field query The [`ExperimentAxisQuery`] object to build indices for.
    query = NULL,

    initialize = function(query) {
      stopifnot(inherits(query, "ExperimentAxisQuery"))
      self$query <- query
    },

    obs_index = function() {
      if (is.null(private$cached_obs)) {
        private$cached_obs <- self$query$obs_joinids()
      }
      private$cached_obs
    },

    var_index = function() {
      if (is.null(private$cached_var)) {
        private$cached_var <- private$query$var_joinids()
      }
      return(private$cached_var)
    },

    by_obs = function(coords) {
      return(private$obs_index()[coords])
    },

    by_var = function(coords) {
      return(private$var_index()[coords])
    }
  ),

  private = list(
    cached_obs = NULL,
    cached_var = NULL
  )
)
