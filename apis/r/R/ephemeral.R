#' Ephemeral Collection Base
#'
#' @description Virtual base class for ephemeral collections; ephemeral
#' collections are equivalent to
#' \link[tiledbsoma:SOMACollection]{SOMA collections} but are stored in memory
#' instead of on disk.
#'
#' @keywords internal
#'
#' @export
#'
#' @seealso Derived classes: \code{\link{EphemeralCollection}},
#' \code{\link{EphemeralMeasurement}},
#' \code{\link{EphemeralExperiment}}
#'
EphemeralCollectionBase <- R6::R6Class(
  classname = "EphemeralCollectionBase",
  inherit = SOMACollectionBase,
  public = list(
    # Override SOMAObject and SOMACollectionBase methods
    #' @description Create an ephemeral collection.
    #'
    #' @template param-dots-ignored
    #'
    initialize = function(...) {
      # Check if any arguments were passed
      # If so, warn about unused arguments for ephemeral objects
      # Python equivalent:
      # def f(*args, **kwargs):
      #   if args or kwargs:
      #     warnings.warn("argumnets passed but unused")
      #   pass
      #
      if (...length()) {
        tryCatch(
          expr = private$.ephemeral_error("custom", "and cannot be customized"),
          error = function(e) {
            warning(conditionMessage(e), call. = FALSE, immediate. = TRUE)
          }
        )
      }
      private$.data <- list()
    },

    #' @description Create a new, empty ephemeral collection.
    #'
    #' @return Returns a new ephemeral collection of class \code{class(self)}.
    #'
    create = function() {
      gen <- getAnywhere(self$class())[["objs"]][[1L]]
      if (!R6::is.R6Class(gen)) {
        stop(
          "Cannot find the class generator for ",
          sQuote(self$class()),
          call. = FALSE
        )
      }
      return(gen$new())
    },

    #' @description Dummy method for ephemeral objects for compatibility with
    #' SOMA collections.
    #'
    #' @param mode Ignored for ephemeral objects.
    #'
    #' @return Throws an error as this method is not supported by ephemeral
    #' objects.
    #'
    open = \(mode) private$.ephemeral_error("opened"),

    #' @description Dummy method for ephemeral objects for compatibility with
    #' SOMA collections.
    #'
    #' @return Invisibly returns \code{NULL}.
    #'
    close = function() {
      tryCatch(
        expr = private$.ephemeral_error("custom", "and cannot be closed"),
        error = function(e) {
          warning(conditionMessage(e), call. = FALSE, immediate. = TRUE)
        }
      )
      return(invisible(NULL))
    },

    #' @description Dummy method for ephemeral objects for compatibility with
    #' SOMA collections.
    #'
    #' @return Returns \code{FALSE} as ephemeral collections do not
    #' exist on-disk.
    #'
    exists = \() FALSE,

    #' @description Special method for printing object representation to
    #' console.
    #'
    #' @return Prints details about the ephemeral collection and invisibly
    #' returns itself.
    #'
    print = function() {
      super$print()
      if (self$length()) {
        f <- vapply(
          self$members,
          FUN = \(x) {
            if (inherits(x, "SOMACollectionBase")) "GROUP" else "ARRAY"
          },
          FUN.VALUE = character(1L),
          USE.NAMES = FALSE
        )
        members <- split(names(self$members), f = f)
        if (!is.null(members$ARRAY)) {
          cat("  arrays: ", string_collapse(sort(members$ARRAY)), "\n")
        }
        if (!is.null(members$GROUP)) {
          cat("  groups:", string_collapse(sort(members$GROUP)), "\n")
        }
      }
      return(invisible(self))
    },

    #' @description Dummy method for ephemeral objects for compatibility with
    #' SOMA collections.
    #'
    #' @param param Ignored for ephemeral objects.
    #'
    #' @return Returns \code{NULL} as ephemeral collections do not have an
    #' on-disk configuration.
    #'
    get_tiledb_config = function(param = NULL) {
      if (!is.null(param)) {
        tryCatch(
          expr = private$.ephemeral_error(
            "custom",
            "and have no TileDB configuration"
          ),
          error = function(e) {
            warning(conditionMessage(e), call. = FALSE, immediate. = TRUE)
          }
        )
      }
      return(NULL)
    },

    #' @description Retrieve the number of items in the collection.
    #'
    #' @return The length of the collection.
    #'
    length = function() {
      length(private$.data)
    },

    #' @description Retrieve the names of members. (lifecycle: maturing).
    #'
    #' @return A \code{character} vector of member names.
    #'
    names = function() {
      names(private$.data) %||% character(length = 0L)
    },

    #' @description Add object to an ephemeral collection.
    #'
    #' @param object A SOMA object (eg. \code{\link{SOMACollection}}) to add
    #' to the collection.
    #' @param name A name to add \code{object} as.
    #' @param relative Ignored for ephemeral objects.
    #'
    #' @return \[chainable] Invisibly returns \code{self} with \code{object}
    #' added as \code{name}.
    #'
    set = function(object, name = NULL, relative = NULL) {
      stopifnot(
        "Only SOMA objects may be added" = inherits(object, "SOMAObject"),
        "'name' must be a single, non-empty string" = is.null(name) ||
          (is_scalar_character(name) && nzchar(name)),
        is.null(relative) || is_scalar_logical(relative)
      )
      if (!is.null(relative)) {
        tryCatch(
          expr = private$.ephemeral_error(
            "custom",
            "so relative has no effect"
          ),
          error = function(e) {
            warning(conditionMessage(e), call. = FALSE, immediate. = TRUE)
          }
        )
      }
      name <- name %||% object$uri
      private$.data[[name]] <- object
      return(invisible(self))
    },

    #' @description Get objects from an ephemeral collection.
    #'
    #' @param name Name of object in the collection to get.
    #'
    #' @return The object named \code{name}.
    #'
    get = function(name) {
      stopifnot(is_scalar_character(name))
      name <- match.arg(arg = name, choices = self$names())
      return(private$.data[[name]])
    },

    #' @description Remove objects from an ephemeral collection.
    #'
    #' @param name Name of object to remove from the collection.
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with the object at
    #' \code{name} removed.
    #'
    remove = function(name) {
      stopifnot(is_scalar_character(name))
      name <- match.arg(arg = name, choices = self$names())
      private$.data[[name]] <- NULL
      return(invisible(self))
    },

    #' @description Dummy method for ephemeral objects for compatibility with
    #' SOMA collections.
    #'
    #' @param metadata Ignored for ephemeral objects.
    #'
    #' @return Throws an error as this method is not supported by ephemeral
    #' objects.
    #'
    set_metadata = \(metadata) private$.ephemeral_error("edited"),

    #' @description Dummy method for ephemeral objects for compatibility with
    #' SOMA collections.
    #'
    #' @param key Ignored for ephemeral objects.
    #'
    #' @return An empty list.
    #'
    get_metadata = function(key = NULL) {
      tryCatch(
        expr = private$.ephemeral_error("custom", "and have no metadata"),
        error = function(e) {
          warning(conditionMessage(e), call. = FALSE, immediate. = TRUE)
        }
      )
      return(list())
    },

    # Override SOMACollectionBase methods
    #' @description Dummy method for ephemeral objects for compatibility with
    #' SOMA collections.
    #'
    #' @param object,key Ignored for ephemeral objects.
    #'
    #' @return Throws an error as this method is not supported by ephemeral
    #' objects.
    #'
    add_new_collection = \(object, key) private$.ephemeral_error(),

    #' @description Dummy method for ephemeral objects for compatibility with
    #' SOMA collections.
    #'
    #' @param key,schema,index_column_names Ignored for ephemeral objects.
    #'
    #' @return Throws an error as this method is not supported by ephemeral
    #' objects.
    #'
    add_new_dataframe = function(key, schema, index_column_names) {
      private$.ephemeral_error()
    },

    #' @description Dummy method for ephemeral objects for compatibility with
    #' SOMA collections.
    #'
    #' @param key,type,shape Ignored for ephemeral objects.
    #'
    #' @return Throws an error as this method is not supported by ephemeral
    #' objects.
    #'
    add_new_dense_ndarray = \(key, type, shape) private$.ephemeral_error(),

    #' @description Dummy method for ephemeral objects for compatibility with
    #' SOMA collections.
    #'
    #' @param key,type,shape Ignored for ephemeral objects.
    #'
    #' @return Throws an error as this method is not supported by ephemeral
    #' objects.
    #'
    add_new_sparse_ndarray = \(key, type, shape) private$.ephemeral_error()
  ),
  active = list(
    #' @field uri \dQuote{\code{ephemeral-collection:<MEMORY_ADDRESS>}}.
    #'
    uri = function(value) {
      if (!missing(value)) {
        private$.read_only_error("uri")
      }
      return(paste0("ephemeral-collection:", data.table::address(self)))
    },

    #' @field members A list with the members of this collection
    #'
    members = function(value) {
      if (!missing(value)) {
        private$.read_only_error("members")
      }
      return(private$.data)
    },

    # Override SOMACollectionBase fields
    #' @field soma_type Dummy field for ephemeral objects for compatibility with
    #' SOMA collections.
    #'
    soma_type = function(value) {
      if (!missing(value)) {
        private$.read_only_error("soma_type")
      }
      private$.ephemeral_error("custom", "and have no SOMA type")
    },

    # Override SOMAObject fields
    #' @field platform_config Dummy field for ephemeral objects for
    #' compatibility with SOMA collections.
    #'
    platform_config = function(value) {
      if (!missing(value)) {
        private$.read_only_error("platform_config")
      }
      private$.ephemeral_error("custom", "and have no configuration")
    },

    #' @field tiledbsoma_ctx Dummy field for ephemeral objects for compatibility
    #' with SOMA collections.
    #'
    tiledbsoma_ctx = function(value) {
      if (!missing(value)) {
        private$.read_only_error("tiledbsoma_ctx")
      }
      private$.ephemeral_error("custom", "and have no context")
    },

    #' @field tiledb_timestamp Dummy field for ephemeral objects for
    #' compatibility with SOMA collections
    #'
    tiledb_timestamp = function(value) {
      if (!missing(value)) {
        private$.read_only_error("tiledb_timestamp")
      }
      private$.ephemeral_error("custom", "and have no timestamp")
    },

    #' @field .tiledb_timestamp_range Dummy field for ephemeral objects for
    #' compatibility with SOMA collections
    #'
    .tiledb_timestamp_range = function(value) {
      if (!missing(value)) {
        private$.read_only_error("tiledb_timestamp_range")
      }
      private$.ephemeral_error("custom", "and have no timestamp")
    }
  ),
  private = list(
    # Override SOMAObject private fields
    .platform_config = NULL,
    .tiledbsoma_ctx = NULL,
    .tiledb_timestamp = NULL,
    .context = NULL,
    .mode = NULL,
    .uri = character(1L),

    # Override SOMAObject private methods
    .check_open_for_read = \() invisible(NULL),
    .check_open_for_write = \() invisible(NULL),
    .check_open_for_delete = \() invisible(NULL),
    .check_open_for_read_or_write = \() invisible(NULL),
    .check_open = \() invisible(NULL),

    # Override SOMACollectionBase private fields
    .member_cache = NULL,
    .update_member_cache = \() invisible(self),

    # Ephemeral fields
    .data = NULL,

    # Ephemeral methods
    .ephemeral_error = function(type = "added", msg = NULL) {
      stopifnot(
        "'type' must be a single character value" = is_scalar_character(type)
      )
      type <- match.arg(
        arg = type,
        choices = c("base", "added", "opened", "edited", "custom")
      )
      if (type == "custom" && !is_scalar_character(msg)) {
        stop("'msg' must be a single character value")
      }
      stop(
        sQuote(self$class()),
        " objects are ephemeral",
        switch(
          EXPR = type,
          added = " and cannot be added to",
          opened = " and cannot be opened",
          edited = " and cannot be edited",
          custom = paste0(" ", trimws(msg))
        ),
        call. = FALSE
      )
    }
  )
)

#' Ephemeral Collections
#'
#' @description Ephemeral version of \code{\link{SOMACollection}s}; ephemeral
#' collections are equivalent to
#' \link[tiledbsoma:SOMACollection]{SOMA collections} but are stored in-memory
#' instead of on-disk.
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' (col <- EphemeralCollection$new())
#' col$soma_type
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' dir <- withr::local_tempfile(pattern = "obs")
#' dir.create(dir, recursive = TRUE)
#'
#' (obs <- load_dataset("soma-dataframe-pbmc3k-processed-obs", dir))
#' col$set(obs, "obs")
#' col$names()
#'
#' \dontshow{
#' obs$close()
#' }
#'
EphemeralCollection <- R6::R6Class(
  classname = "EphemeralCollection",
  inherit = EphemeralCollectionBase,
  active = list(
    #' @field soma_type The SOMA object type.
    #'
    soma_type = function(value) {
      if (!missing(value)) {
        private$.read_only_error("soma_type")
      }
      return("SOMACollection")
    }
  )
)

#' Ephemeral SOMA Measurement
#'
#' @description Ephemeral version of \code{\link{SOMAMeasurement}s}; ephemeral
#' measurements are equivalent to
#' \link[tiledbsoma:SOMAMeasurement]{SOMA measurements} but are stored in-memory
#' instead of on-disk.
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' (ms <- EphemeralMeasurement$new())
#' ms$soma_type
#'
#' ms$set(EphemeralCollection$new(), "X")
#' ms$X
#'
EphemeralMeasurement <- R6::R6Class(
  classname = "EphemeralMeasurement",
  inherit = EphemeralCollectionBase,
  active = list(
    #' @field var A \code{\link{SOMADataFrame}} containing primary annotations
    #' on the variable axis, for variables in this measurement (i.e., annotates
    #' columns of \code{X}). The contents of the \code{soma_joinid} column
    #' define the variable index domain, \code{var_id}. All variables for this
    #' measurement must be defined in this data frame.
    #'
    var = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = "var",
        expected_class = "SOMADataFrame"
      )
    },

    #' @field X A \code{\link{SOMACollection}} of
    #' \code{\link{SOMASparseNDArray}}s; each contain measured feature values
    #' indexed by \code{[obsid, varid]}.
    #'
    X = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = "X",
        expected_class = c("EphemeralCollection", "SOMACollection")
      )
    },

    #' @field obsm A \code{\link{SOMACollection}} of
    #' \code{\link{SOMADenseNDArray}}s containing annotations on the observation
    #' axis. Each array is indexed by \code{obsid} and has the same shape as
    #' \code{obs}.
    #'
    obsm = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = "obsm",
        expected_class = c("EphemeralCollection", "SOMACollection")
      )
    },

    #' @field obsp A \code{\link{SOMACollection}} of
    #' \code{\link{SOMASparseNDArray}}s containing pairwise annotations on the
    #' observation axis and indexed with \code{[obsid_1, obsid_2]}.
    #'
    obsp = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = "obsp",
        expected_class = c("EphemeralCollection", "SOMACollection")
      )
    },

    #' @field varm A \code{\link{SOMACollection}} of
    #' \code{\link{SOMADenseNDArray}}s containing annotations on the variable
    #' axis. Each array is indexed by \code{varid} and has the same shape as
    #' \code{var}.
    #'
    varm = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = "varm",
        expected_class = c("EphemeralCollection", "SOMACollection")
      )
    },

    #' @field varp A \code{\link{SOMACollection}} of
    #' \code{\link{SOMASparseNDArray}}s containing pairwise annotations on the
    #' variable axis and indexed with \code{[varid_1, varid_2]}.
    #'
    varp = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = "varp",
        expected_class = c("EphemeralCollection", "SOMACollection")
      )
    },

    #' @field soma_type The SOMA object type.
    #'
    soma_type = function(value) {
      if (!missing(value)) {
        private$.read_only_error("soma_type")
      }
      return("SOMAMeasurement")
    }
  )
)

#' Ephemeral SOMA Experiment
#'
#' @description Ephemeral version of \code{\link{SOMAExperiment}s}; ephemeral
#' experiments are equivalent to
#' \link[tiledbsoma:SOMAExperiment]{SOMA experiments} but are stored in-memory
#' instead of on-disk.
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' (exp <- EphemeralExperiment$new())
#' exp$soma_type
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' dir <- withr::local_tempfile(pattern = "obs")
#' dir.create(dir, recursive = TRUE)
#'
#' (obs <- load_dataset("soma-dataframe-pbmc3k-processed-obs", dir))
#' exp$set(obs, "obs")
#' exp$obs
#'
#' \dontshow{
#' obs$close()
#' }
#'
EphemeralExperiment <- R6::R6Class(
  classname = "EphemeralExperiment",
  inherit = EphemeralCollectionBase,
  active = list(
    #' @field obs A \code{\link{SOMADataFrame}} containing the annotations on
    #' the observation axis. The contents of the \code{soma_joinid} column
    #' define the observation index domain \code{obs_id}. All observations for
    #' the \code{SOMAExperiment} must be defined in this data frame.
    #'
    obs = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = "obs",
        expected_class = "SOMADataFrame"
      )
    },

    #' @field ms A \code{\link{SOMACollection}} of named
    #' \code{\link{SOMAMeasurement}}s.
    #'
    ms = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = "ms",
        expected_class = c("EphemeralCollection", "SOMACollection")
      )
    },

    #' @field soma_type The SOMA object type.
    #'
    soma_type = function(value) {
      if (!missing(value)) {
        private$.read_only_error("soma_type")
      }
      return("SOMAExperiment")
    }
  )
)
