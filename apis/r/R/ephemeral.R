#' Ephemeral Collection Base
#'
#' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_cls("collection", TRUE)}
#'
#' @keywords internal
#' @importFrom data.table address
#'
EphemeralCollectionBase <- R6::R6Class(
  classname = 'EphemeralCollectionBase',
  inherit = SOMACollectionBase,

  public = list(
    # Override TileDBObject methods
    #' @description Create an ephemeral collection
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
      if (rlang::dots_n(...)) {
        tryCatch(
          expr = private$.ephemeral_error('custom', 'and cannot be customized'),
          error = function(e) {
            warning(conditionMessage(e), call. = FALSE, immediate. = TRUE)
          }
        )
      }
      private$.data <- list()
    },

    #' @description Create a new, empty ephemeral collection
    #'
    #' @return Returns a new ephemeral collection of class \code{class(self)}
    #'
    create = function() {
      gen <- getAnywhere(self$class())[['objs']][[1L]]
      if (!R6::is.R6Class(gen)) {
        stop(
          "Cannot find the class generator for ",
          sQuote(self$class()),
          call. = FALSE
        )
      }
      return (gen$new())
    },

    # Override TileDBGroup private methods
    #' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_desc()}
    #'
    #' @param mode \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_param()}
    #'
    #' @return \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_error()}
    #'
    open = function(mode) {
      private$.ephemeral_error('opened')
    },

    #' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_desc()}
    #'
    #' @return Invisibly returns \code{NULL}
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

    #' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_desc()}
    #'
    #' @return Returns \code{FALSE} as ephemeral collections do not
    #' exist on-disk
    #'
    exists = function() {
      return(FALSE)
    },

    # Override TileDBGroup methods
    #' @description Special method for printing object representation to console
    #'
    #' @return Prints details about the ephemeral collection and invisibly
    #' returns itself
    #'
    print = function() {
      super$print()
      private$format_members()
      return(invisible(self))
    },

    #' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_desc()}
    #'
    #' @param param \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_param()}
    #'
    #' @return Returns \code{NULL} as ephemeral collections do not have an
    #' on-disk configuration
    #'
    get_tiledb_config = function(param = NULL) {
      if (!is.null(param)) {
        tryCatch(
          expr = private$.ephemeral_error('custom', 'and have no TileDB configuration'),
          error = function(e) {
            warning(conditionMessage(e), call. = FALSE, immediate. = TRUE)
          }
        )
      }
      return(NULL)
    },

    #' @description Retrieve the number of items in the collection
    #'
    #' @return The length of the collection
    length = function() {
      length(private$.data)
    },

    #' @description Retrieve the names of members. (lifecycle: experimental)
    #' @return A \code{character} vector of member names.
    names = function() {
      names(private$.data) %||% character(length = 0L)
    },

    #' @description Add object to an ephemeral collection
    #'
    #' @param object A TileDB object (eg. \code{\link{TileDBGroup}}) to add
    #' to the collection
    #' @param name A name to add \code{object} as
    #' @param relative \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_param()}
    #'
    #' @return \[chainable] Invisibly returns \code{self} with \code{object}
    #' added as \code{name}
    #'
    set = function(object, name = NULL, relative = NULL) {
      stopifnot(
        "Only 'TileDBArray' or 'TileDBGroup' objects can be added" =
          inherits(object, "TileDBGroup") || inherits(object, "TileDBArray"),
        is.null(name) || is_scalar_character(name),
        is.null(relative) || is_scalar_logical(relative)
      )
      if (!is.null(relative)) {
        tryCatch(
          expr = private$.ephemeral_error('custom', 'so relative has no effect'),
          error = function(e) {
            warning(conditionMessage(e), call. = FALSE, immediate. = TRUE)
          }
        )
      }
      name <- name %||% object$uri
      private$.data[[name]] <- object
      return(invisible(self))
    },

    #' @description Get objects from an ephemeral collection
    #'
    #' @param name Name of object in the collection to get
    #'
    #' @return The object named \code{name}
    #'
    get = function(name) {
      stopifnot(is_scalar_character(name))
      name <- match.arg(arg = name, choices = self$names())
      return(private$.data[[name]])
    },

    #' @description Remove objects from an ephemeral collection
    #'
    #' @param name Name of object to remove from the collection
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with the object at
    #' \code{name} removed
    remove = function(name) {
      stopifnot(is_scalar_character(name))
      name <- match.arg(arg = name, choices = self$names())
      private$.data[[name]] <- NULL
      return(invisible(self))
    },

    #' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_desc()}
    #'
    #' @param metadata \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_param()}
    #'
    #' @return \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_error()}
    #'
    set_metadata = function(metadata) {
      private$.ephemeral_error('edited')
    },

    #' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_desc()}
    #'
    #' @param key \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_param()}
    #'
    #' @return An empty list
    #'
    get_metadata = function(key = NULL) {
      tryCatch(
        expr = private$.ephemeral_error('custom', 'and have no metadata'),
        error = function(e) {
          warning(conditionMessage(e), call. = FALSE, immediate. = TRUE)
        }
      )
      return(list())
    },

    # Override SOMACollectionBase methods
    #' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_desc()}
    #'
    #' @param object,key \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_param()}
    #'
    #' @return \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_error()}
    #'
    add_new_collection = function(object, key) {
      private$.ephemeral_error()
    },

    #' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_desc()}
    #'
    #' @param key,schema,index_column_names \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_param()}
    #'
    #' @return \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_error()}
    #'
    add_new_dataframe = function(key, schema, index_column_names) {
      private$.ephemeral_error()
    },

    #' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_desc()}
    #'
    #' @param key,type,shape \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_param()}
    #'
    #' @return \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_error()}
    #'
    add_new_dense_ndarray = function(key, type, shape) {
      private$.ephemeral_error()
    },

    #' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_desc()}
    #'
    #' @param key,type,shape \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_param()}
    #'
    #' @return \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_error()}
    #'
    add_new_sparse_ndarray = function(key, type, shape) {
      private$.ephemeral_error()
    }
  ),

  active = list(

    #' @field uri \dQuote{\code{ephemeral-collection:<MEMORY_ADDRESS>}}
    uri = function(value) {
      if (!missing(value)) {
        private$.read_only_error('uri')
      }
      return(paste0('ephemeral-collection:', data.table::address(self)))
    },

    # Override SOMACollectionBase fields
    #' @field soma_type \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_field()}
    soma_type = function(value) {
      if (!missing(value)) {
        private$.read_only_error('soma_type')
      }
      private$.ephemeral_error('custom', 'and have no SOMA type')
    },

    # Override TileDBObject fields
    #' @field platform_config \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_field()}
    platform_config = function(value) {
      if (!missing(value)) {
        private$.read_only_error('platform_config')
      }
      private$.ephemeral_error('custom', 'and have no configuration')
    },

    #' @field tiledbsoma_ctx \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_field()}
    tiledbsoma_ctx = function(value) {
      if (!missing(value)) {
        private$.read_only_error('tiledbsoma_ctx')
      }
      private$.ephemeral_error('custom', 'and have no context')
    },

    #' @field object \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_field()}
    object = function(value) {
      if (!missing(value)) {
        private$.read_only_error('object')
      }
      private$.ephemeral_error('custom', 'and have no underlying object')
    }
  ),

  private = list(
    # Override SOMACollectionBase private fields
    tiledb_uri = NULL,
    tiledb_platform_config = NULL,
    .tiledbsoma_ctx = NULL,
    mode = NULL,

    # Override TileDBGroup private fields
    member_cache = NULL,

    check_open_for_read = function() { },
    check_open_for_write = function() { },
    check_open_for_read_or_write = function() { },

    fill_member_cache_if_null = function() { },
    update_member_cache = function() { },

    # Override SOMACollectionBase private fields
    soma_type_cache = NULL,

    initialize_object = function() {
      private$.ephemeral_error('custom', 'and cannot be initialized')
    },

    get_all_members = function() {
      if (!length(private$.data)) {
        return(list())
      }
      members <- vector(mode = 'list', length = length(private$.data))
      names(members) <- names(private$.data)
      for (i in seq_along(members)) {
        members[[i]] <- list(
          type = tiledb::tiledb_object_type(private$.data[[i]]$uri),
          uri = private$.data[[i]]$uri,
          name = names(private$.data)[i]
        )
      }
      return(members)
    },

    # Ephemeral fields
    .data = NULL,

    # Ephemeral methods
    .ephemeral_error = function(type = 'added', msg = NULL) {
      stopifnot("'type' must be a single character value" = is_scalar_character(type))
      type <- match.arg(
        arg = type,
        choices = c(
          'base',
          'added',
          'opened',
          'edited',
          'custom'
        )
      )
      if (type == 'custom' && !is_scalar_character(msg)) {
        stop("'msg' must be a single character value")
      }
      stop(
        sQuote(self$class()),
        " objects are ephemeral",
        switch(
          EXPR = type,
          added = ' and cannot be added to',
          opened = ' and cannot be opened',
          edited = ' and cannot be edited',
          custom = paste0(' ', trimws(msg))
        ),
        call. = FALSE
      )
    }
  )
)

#' Ephemeral Collections
#'
#' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_cls("collection")}
#'
#' @keywords internal
#'
#' @export
#'
EphemeralCollection <- R6::R6Class(
  classname = 'EphemeralCollection',
  inherit = EphemeralCollectionBase,

  active = list(
    #' @field soma_type The SOMA object type
    soma_type = function(value) {
      if (!missing(value)) {
        private$.read_only_error('soma_type')
      }
      return('SOMACollection')
    }
  )
)

#' Ephemeral SOMA Measurement
#'
#' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_cls("measurement")}
#'
#' @keywords internal
#'
#' @export
#'
EphemeralMeasurement <- R6::R6Class(
  classname = 'EphemeralMeasurement',
  inherit = EphemeralCollectionBase,

  active = list(
    #' @field var \Sexpr[results=rd]{tiledbsoma:::rd_soma_field("var")}
    var = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = 'var',
        expected_class = 'SOMADataFrame'
      )
    },

    #' @field X \Sexpr[results=rd]{tiledbsoma:::rd_soma_field("X")}
    X = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = 'X',
        expected_class = c('EphemeralCollection', 'SOMACollection')
      )
    },

    #' @field obsm \Sexpr[results=rd]{tiledbsoma:::rd_soma_field("obsm")}
    obsm = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = 'obsm',
        expected_class = c('EphemeralCollection', 'SOMACollection')
      )
    },

    #' @field obsp \Sexpr[results=rd]{tiledbsoma:::rd_soma_field("obsp")}
    obsp = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = 'obsp',
        expected_class = c('EphemeralCollection', 'SOMACollection')
      )
    },

    #' @field varm \Sexpr[results=rd]{tiledbsoma:::rd_soma_field("varm")}
    varm = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = 'varm',
        expected_class = c('EphemeralCollection', 'SOMACollection')
      )
    },

    #' @field varp \Sexpr[results=rd]{tiledbsoma:::rd_soma_field("varp")}
    varp = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = 'varp',
        expected_class = c('EphemeralCollection', 'SOMACollection')
      )
    },

    #' @field soma_type The SOMA object type
    soma_type = function(value) {
      if (!missing(value)) {
        private$.read_only_error('soma_type')
      }
      return('SOMAMeasurement')
    }
  )
)

#' Ephemeral SOMA Experiment
#'
#' @description \Sexpr[results=rd]{tiledbsoma:::rd_ephemeral_cls("experiment")}
#'
#' @keywords internal
#'
#' @export
#'
EphemeralExperiment <- R6::R6Class(
  classname = 'EphemeralExperiment',
  inherit = EphemeralCollectionBase,

  active = list(
    #' @field obs \Sexpr[results=rd]{tiledbsoma:::rd_soma_field("obs")}
    obs = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = 'obs',
        expected_class = 'SOMADataFrame'
      )
    },

    #' @field ms \Sexpr[results=rd]{tiledbsoma:::rd_soma_field("ms")}
    ms = function(value) {
      private$get_or_set_soma_field(
        value = value,
        name = 'ms',
        expected_class = c('EphemeralCollection', 'SOMACollection')
      )
    },

    #' @field soma_type The SOMA object type
    soma_type = function(value) {
      if (!missing(value)) {
        private$.read_only_error('soma_type')
      }
      return('SOMAExperiment')
    }
  )
)
