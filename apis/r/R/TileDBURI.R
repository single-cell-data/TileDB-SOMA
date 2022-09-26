#' TileDB URI Parser
#' @importFrom urltools url_parse url_compose
#' @noRd
TileDBURI <- R6::R6Class(
  classname = "TileDBURI",
  public = list(

    initialize = function(uri) {
      stopifnot(is_scalar_character(uri))
      private$.uri <- uri
      private$.object_uri <- uri
      private$.pieces <- private$url_parse(uri)

      if (self$is_tiledb_cloud_uri()) {
        private$.tiledb_cloud_uri <- private$url_compose(
          scheme = private$.pieces$scheme,
          domain = private$.pieces$domain,
          path = basename(private$.pieces$path)
        )

        if (self$is_tiledb_cloud_creation_uri()) {
          private$.object_uri <- private$.pieces$path
        }
      }
    },

    #' @description Does the URI point to a remote object?
    is_remote_uri = function() {
      private$.pieces$scheme %||% "" %in% c("tiledb", "s3")
    },

    #' @description Is the URI prefix tiledb://?
    is_tiledb_cloud_uri = function() {
      private$.pieces$scheme %||% "" == "tiledb"
    },

    #' @description Is this a TileDB Cloud creation URI?
    #' Example: "tiledb://namespace/s3://bucket"
    is_tiledb_cloud_creation_uri = function() {
      self$is_tiledb_cloud_uri() && startsWith(private$.pieces$path, "s3://")
    }
  ),

  active = list(
    #' @field uri The original URI
    uri = function(value) {
      if (missing(value)) return(private$.uri)
      self$initialize(value)
    },

    #' @field tiledb_cloud_uri URI for TileDB Cloud access
    tiledb_cloud_uri = function(value) {
      if (missing(value)) return(private$.tiledb_cloud_uri)
      read_only_error("tiledb_cloud_uri")
    },

    #' @field object_uri URI for direct access to the object
    object_uri = function(value) {
      if (missing(value)) return(private$.object_uri)
      read_only_error("object_uri")
    }

  ),

  private = list(
    .uri = character(),
    # Named list of individual URL pieces
    .pieces = list(),
    .tiledb_cloud_uri = NULL,
    .object_uri = NULL,

    # Parse URL and remove empty pieces rather than keeping them as NA
    url_parse = function(url) {
      Filter(
        Negate(is.na),
        as.list(urltools::url_parse(url))
      )
    },

    # Compose URLs using only a subset of pieces
    url_compose = function(scheme = NULL, domain = NULL, path = NULL) {
      url <- urltools::url_compose(
        data.frame(
          scheme = scheme %||% NA_character_,
          domain = domain %||% NA_character_,
          path = path %||% NA_character_,
          port = NA_character_,
          parameter = NA_character_,
          fragment = NA_character_
        )
      )
      # remove trailing slash
      sub("\\/$", "", url)
    },

    read_only_error = function(field_name) {
      stop(glue::glue("'{field_name}' is a read-only field."))
    }
  )
)
