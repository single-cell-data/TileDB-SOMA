#' SOMA Array Base Class
#'
#' Virtual base class to add SOMA-specific functionality to the
#' \code{\link{TileDBArray}} class (lifecycle: maturing).
#'
#' @keywords internal
#'
#' @export
#'
#' @seealso Derived classes: \code{\link{SOMADataFrame}},
#' \code{\link{SOMANDArrayBase}}
#'
SOMAArrayBase <- R6::R6Class(
  classname = "SOMAArrayBase",
  inherit = TileDBArray,
  public = list(

    #' @description Does an array allow duplicates?
    #'
    #' @return \code{TRUE} if the underlying TileDB array allows duplicates;
    #' otherwise \code{FALSE}.
    #'
    allows_duplicates = \() c_allows_dups(self$uri, private$.soma_context),

    #' @description Retrieve the array attributes.
    #'
    #' @return A named list of array attributes; each entry contains the
    #' following named entries:
    #' \itemize{
    #'  \item \dQuote{\code{name}}: name of the attribute.
    #'  \item \dQuote{\code{type}}: datatype of the attribute.
    #'  \item \dQuote{\code{ncells}}: number of values per dimension cell.
    #'  \item \dQuote{\code{nullable}}: is the attribute nullable.
    #'  \item \dQuote{\code{filter_list}}: a list with filter information; this
    #'   list contains the following entries:
    #'   \itemize{
    #'    \item \dQuote{\code{filter_type}}
    #'    \item \dQuote{\code{compression_level}}
    #'    \item \dQuote{\code{bit_width}}
    #'    \item \dQuote{\code{positive_delta}}
    #'    \item \dQuote{\code{float_bytewidth}}
    #'    \item \dQuote{\code{float_factor}}
    #'    \item \dQuote{\code{float_offset}}
    #'   }
    #' }
    #'
    attributes = \() c_attributes(self$uri, private$.soma_context),

    #' @description Is an array sparse?
    #'
    #' @return \code{TRUE} if the underlying TileDB array is sparse;
    #' otherwise \code{FALSE}.
    #'
    is_sparse = \() c_is_sparse(self$uri, private$.soma_context)

  ),
  active = list(
    #' @field soma_type Retrieve the SOMA object type.
    #'
    soma_type = function(value) {
      stopifnot("'soma_type' is a read-only field" = missing(value))
      if (is.null(private$soma_type_cache)) {
        private$update_soma_type_cache()
      }
      private$soma_type_cache
    }
  ),
  private = list(

    # Cache object's SOMA_OBJECT_TYPE_METADATA_KEY
    soma_type_cache = NULL,
    update_soma_type_cache = function() {
      private$soma_type_cache <- self$get_metadata(SOMA_OBJECT_TYPE_METADATA_KEY)
    },
    write_object_type_metadata = function() {
      # private$check_open_for_write()

      meta <- list()
      meta[[SOMA_OBJECT_TYPE_METADATA_KEY]] <- self$class()
      meta[[SOMA_ENCODING_VERSION_METADATA_KEY]] <- SOMA_ENCODING_VERSION
      spdl::debug("[SOMAArrayBase::write_object_metadata] calling set metadata")
      self$set_metadata(meta)
    }
  )
)
