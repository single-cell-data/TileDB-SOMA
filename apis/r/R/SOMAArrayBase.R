#' SOMA Array Base Class
#'
#' Virtual base class to add array-specific functionality to the
#' \code{\link{SOMAObject}} class (lifecycle: maturing).
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
  inherit = SOMAObject,
  public = list(
    #' @description Open the SOMA object for read or write.\cr
    #' \cr
    #' \strong{Note}: \code{open()} is considered internal and should not be
    #' called directly; use factory functions
    #' (eg. \code{\link{SOMASparseNDArrayOpen}()}) instead.
    #'
    #' @param mode Mode to open the object in.
    #'
    #' @return Return s\code{self}.
    #'
    open = function(mode = c("READ", "WRITE", "DELETE")) {
      envs <- unique(vapply(
        X = unique(sys.parents()),
        FUN = function(n) environmentName(environment(sys.function(n))),
        FUN.VALUE = character(1L)
      ))
      if (sys.parent()) {
        if (
          inherits(
            environment(sys.function(sys.parent()))$self,
            what = "SOMAObject"
          )
        ) {
          envs <- union(envs, "tiledbsoma")
        }
      }
      if (!"tiledbsoma" %in% envs) {
        stop(
          paste(
            strwrap(private$.internal_use_only("open", "collection")),
            collapse = '\n'
          ),
          call. = FALSE
        )
      }

      # Set the mode of the array
      private$.mode <- match.arg(mode)

      if (is.null(self$tiledb_timestamp)) {
        soma_debug(sprintf(
          "[SOMAArrayBase$open] Opening %s '%s' in %s mode",
          self$class(),
          self$uri,
          self$mode()
        ))
      } else {
        soma_debug(sprintf(
          "[SOMAArrayBase$open] Opening %s '%s' in %s mode at (%s)",
          self$class(),
          self$uri,
          mode,
          self$tiledb_timestamp %||% "now"
        ))
      }

      private$.update_metadata_cache(TRUE)

      return(self)
    },

    #' @description Close the SOMA array.
    #'
    #' @return Invisibly returns \code{self}.
    #'
    close = function() {
      soma_debug(sprintf(
        "[SOMAArrayBase$close] Closing %s '%s'",
        self$class(),
        self$uri
      ))
      private$.mode <- NULL

      return(invisible(self))
    },

    #' @description Does an array allow duplicates?
    #'
    #' @return \code{TRUE} if the underlying TileDB array allows duplicates;
    #' otherwise \code{FALSE}.
    #'
    allows_duplicates = \() c_allows_dups(self$uri, private$.soma_context$handle),

    #' @description Is an array sparse?
    #'
    #' @return \code{TRUE} if the underlying TileDB array is sparse;
    #' otherwise \code{FALSE}.
    #'
    is_sparse = \() c_is_sparse(self$uri, private$.soma_context$handle),

    #' @description Retrieve the array schema as an Arrow schema
    #' (lifecycle: maturing).
    #'
    #' @return An Arrow \code{\link[arrow:Schema]{Schema}} object.
    #'
    schema = \() arrow::as_schema(c_schema(self$uri, private$.soma_context$handle)),

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
    attributes = \() c_attributes(self$uri, private$.soma_context$handle),

    #' @description Retrieve attribute names (lifecycle: maturing).
    #'
    #' @return A character vector with the array's attribute names.
    #'
    attrnames = \() c_attrnames(self$uri, private$.soma_context$handle),

    #' @description Retrieve the array dimensions (lifecycle: maturing)
    #'
    #' @return A named list of array dimensions; each entry contains the
    #' following named entries:
    #' \itemize{
    #'  \item \dQuote{\code{name}}: name of the dimension.
    #'  \item \dQuote{\code{type}}: datatype of the dimension.
    #'  \item \dQuote{\code{ncells}}: number of values per dimension cell.
    #'  \item \dQuote{\code{domain}}: domain of the dimension.
    #'  \item \dQuote{\code{tile}}: tile of the dimension.
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
    dimensions = \() c_domain(self$uri, private$.soma_context$handle),

    #' @description Retrieve dimension names (lifecycle: maturing).
    #'
    #' @return A character vector with the array's dimension names.
    #'
    dimnames = \() c_dimnames(self$uri, private$.soma_context$handle),

    #' @description Retrieve the names of all columns, including dimensions and
    #' attributes (lifecycle: maturing).
    #'
    #' @return A character vector with the array's column names.
    #'
    colnames = \() c(self$dimnames(), self$attrnames()),

    #' @description Retrieve names of index (dimension) columns (lifecycle: maturing)
    #'
    #' @return A character vector with the array index (dimension) names
    #'
    index_column_names = \() self$dimnames(),

    #' @description Retrieve the shape, i.e. the capacity of each dimension
    #' Attempted reads and writes outside the \code{shape} will result in a
    #' run-time error: this is the purpose of \code{shape}. This will not
    #' necessarily match the bounds of occupied cells within the array.
    #' Using \code{$resize()}, this may be increased up to the hard limit which
    #' \code{maxshape()} reports (lifecycle: maturing).
    #'
    #' @return A named vector of dimension length and of the same type as
    #' the dimension.
    #'
    shape = \() bit64::as.integer64(shape(self$uri, private$.soma_context$handle)),

    #' @description Retrieve the hard limit up to which the array may be resized
    #' using the \code{$resize()} method (lifecycle: maturing).
    #'
    #' @return A named vector of dimension length and of the same type as
    #' the dimension.
    #'
    maxshape = \() {
      bit64::as.integer64(maxshape(self$uri, private$.soma_context$handle))
    },

    #' @description Returns a named list of minimum/maximum pairs, one per index
    #' column, which are the smallest and largest values written on that
    #' index column.
    #'
    #' @param index1 Return the non-empty domain with 1-based indices
    #' @param max_only Return only the max value per dimension, and return
    #' this as a vector. Names are dropped (lifecycle: maturing).
    #'
    #' @return Named list of minimum/maximum values, or integer vector
    #' of maximum values.
    #'
    non_empty_domain = function(index1 = FALSE, max_only = FALSE) {
      retval <- as.list(
        arrow::as_record_batch(
          arrow::as_arrow_table(
            non_empty_domain(self$uri, private$.soma_context$handle)
          )
        )
      )
      if (index1) {
        retval <- lapply(retval, FUN = \(c) c + 1)
      }
      if (max_only) {
        # # No vapply options since SOMADataFrame can have varying types.
        retval <- unname(sapply(retval, FUN = max, USE.NAMES = FALSE))
      }
      return(retval)
    },

    #' @description Retrieve number of dimensions (lifecycle: maturing).
    #'
    #' @return A scalar with the number of dimensions.
    #'
    ndim = \() ndim(self$uri, private$.soma_context$handle),

    #' @description Print-friendly representation of the object.
    #'
    #' @return Invisibly returns \code{self}.
    #'
    print = function() {
      super$print()
      if (self$exists()) {
        cat("  dimensions:", string_collapse(self$dimnames()), "\n")
        cat("  attributes:", string_collapse(self$attrnames()), "\n")
      }
      return(invisible(self))
    }
  ),
  private = list(
    write_object_type_metadata = function() {
      # private$.check_open_for_write()

      meta <- list()
      meta[[SOMA_OBJECT_TYPE_METADATA_KEY]] <- self$class()
      meta[[SOMA_ENCODING_VERSION_METADATA_KEY]] <- SOMA_ENCODING_VERSION
      soma_debug("[SOMAArrayBase::write_object_metadata] calling set metadata")
      self$set_metadata(meta)
    }
  )
)
