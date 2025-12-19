#' SOMA Sparse Nd-Array
#'
#' @description \code{SOMASparseNDArray} is a sparse, N-dimensional array with
#' offset (zero-based) integer indexing on each dimension. The
#' \code{SOMASparseNDArray} has a user-defined schema, which includes:
#' \itemize{
#'  \item \code{type}: a \code{primitive} type, expressed as an Arrow type
#'   (e.g., \code{\link[arrow]{int64}}, \code{\link[arrow]{float32}}, etc),
#'   indicating the type of data contained within the array.
#'  \item \code{shape}: the shape of the array, i.e., number and length of each
#'   dimension. This is a soft limit which can be increased using
#'   \code{$resize()} up to the \code{maxshape}.
#'  \item \code{maxshape}: the hard limit up to which \code{shape} may be
#'   increased using \code{$resize()}.
#' }
#' All dimensions must have a positive, non-zero length.
#'
#' @note In TileDB this is an sparse array with \code{N} \code{int64} dimensions
#' of domain \code{[0, maxInt64)} and a single attribute.
#'
#' @section Duplicate Writes:
#' As duplicate index values are not allowed, index values already present in
#' the object are overwritten and new index values are added
#' (lifecycle: maturing).
#'
#' @export
#'
#' @inherit SOMASparseNDArrayCreate examples
#'
SOMASparseNDArray <- R6::R6Class(
  classname = "SOMASparseNDArray",
  inherit = SOMANDArrayBase,
  public = list(
    #' @description Reads a user-defined slice of the \code{SOMASparseNDArray}.
    #'
    #' @param coords Optional \code{list} of integer vectors, one for each
    #' dimension, with a length equal to the number of values to read. If
    #' \code{NULL}, all values are read. List elements can be named when
    #' specifying a subset of dimensions.
    #' @template param-result-order
    #' @param log_level Optional logging level with default value of
    #' \dQuote{\code{warn}}.
    #'
    #' @return A \link{SOMASparseNDArrayRead}.
    #'
    read = function(coords = NULL, result_order = "auto", log_level = "auto") {
      private$.check_open_for_read()
      result_order <- map_query_layout(match_query_layout(result_order))

      if (!is.null(coords)) {
        coords <- private$.convert_coords(coords)
      }

      sr <- mq_setup(
        uri = self$uri,
        private$.soma_context$handle,
        dim_points = coords,
        result_order = result_order,
        timestamprange = self$.tiledb_timestamp_range,
        loglevel = log_level
      )

      return(SOMASparseNDArrayRead$new(sr, self, coords))
    },

    #' @description Write matrix-like data to the array (lifecycle: maturing).
    #'
    #' @param values Any \code{matrix}-like object coercible to a
    #' \code{\link[Matrix:TsparseMatrix-class]{TsparseMatrix}}. Character
    #' dimension names are ignored because \code{SOMANDArray}s use integer
    #' indexing.
    #' @param bbox A vector of integers describing the upper bounds of each
    #' dimension of \code{values}. Generally should be \code{NULL}.
    #'
    #' @return Invisibly returns \code{self}.
    #'
    write = function(values, bbox = NULL) {
      stopifnot(
        "'values' must be a matrix" = is_matrix(values),
        "'bbox' must contain two entries" = is.null(bbox) ||
          length(bbox) == length(dim(values)),
        "'bbox' must be a vector of two integers or a list with each entry containg two integers" = is.null(
          bbox
        ) ||
          .is_integerish(bbox) ||
          (is.list(bbox) &&
            all(vapply_lgl(bbox, function(x, n) length(x) == 2L)))
      )
      # coerce to a TsparseMatrix, which uses 0-based COO indexing
      values <- as(values, Class = "TsparseMatrix")
      coo <- data.frame(
        i = bit64::as.integer64(values@i),
        j = bit64::as.integer64(values@j),
        x = values@x
      )
      if (!is.null(private$.type)) {
        rt <- r_type_from_arrow_type(private$.type)
        if (rt == "integer" && rlang::is_integerish(coo$x)) {
          coo$x <- as.integer(coo$x)
        }
      }
      dnames <- self$dimnames()
      colnames(coo) <- c(dnames, self$attrnames())
      ranges <- sapply(
        X = dnames,
        FUN = function(x) {
          return(range(coo[[x]]))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      )
      bbox <- bbox %||%
        stats::setNames(
          lapply(
            X = dim(x = values) - 1L,
            FUN = function(x) {
              bit64::as.integer64(c(0L, x))
            }
          ),
          nm = dnames
        )
      if (is.null(names(bbox))) {
        names(bbox) <- dnames
      }
      if (!is_named(bbox, allow_empty = FALSE)) {
        # Determine which indexes of `bbox` are unnamed (incl empty strings "")
        # Python equivalent:
        # bbox = pandas.Series([[0, 99], [0, 299]], index=["soma_dim_0", ""])
        # [i for i, key in enumerate(bbox.keys()) if not len(key)]
        idx <- which(!nzchar(names(bbox)))
        names(bbox)[idx] <- dnames[idx]
      }
      if (!identical(sort(names(bbox)), sort(dnames))) {
        stop("The names of 'bbox' must be the names of the array")
      }
      if (rlang::is_integerish(bbox) || bit64::is.integer64(bbox)) {
        bbox <- sapply(
          X = names(bbox),
          FUN = function(x) {
            return(sort(c(min(ranges[[x]]), bbox[[x]])))
          },
          simplify = FALSE
        )
      }
      for (x in dnames) {
        xrange <- bbox[[x]]
        if (any(is.na(xrange))) {
          stop(
            "Ranges in the bounding box must be finite (offending: ",
            sQuote(x),
            ")",
            call. = FALSE
          )
        }
        if (!(rlang::is_integerish(xrange) || bit64::is.integer64(xrange))) {
          stop(
            "Ranges in the bounding box must be integers (offending: ",
            sQuote(x),
            ")",
            call. = FALSE
          )
        }
        xrange <- sort(bit64::as.integer64(xrange))
        if (length(xrange) != 2L) {
          stop(
            "Ranges in the bounding box must consist of two integerish values",
            call. = FALSE
          )
        }
        if (xrange[1L] < 0 || xrange[1L] > min(ranges[[x]])) {
          stop(
            "Ranges in the bounding box must be greater than zero and less than the lowest value being added (offending: ",
            sQuote(x),
            ")",
            call. = FALSE
          )
        }
        if (xrange[2L] < max(ranges[[x]])) {
          stop(
            "Ranges in the bounding box must be greater than the largest value being added (offending: ",
            sQuote(x),
            ")",
            call. = FALSE
          )
        }
        bbox[[x]] <- xrange
      }
      names(bbox) <- paste0(names(bbox), "_domain")
      bbox_flat <- vector(mode = "list", length = length(x = bbox) * 2L)
      index <- 1L
      for (i in seq_along(bbox)) {
        bbox_flat[[index]] <- bbox[[i]][1L]
        bbox_flat[[index + 1L]] <- bbox[[i]][2L]
        names(bbox_flat)[index:(index + 1L)] <- paste0(
          names(bbox)[i],
          c("_lower", "_upper")
        )
        index <- index + 2L
      }
      self$set_metadata(bbox_flat)
      soma_debug(sprintf(
        "[SOMASparseNDArray$write] Calling .write_coo_df (%s)",
        self$tiledb_timestamp %||% "now"
      ))

      self$.write_coordinates(coo)

      return(invisible(self))
    },

    #' @description Retrieve number of non-zero elements (lifecycle: maturing).
    #'
    #' @return A scalar with the number of non-zero elements.
    #'
    nnz = function() {
      nnz(self$uri, private$.soma_context$handle)
    },

    #' @description Write a COO table to the array.
    #'
    #' @param values A \code{data.frame} or \link[arrow:Table]{Arrow table}
    #' with data in COO format; must be named with the dimension and attribute
    #' labels of the array.
    #'
    #' @return Invisibly returns \code{self}.
    #'
    .write_coordinates = function(values) {
      private$.check_open_for_write()
      dnames <- self$dimnames()
      attrn <- self$attrnames()

      stopifnot(
        "'values' must be a data frame or Arrow Table" = is.data.frame(
          values
        ) ||
          inherits(values, what = "Table"),
        "'values' must have one column for each dimension and the data" = ncol(
          values
        ) ==
          length(dnames) + 1L,
        "'values' must be named with the dimension and attribute labels" = is.null(names(
          values
        )) ||
          identical(names(values), c(dnames, attrn))
      )

      # Arrow Tables cannot have NULL names, so this only applies to dataframes
      if (is.null(names(values))) {
        soma_warn(
          "[SOMASparseNDArray$.write_coordinates] no names on input data frame, assuming <dimensions[...], data> order"
        )
        names(values) <- c(dnames, attrn)
      }

      # Check dimensions
      soma_debug(
        "[SOMASparseNDArray$.write_coordinates] checking dimension values"
      )
      for (i in seq_along(dnames)) {
        dn <- dnames[i]
        offending <- sprintf("(offending column: '%s')", dn)
        if (!.is_integerish(values[[dn]])) {
          stop("All dimension columns must be integerish ", offending)
        }
        if (as.logical(min(values[[dn]]) < 0L)) {
          stop("Dimension columns cannot contain negative values ", offending)
        }
        if (as.logical(max(values[[dn]]) >= as.numeric(self$shape()[i]))) {
          stop(
            "Dimension columns cannot exceed the shape of the array ",
            offending
          )
        }
      }

      # Check attribute
      soma_debug("[SOMASparseNDArray$.write_coordinates] checking data values")
      if (is.null(private$.type)) {
        tt <- self$schema()[attrn]$type
        if (is.null(tt)) {
          tt <- if (is.data.frame(values)) {
            arrow::infer_type(values[[attrn]])
          } else {
            values[[attrn]]$type
          }
        }
        private$.type <- tt
      }
      vt <- if (is.data.frame(values)) {
        arrow::infer_type(values[[attrn]])
      } else {
        values[[attrn]]$type
      }
      if (
        (vrt <- r_type_from_arrow_type(vt)) !=
          (rt <- r_type_from_arrow_type(private$.type))
      ) {
        stop("The data column must be of type '", rt, "', got '", vrt, "'")
      }

      # Build our Arrow table and schema
      fields <- c(
        lapply(dnames, arrow::field, type = arrow::int64()),
        arrow::field(attrn, private$.type)
      )
      sch <- do.call(arrow::schema, fields)
      tbl <- arrow::as_arrow_table(values, schema = sch)

      # Write via libtiledbsoma
      soma_debug("[SOMASparseNDArray$.write_coordinates] writing arrow table")
      naap <- nanoarrow::nanoarrow_allocate_array()
      nasp <- nanoarrow::nanoarrow_allocate_schema()
      arrow::as_record_batch(tbl)$export_to_c(naap, nasp)
      writeArrayFromArrow(
        uri = self$uri,
        naap = naap,
        nasp = nasp,
        ctxxp = private$.soma_context$handle,
        arraytype = "SOMASparseNDArray",
        config = NULL,
        tsvec = self$.tiledb_timestamp_range
      )
      return(invisible(self))
    }
  ),
  private = list(
    # Given a user-specified shape along a particular dimension, returns a named
    # list containing name, capacity, and extent elements. If no shape is
    # provided the .Machine$integer.max - 1 is used.
    .dim_capacity_and_extent = function(name, shape = NULL, create_options) {
      out <- list(name = name, capacity = NULL, extent = NULL)

      if (is.null(shape)) {
        out$capacity <- .Machine$integer.max - 1
        out$extent <- min(out$capacity, create_options$dim_tile(name))
      } else {
        stopifnot(
          "'shape' must be a positive scalar integer" = rlang::is_scalar_integerish(
            shape
          ) &&
            shape > 0
        )
        out$capacity <- shape
        out$extent <- min(shape, create_options$dim_tile(name))
      }

      out
    }
  )
)
