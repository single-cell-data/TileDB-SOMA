#' Read SOMA Data From a Given URI
#'
#' This functions access a given SOMA URI and returns a complete data.frame. It does
#' not iterate; if your data is large than the initial read size consider the \code{sr_*}
#' functions.
#'
#' @title soma_array_reader
#' @param uri Character value with URI path to a SOMA data set
#' @param colnames Optional vector of character value with the name of the columns to retrieve
#' @param qc Optional external Pointer object to TileDB Query Condition, defaults to \sQuote{NULL} i.e.
#' no query condition
#' @param dim_points Optional named list with vector of data points to select on the given
#' dimension(s). Each dimension can be one entry in the list.
#' @param dim_ranges Optional named list with two-column matrix where each row select a range
#' for the given dimension. Each dimension can be one entry in the list.
#' @param batch_size Character value with the desired batch size, defaults to \sQuote{auto}
#' @param result_order Character value with the desired result order, defaults to \sQuote{auto}
#' @param loglevel Character value with the desired logging level, defaults to \sQuote{auto}
#' which lets prior setting prevail, any other value is set as new logging level.
#' @param config Optional character vector containing TileDB config.
#' @return A List object with two pointers to Arrow array data and schema is returned
#' @examples
#' \dontrun{
#' uri <- extract_dataset("soma-dataframe-pbmc3k-processed-obs")
#' z <- soma_array_reader(uri)
#' arrow::RecordBatch$import_from_c(z$array_data, z$schema)
#' }
#' @noRd
soma_array_reader <- function(uri, colnames = NULL, qc = NULL, dim_points = NULL, dim_ranges = NULL,
                              batch_size = "auto", result_order = "auto", loglevel = "auto",
                              config = NULL) {

    stopifnot("'uri' must be character" = is.character(uri),
              "'colnames' must be character or NULL" = is_character_or_null(colnames),
              "'qc' must be a query condition object pointer or NULL" =
                  #inherits(qc,"tiledb_query_condition") || is.null(qc),
                  is(qc, "externalptr") || is.null(qc),
              "'dim_points' must be a named list or NULL" =
                  is_named_list(dim_points) || is.null(dim_points),
              "'dim_ranges' must be a named list or NULL" =
                  is_named_list(dim_ranges) || is.null(dim_ranges),
              "'batch_size' must be character" = is.character(batch_size),
              "'result_order' must be character" = is.character(result_order),
              "'loglevel' must be character" = is.character(loglevel),
              "'config' must be character or NULL" = is_character_or_null(config))

    if (!is.null(dim_points)) {
        for (i in seq_along(dim_points)) {
            if (is_arrow_array(dim_points[[i]])) {
                obj <- dim_points[[i]]
                dim_points[[i]] <- obj$as_vector()
            }
        }
    }

    soma_array_reader_impl(uri, colnames, qc, dim_points, dim_ranges, batch_size, result_order,
                           loglevel, config)
}
