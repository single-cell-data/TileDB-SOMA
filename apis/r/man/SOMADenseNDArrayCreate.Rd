% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Factory.R
\name{SOMADenseNDArrayCreate}
\alias{SOMADenseNDArrayCreate}
\title{Create a SOMA Dense ND Array}
\usage{
SOMADenseNDArrayCreate(
  uri,
  type,
  shape,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
)
}
\arguments{
\item{uri}{URI for the TileDB object.}

\item{type}{An \link[arrow:data-type]{Arrow type} defining the type
of each element in the array.}

\item{shape}{A vector of integers defining the shape of the array.}

\item{platform_config}{Optional platform configuration.}

\item{tiledbsoma_ctx}{Optional SOMATileDBContext.}

\item{tiledb_timestamp}{Optional Datetime (POSIXct) for TileDB timestamp.}
}
\value{
A new \link[tiledbsoma:SOMADenseNDArray]{SOMA dense ND array}
stored at \code{uri} opened for writing.
}
\description{
Factory function to create a
\link[tiledbsoma:SOMADenseNDArray]{SOMA dense ND array} for writing
(lifecycle: maturing).
}
\examples{
\dontshow{if (requireNamespace("withr", quietly = TRUE)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
uri <- withr::local_tempfile(pattern = "soma-dense-array")
mat <- matrix(stats::rnorm(100L ^ 2L), nrow = 100L, ncol = 100L)
mat[1:3, 1:5]

(arr <- SOMADenseNDArrayCreate(uri, arrow::float64(), shape = dim(mat)))
arr$write(mat)
arr$close()

(arr <- SOMADenseNDArrayOpen(uri))
arr$read_arrow_table()
\dontshow{
arr$close()
}
\dontshow{\}) # examplesIf}
}
