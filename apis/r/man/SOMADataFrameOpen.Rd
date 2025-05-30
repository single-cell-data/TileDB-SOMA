% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Factory.R
\name{SOMADataFrameOpen}
\alias{SOMADataFrameOpen}
\title{Open a SOMA Data Frame}
\usage{
SOMADataFrameOpen(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  soma_context = NULL
)
}
\arguments{
\item{uri}{URI for the TileDB object.}

\item{mode}{One of \dQuote{\code{READ}} or \dQuote{\code{WRITE}}.}

\item{platform_config}{Optional platform configuration.}

\item{tiledbsoma_ctx}{Optional SOMATileDBContext.}

\item{tiledb_timestamp}{Optional Datetime (POSIXct) for TileDB timestamp;
defaults to the current time.}

\item{soma_context}{A SOMA context as created by
\code{\link{soma_context}()}.}
}
\value{
A \link[tiledbsoma:SOMADataFrame]{SOMA data frame} stored at
\code{uri} opened in mode \code{mode}.
}
\description{
Factory function to open a \link[tiledbsoma:SOMADataFrame]{SOMA data frame}
for reading (lifecycle: maturing).
}
\examples{
\dontshow{if (requireNamespace("withr", quietly = TRUE)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
uri <- withr::local_tempfile(pattern = "soma-data-frame")
df <- data.frame(
  soma_joinid = bit64::seq.integer64(0L, 99L),
  group = sample(factor(c("g1", "g2")), size = 100L, replace = TRUE),
  nCount = stats::rbinom(100L, 10L, 0.3)
)
(sch <- arrow::infer_schema(df))
(sdf <- SOMADataFrameCreate(uri, sch, domain = list(soma_joinid = c(0, 100))))
sdf$write(arrow::as_arrow_table(df, schema = sch))
sdf$close()

(sdf <- SOMADataFrameOpen(uri))
head(as.data.frame(sdf$read()$concat()))
\dontshow{
sdf$close()
}
\dontshow{\}) # examplesIf}
}
