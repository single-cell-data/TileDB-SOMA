% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SOMAMeasurement.R
\name{SOMAMeasurement}
\alias{SOMAMeasurement}
\title{SOMA Measurement}
\description{
A \code{SOMAMeasurement} is a sub-element of a
\code{\link{SOMAExperiment}}, and is otherwise a specialized
\code{\link{SOMACollection}} with pre-defined fields: \code{X}, \code{var},
\code{obsm}/\code{varm}, and \code{obsp}/\code{varp} (see
\emph{Active Bindings} below for details) (lifecycle: maturing).
}
\section{Adding new objects to a collection}{


The \code{\link{SOMAMeasurement}} class provides a number of type-specific methods for
adding new a object to the collection, such as \code{add_new_sparse_ndarray()} and
\code{add_new_dataframe()}. These methods will create the new object and add it as
member of the \code{SOMAMeasurement}. The new object will always inherit the parent
context (see \code{\link{SOMATileDBContext}}) and, by default, its platform
configuration (see \code{\link{PlatformConfig}}). However, the user can override the
default platform configuration by passing a custom configuration to the
\code{platform_config} argument.
}

\examples{
\dontshow{if (requireNamespace("withr", quietly = TRUE)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
uri <- withr::local_tempfile(pattern = "soma-measurement")
var <- data.frame(
  soma_joinid = bit64::seq.integer64(0L, 99L),
  var_id = paste0("feature_", seq_len(100L))
)
sch <- arrow::infer_schema(var)

(ms <- SOMAMeasurementCreate(uri))
sdf <- ms$add_new_dataframe(
  "var",
  sch,
  "soma_joinid",
  list(soma_joinid = c(0, 100))
)
sdf$write(arrow::as_arrow_table(var, schema = sch))
sdf$close()
ms$close()

(ms <- SOMAMeasurementOpen(uri))
ms$var

\dontshow{
ms$close()
}
\dontshow{\}) # examplesIf}
}
\section{Super classes}{
\code{\link[tiledbsoma:SOMAObject]{tiledbsoma::SOMAObject}} -> \code{\link[tiledbsoma:SOMACollectionBase]{tiledbsoma::SOMACollectionBase}} -> \code{SOMAMeasurement}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{var}}{A \code{\link{SOMADataFrame}} containing primary annotations
on the variable axis, for variables in this measurement
(i.e., annotates columns of \code{X}). The contents of the
\code{soma_joinid} column define the variable index domain,
\code{var_id}. All variables for this measurement must be defined in this
data frame.}

\item{\code{X}}{A \code{\link{SOMACollection}} of
\code{\link{SOMASparseNDArray}s}, each contains measured feature values
indexed by \code{[obsid, varid]}.}

\item{\code{obsm}}{A \code{\link{SOMACollection}} of
\code{\link{SOMADenseNDArray}s} containing annotations on the observation
axis. Each array is indexed by \code{obsid} and has the same shape as
\code{obs}.}

\item{\code{obsp}}{A \code{\link{SOMACollection}} of
\code{\link{SOMASparseNDArray}s} containing pairwise annotations on the
observation axis and indexed with \code{[obsid_1, obsid_2]}.}

\item{\code{varm}}{A \code{\link{SOMACollection}} of
\code{\link{SOMADenseNDArray}s} containing annotations on the variable
axis. Each array is indexed by \code{varid} and has the same shape as
\code{var}.}

\item{\code{varp}}{A \code{\link{SOMACollection}} of
\code{\link{SOMASparseNDArray}s} containing pairwise annotations on the
variable axis and indexed with \code{[varid_1, varid_2]}.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-SOMAMeasurement-clone}{\code{SOMAMeasurement$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMAObject" data-id="class"><a href='../../tiledbsoma/html/SOMAObject.html#method-SOMAObject-class'><code>tiledbsoma::SOMAObject$class()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMAObject" data-id="exists"><a href='../../tiledbsoma/html/SOMAObject.html#method-SOMAObject-exists'><code>tiledbsoma::SOMAObject$exists()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMAObject" data-id="get_metadata"><a href='../../tiledbsoma/html/SOMAObject.html#method-SOMAObject-get_metadata'><code>tiledbsoma::SOMAObject$get_metadata()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMAObject" data-id="initialize"><a href='../../tiledbsoma/html/SOMAObject.html#method-SOMAObject-initialize'><code>tiledbsoma::SOMAObject$initialize()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMAObject" data-id="is_open"><a href='../../tiledbsoma/html/SOMAObject.html#method-SOMAObject-is_open'><code>tiledbsoma::SOMAObject$is_open()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMAObject" data-id="mode"><a href='../../tiledbsoma/html/SOMAObject.html#method-SOMAObject-mode'><code>tiledbsoma::SOMAObject$mode()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMAObject" data-id="reopen"><a href='../../tiledbsoma/html/SOMAObject.html#method-SOMAObject-reopen'><code>tiledbsoma::SOMAObject$reopen()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="add_new_collection"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-add_new_collection'><code>tiledbsoma::SOMACollectionBase$add_new_collection()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="add_new_dataframe"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-add_new_dataframe'><code>tiledbsoma::SOMACollectionBase$add_new_dataframe()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="add_new_dense_ndarray"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-add_new_dense_ndarray'><code>tiledbsoma::SOMACollectionBase$add_new_dense_ndarray()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="add_new_sparse_ndarray"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-add_new_sparse_ndarray'><code>tiledbsoma::SOMACollectionBase$add_new_sparse_ndarray()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="close"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-close'><code>tiledbsoma::SOMACollectionBase$close()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="create"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-create'><code>tiledbsoma::SOMACollectionBase$create()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="get"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-get'><code>tiledbsoma::SOMACollectionBase$get()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="length"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-length'><code>tiledbsoma::SOMACollectionBase$length()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="names"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-names'><code>tiledbsoma::SOMACollectionBase$names()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="open"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-open'><code>tiledbsoma::SOMACollectionBase$open()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="print"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-print'><code>tiledbsoma::SOMACollectionBase$print()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="remove"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-remove'><code>tiledbsoma::SOMACollectionBase$remove()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="set"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-set'><code>tiledbsoma::SOMACollectionBase$set()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="SOMACollectionBase" data-id="set_metadata"><a href='../../tiledbsoma/html/SOMACollectionBase.html#method-SOMACollectionBase-set_metadata'><code>tiledbsoma::SOMACollectionBase$set_metadata()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SOMAMeasurement-clone"></a>}}
\if{latex}{\out{\hypertarget{method-SOMAMeasurement-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SOMAMeasurement$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
