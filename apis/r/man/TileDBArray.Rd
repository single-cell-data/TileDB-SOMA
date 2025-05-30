% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TileDBArray.R
\name{TileDBArray}
\alias{TileDBArray}
\title{TileDB Array Base Class}
\description{
Virtual base class for representing an individual TileDB array
(lifecycle: maturing).
}
\seealso{
Derived classes: \code{\link{SOMAArray}}
}
\keyword{internal}
\section{Super class}{
\code{\link[tiledbsoma:TileDBObject]{tiledbsoma::TileDBObject}} -> \code{TileDBArray}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{object}}{Access the underlying TileB object directly (either a
\code{\link[tiledb]{tiledb_array}} or \code{\link[tiledb]{tiledb_group}}).}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-TileDBArray-open}{\code{TileDBArray$open()}}
\item \href{#method-TileDBArray-close}{\code{TileDBArray$close()}}
\item \href{#method-TileDBArray-print}{\code{TileDBArray$print()}}
\item \href{#method-TileDBArray-tiledb_array}{\code{TileDBArray$tiledb_array()}}
\item \href{#method-TileDBArray-get_metadata}{\code{TileDBArray$get_metadata()}}
\item \href{#method-TileDBArray-set_metadata}{\code{TileDBArray$set_metadata()}}
\item \href{#method-TileDBArray-schema}{\code{TileDBArray$schema()}}
\item \href{#method-TileDBArray-tiledb_schema}{\code{TileDBArray$tiledb_schema()}}
\item \href{#method-TileDBArray-dimensions}{\code{TileDBArray$dimensions()}}
\item \href{#method-TileDBArray-shape}{\code{TileDBArray$shape()}}
\item \href{#method-TileDBArray-maxshape}{\code{TileDBArray$maxshape()}}
\item \href{#method-TileDBArray-non_empty_domain}{\code{TileDBArray$non_empty_domain()}}
\item \href{#method-TileDBArray-ndim}{\code{TileDBArray$ndim()}}
\item \href{#method-TileDBArray-attributes}{\code{TileDBArray$attributes()}}
\item \href{#method-TileDBArray-dimnames}{\code{TileDBArray$dimnames()}}
\item \href{#method-TileDBArray-attrnames}{\code{TileDBArray$attrnames()}}
\item \href{#method-TileDBArray-colnames}{\code{TileDBArray$colnames()}}
\item \href{#method-TileDBArray-index_column_names}{\code{TileDBArray$index_column_names()}}
\item \href{#method-TileDBArray-clone}{\code{TileDBArray$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="TileDBObject" data-id="class"><a href='../../tiledbsoma/html/TileDBObject.html#method-TileDBObject-class'><code>tiledbsoma::TileDBObject$class()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="TileDBObject" data-id="exists"><a href='../../tiledbsoma/html/TileDBObject.html#method-TileDBObject-exists'><code>tiledbsoma::TileDBObject$exists()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="TileDBObject" data-id="initialize"><a href='../../tiledbsoma/html/TileDBObject.html#method-TileDBObject-initialize'><code>tiledbsoma::TileDBObject$initialize()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="TileDBObject" data-id="is_open"><a href='../../tiledbsoma/html/TileDBObject.html#method-TileDBObject-is_open'><code>tiledbsoma::TileDBObject$is_open()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="TileDBObject" data-id="mode"><a href='../../tiledbsoma/html/TileDBObject.html#method-TileDBObject-mode'><code>tiledbsoma::TileDBObject$mode()</code></a></span></li>
<li><span class="pkg-link" data-pkg="tiledbsoma" data-topic="TileDBObject" data-id="reopen"><a href='../../tiledbsoma/html/TileDBObject.html#method-TileDBObject-reopen'><code>tiledbsoma::TileDBObject$reopen()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-open"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-open}{}}}
\subsection{Method \code{open()}}{
Open the SOMA object for read or write.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$open(mode = c("READ", "WRITE"), internal_use_only = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{mode}}{Mode to open in; defaults to \code{READ}.}

\item{\code{internal_use_only}}{Character value to signal this is a 'permitted'
call, as \code{open()} is considered internal and should not be called
directly.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Returns \code{self}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-close"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-close}{}}}
\subsection{Method \code{close()}}{
Close the SOMA object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$close()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Invisibly returns \code{self}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-print"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-print}{}}}
\subsection{Method \code{print()}}{
Print summary of the array (lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-tiledb_array"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-tiledb_array}{}}}
\subsection{Method \code{tiledb_array()}}{
Return a \code{\link[tiledb]{tiledb_array}} object
(lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$tiledb_array(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Optional arguments to pass to
\code{\link[tiledb:tiledb_array]{tiledb::tiledb_array}()}.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A \code{\link[tiledb]{tiledb_array}} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-get_metadata"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-get_metadata}{}}}
\subsection{Method \code{get_metadata()}}{
Retrieve metadata from the TileDB array
(lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$get_metadata(key = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{key}}{The name of the metadata attribute to retrieve.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list of metadata values.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-set_metadata"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-set_metadata}{}}}
\subsection{Method \code{set_metadata()}}{
Add list of metadata to the specified TileDB array
(lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$set_metadata(metadata)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{metadata}}{Named list of metadata to add.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Invisibly returns \code{NULL}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-schema"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-schema}{}}}
\subsection{Method \code{schema()}}{
Retrieve the array schema as an Arrow schema
(lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$schema()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
An \link[arrow:Schema]{Arrow schema} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-tiledb_schema"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-tiledb_schema}{}}}
\subsection{Method \code{tiledb_schema()}}{
Retrieve the array schema as TileDB schema
(lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$tiledb_schema()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A \code{\link[tiledb]{tiledb_array_schema}} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-dimensions"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-dimensions}{}}}
\subsection{Method \code{dimensions()}}{
Retrieve the array dimensions (lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$dimensions()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A named list of \code{\link[tiledb]{tiledb_dim}} objects.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-shape"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-shape}{}}}
\subsection{Method \code{shape()}}{
Retrieve the shape, i.e. the capacity of each dimension.
Attempted reads and writes outside the \code{shape} will result in a
runtime error: this is the purpose of \code{shape}. This will not
necessarily match the bounds of occupied cells within the array.
Using \code{$resize()}, this may be increased up to the hard limit which
\code{$maxshape()} reports (lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$shape()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A named vector of dimension length
(and the same type as the dimension).
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-maxshape"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-maxshape}{}}}
\subsection{Method \code{maxshape()}}{
Retrieve the hard limit up to which the array may be resized
using the \code{$resize()} method (lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$maxshape()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A named vector of dimension length
(and the same type as the dimension).
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-non_empty_domain"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-non_empty_domain}{}}}
\subsection{Method \code{non_empty_domain()}}{
Returns a named list of minimum/maximum pairs, one per
index column, which are the smallest and largest values written on that
index column.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$non_empty_domain(index1 = FALSE, max_only = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{index1}}{Return the non-empty domain with 1-based indices.}

\item{\code{max_only}}{Return only the max value per dimension, and return
this as a vector. Names are dropped (lifecycle: maturing).}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Named list of minimum/maximum values, or integer vector
of maximum values.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-ndim"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-ndim}{}}}
\subsection{Method \code{ndim()}}{
Retrieve number of dimensions (lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$ndim()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A scalar with the number of dimensions.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-attributes"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-attributes}{}}}
\subsection{Method \code{attributes()}}{
Retrieve the array attributes (lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$attributes()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A list of \code{\link[tiledb]{tiledb_attr}} objects.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-dimnames"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-dimnames}{}}}
\subsection{Method \code{dimnames()}}{
Retrieve dimension names (lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$dimnames()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A character vector with the array's dimension names.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-attrnames"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-attrnames}{}}}
\subsection{Method \code{attrnames()}}{
Retrieve attribute names (lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$attrnames()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A character vector with the array's attribute names.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-colnames"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-colnames}{}}}
\subsection{Method \code{colnames()}}{
Retrieve the names of all columns, including dimensions and
attributes (lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$colnames()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A character vector with the array's column names.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-index_column_names"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-index_column_names}{}}}
\subsection{Method \code{index_column_names()}}{
Retrieve names of index (dimension) columns
(lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$index_column_names()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A character vector with the array index (dimension) names.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBArray-clone"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBArray-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBArray$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
