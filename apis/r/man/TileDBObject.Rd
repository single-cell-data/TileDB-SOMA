% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TileDBObject.R
\name{TileDBObject}
\alias{TileDBObject}
\title{TileDB Object Base Class}
\description{
Virtual base class to implement shared functionality across the
TileDBArray and TileDBGroup classes (lifecycle: maturing).
}
\seealso{
Derived classes: \code{\link{TileDBArray}},
\code{\link{TileDBGroup}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{platform_config}}{Platform configuration.}

\item{\code{tiledbsoma_ctx}}{SOMATileDBContext.}

\item{\code{tiledb_timestamp}}{Time that object was opened at.}

\item{\code{uri}}{The URI of the TileDB object.}

\item{\code{.tiledb_timestamp_range}}{Time range for libtiledbsoma.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-TileDBObject-new}{\code{TileDBObject$new()}}
\item \href{#method-TileDBObject-class}{\code{TileDBObject$class()}}
\item \href{#method-TileDBObject-is_open}{\code{TileDBObject$is_open()}}
\item \href{#method-TileDBObject-mode}{\code{TileDBObject$mode()}}
\item \href{#method-TileDBObject-reopen}{\code{TileDBObject$reopen()}}
\item \href{#method-TileDBObject-print}{\code{TileDBObject$print()}}
\item \href{#method-TileDBObject-exists}{\code{TileDBObject$exists()}}
\item \href{#method-TileDBObject-clone}{\code{TileDBObject$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBObject-new"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBObject-new}{}}}
\subsection{Method \code{new()}}{
Create a new TileDB object (lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBObject$new(
  uri,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  internal_use_only = NULL,
  soma_context = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{uri}}{URI for the TileDB object.}

\item{\code{platform_config}}{Optional platform configuration.}

\item{\code{tiledbsoma_ctx}}{Optional SOMATileDBContext.}

\item{\code{tiledb_timestamp}}{Optional Datetime (POSIXct) with TileDB
timestamp.}

\item{\code{internal_use_only}}{Character value to signal this is a 'permitted'
call, as \code{new()} is considered internal and should not be called
directly.}

\item{\code{soma_context}}{A SOMA context as created by
\code{\link{soma_context}()}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBObject-class"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBObject-class}{}}}
\subsection{Method \code{class()}}{
Print the name of the R6 class.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBObject$class()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBObject-is_open"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBObject-is_open}{}}}
\subsection{Method \code{is_open()}}{
Determine if the object is open for reading or writing.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBObject$is_open()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
\code{TRUE} if the object is open, otherwise \code{FALSE}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBObject-mode"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBObject-mode}{}}}
\subsection{Method \code{mode()}}{
Get the mode of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBObject$mode()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
If the object is closed, returns \dQuote{\code{CLOSED}};
otherwise returns the mode (eg. \dQuote{\code{READ}}) of the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBObject-reopen"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBObject-reopen}{}}}
\subsection{Method \code{reopen()}}{
Close and reopen the TileDB object in a new mode.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBObject$reopen(mode, tiledb_timestamp = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{mode}}{New mode to open the object in; choose from:
\itemize{
\item \dQuote{\code{READ}}
\item \dQuote{\code{WRITE}}
}}

\item{\code{tiledb_timestamp}}{Optional Datetime (POSIXct) with TileDB
timestamp.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Invisibly returns \code{self} opened in \code{mode}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBObject-print"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBObject-print}{}}}
\subsection{Method \code{print()}}{
Print-friendly representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBObject$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBObject-exists"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBObject-exists}{}}}
\subsection{Method \code{exists()}}{
Check if the object exists (lifecycle: maturing).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBObject$exists()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
\code{TRUE} if the object exists, otherwise \code{FALSE}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TileDBObject-clone"></a>}}
\if{latex}{\out{\hypertarget{method-TileDBObject-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TileDBObject$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
