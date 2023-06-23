#' Coercion methods for SOMA classes

# Set R6 methods formally as S3 to create apply generics
methods::setOldClass("SOMASparseNDArrayRead")
methods::setOldClass("TableReadIter")

#' @importClassesFrom Matrix TsparseMatrix CsparseMatrix RsparseMatrix 
NULL

#' @importFrom arrow as_arrow_table
#' @export
arrow::as_arrow_table


#' @importFrom arrow as_arrow_table
#' @export
as_arrow_table.SOMASparseNDArrayRead <- function(x){
  x$tables()$concat()
}

#' @export
as.data.frame.SOMASparseNDArrayRead <- function(x, ...){ 
   as.data.frame(x$tables()$concat(), ...)
}

# Coerce \link[tiledbsoma]{SOMASparseNDArrayRead} to Matrix::\link[Matrix]{dgTMatrix}
setAs(from = "SOMASparseNDArrayRead", 
      to = "TsparseMatrix", 
      def = function(from) from$sparse_matrix()$concat()
      )

# Coerce \link[tiledbsoma]{SOMASparseNDArrayRead} to Matrix::\link[Matrix]{dgCMatrix}
setAs(from = "SOMASparseNDArrayRead", 
      to = "CsparseMatrix", 
      def = function(from) as(as(from, "TsparseMatrix"), "CsparseMatrix")
      )

# Coerce \link[tiledbsoma]{SOMASparseNDArrayRead} to Matrix::\link[Matrix]{dgRMatrix}
setAs(from = "SOMASparseNDArrayRead", 
      to = "RsparseMatrix", 
      def = function(from) as(as(from, "TsparseMatrix"), "RsparseMatrix")
      )

#' @export
as_arrow_table.TableReadIter <- function(x) x$concat()

#' @export
as.data.frame.TableReadIter <- function(x, row.names = NULL, optional = FALSE, ...){ 
   as.data.frame(x$concat(), row.names = row.names, optional = optional, ...)
}
