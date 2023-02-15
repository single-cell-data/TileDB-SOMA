
#' @export
#'
rd_return_virtual <- function() {
  return("This is a \\strong{virtual} class and cannot be directly instantiated")
}

#' @export
#'
rd_atomic <- function() {
  return(paste(
    '\\itemize{',
    paste0('\\item \\dQuote{\\code{', .SCALAR_TYPES(), '}}', collapse = '\n'),
    '}',
    sep = '\n'
  ))
}
