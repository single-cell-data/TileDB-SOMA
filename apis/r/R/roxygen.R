rd_return_virtual <- function() {
  return("This is a \\strong{virtual} class and cannot be directly instantiated")
}

rd_atomic <- function() {
  return(paste(
    '\\itemize{',
    paste0('\\item \\dQuote{\\code{', .SCALAR_TYPES(), '}}', collapse = '\n'),
    '}',
    sep = '\n'
  ))
}

rd_ephemeral_desc <- function() {
  return("Dummy method for ephemeral cobjects for compatibility with SOMA collections")
}

rd_ephemeral_error <- function() {
  return("Throws an error as this method is not supported by ephemeral objects")
}

rd_ephemeral_field <- function() {
  return("Dummy field for ephemeral objects for compatibility with SOMA collections")
}

rd_ephemeral_param <- function() {
  return("Ignored for ephemeral objects")
}
