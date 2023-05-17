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

rd_ephemeral_cls <- function(cls, base = FALSE) {
  stopifnot(is_scalar_character(cls), is_scalar_logical(base))
  cls <- match.arg(
    arg = cls,
    choices = c('collection', 'experiment', 'measurement')
  )
  switch(
    EXPR = cls,
    collection = {
      first <- ifelse(
        test = isTRUE(base),
        yes = 'Base class for ephemeral collections',
        no = 'Ephemeral version of \\code{\\link{SOMACollection}s}'
      )
      link <- '\\link[tiledbsoma:SOMACollection]{SOMA collections}'
    },
    experiment = {
      first <- ifelse(
        test = isTRUE(base),
        yes = 'Base class for ephemeral experiments',
        no = 'Ephemeral version of \\code{\\link{SOMAExperiment}s}'
      )
      link <- '\\link[tiledbsoma:SOMAExperiment]{SOMA experiments}'
    },
    measurement = {
      first <- ifelse(
        test = isTRUE(base),
        yes = 'Base class for ephemeral measurements',
        no = 'Ephemeral version of \\code{\\link{SOMAMeasurement}s}'
      )
      link <- '\\link[tiledbsoma:SOMAMeasurement]{SOMA measurements}'
    },
  )
  return(paste0(
    first,
    '; ephemeral ',
    cls,
    's are equivalent to ',
    link,
    ' but are stored in-memory instead of on-disk'
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

rd_soma_field <- function(field) {
  stopifnot(is_scalar_character(field))
  field <- match.arg(
    arg = field,
    choices = c(
      'X',
      'ms',
      'obs',
      'obsm',
      'obsp',
      'var',
      'varm',
      'varp'
    )
  )
  cd <- function(x) {
    return(paste0('\\code{', x, '}'))
  }
  cl <- function(x) {
    return(cd(paste0('\\link{', x, '}')))
  }
  return(switch(
    EXPR = field,
    X = paste(
      'A',
      cl('SOMACollection'),
      'of',
      paste0(cl('SOMASparseNDArray'), 's;'),
      'each contain measured feature values indexed by',
      cd('[obsid, varid]')
    ),
    ms = paste(
      'A',
      cl('SOMACollection'),
      'of named',
      paste0(cl('SOMAMeasurement'), 's')
    ),
    obs = paste(
      'A',
      cl('SOMADataFrame'),
      'containing the annotations on the observation axis.',
      'The contents of the',
      cd('soma_joinid'),
      'column define the observation index domain',
      paste0(cd('obs_id'), '.'),
      'All observations for the',
      cd('SOMAExperiment'),
      'must be defined in this data frame'
    ),
    obsm = paste(
      'A',
      cl('SOMACollection'),
      'of',
      paste0(cl('SOMADenseNDArray'), 's'),
      'containing annotations on the observation axis. Each array is indexed by',
      cd('obsid'),
      'and has the same shape as',
      cd('obs')
    ),
    obsp = paste(
      'A',
      cl('SOMACollection'),
      'of',
      paste0(cl('SOMASparseNDArray'), 's'),
      'containing pairwise annotations on the observation axis and indexed with',
      cd('[obsid_1, obsid_2]')
    ),
    var = paste(
      'A',
      cl('SOMADataFrame'),
      'containing primary annotations on the variable axis,',
      'for variables in this measurement (i.e., annotates columns of',
      paste0(cd('X'), ').'),
      'The contents of the',
      cd('soma_joinid'),
      'column define the variable index domain,',
      paste0(cd('var_id'), '.'),
      'All variables for this measurement must be defined in this data frame'
    ),
    varm = paste(
      'A',
      cl('SOMACollection'),
      'of',
      paste0(cl('SOMADenseNDArray'), 's'),
      'containing annotations on the variable axis. Each array is indexed by',
      cd('varid'),
      'and has the same shape as',
      cd('var')
    ),
    varp = paste(
      'A',
      cl('SOMACollection'),
      'of',
      paste0(cl('SOMASparseNDArray'), 's'),
      'containing pairwise annotations on the variable axis and indexed with',
      cd('[varid_1, varid_2]')
    )
  ))
}
