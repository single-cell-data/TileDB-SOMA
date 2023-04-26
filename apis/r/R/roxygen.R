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

rd_outgest_index <- function(type = 'v3', axis = 'obs') {
  type <- type[1L]
  type <- match.arg(arg = type, choices = c('v3', 'sce'))
  axis <- axis[1L]
  axis <- match.arg(arg = axis, choices = c('obs', 'var'))
  label <- switch(EXPR = axis, obs = 'cell', var = 'feature')
  return(paste0(
    'Name of column in \\code{',
    axis,
    '} to add as ',
    label,
    ' names; uses \\code{paste0("',
    label,
    '", ',
    axis,
    '_joinids())} by default'
  ))
}

rd_outgest_mdnames <- function(type = 'v3', axis = 'obs') {
  type <- type[1L]
  type <- match.arg(arg = type, choices = c('v3', 'sce'))
  axis <- axis[1L]
  axis <- match.arg(arg = axis, choices = c('obs', 'var'))
  return(paste0(
    'Names of columns in \\code{',
    axis,
    '} to add as ',
    switch(
      EXPR = type,
      v3 = paste0(
        switch(EXPR = axis, obs = 'cell', var = 'feature'),
        '-level meta data'
      ),
      sce = switch(EXPR = axis, obs = '\\code{colData}', var = '\\code{rowData}')
    ),
    '; by default, loads all columns'
  ))
}

rd_outgest_mlayers <- function(type = 'v3', axis = 'obsm') {
  type <- type[1L]
  type <- match.arg(arg = type, choices = c('v3', 'sce'))
  axis <- axis[1L]
  axis <- match.arg(arg = axis, choices = c('obsm', 'varm'))
  dr <- switch(EXPR = type, v3 = 'dimensional reduction', sce = 'reduced dimension')
  label <- switch(
    EXPR = type,
    v3 = switch(EXPR = axis, obsm = 'cell embeddings', varm = 'feature loadings'),
    sce = paste0(dr, 's')
  )
  unnamed <- paste0(
    'Names of arrays in \\code{',
    axis,
    '} to add as the ',
    label
  )
  intro <- switch(
    EXPR = axis,
    obsm = unnamed,
    varm = switch(
      EXPR = type,
      v3 = paste0(
        'Named vector of arrays in \\code{',
        axis,
        '} to load in as the feature loadings; names must be names of arrays in',
        ' \\code{obsm} (eg. \\code{varm_layers = c(X_pca = "PCs")})'
      ),
      unnamed
    )
  )
  suppress <- paste(
    'pass \\code{FALSE} to suppress loading in any',
    switch(EXPR = axis, obsm = paste0(dr, 's'), varm = label)
  )
  default <- switch(
    EXPR = axis,
    obsm = paste0(
      'by default, loads all ',
      dr,
      switch(EXPR = type, v3 = ' information', sce = 's')
    ),
    varm = 'will try to determine \\code{varm_layers} from \\code{obsm_layers}'
  )
  return(paste(intro, suppress, default, sep = '; '))
}

rd_outgest_players <- function(type = 'v3', axis = 'obsp') {
  type <- type[1L]
  type <- match.arg(arg = type, choices = c('v3', 'sce'))
  axis <- axis[1L]
  axis <- match.arg(arg = axis, choices = c('obsp', 'varp'))
  return(paste0(
    'Names of arrays in \\code{',
    axis,
    '} to load in as \\code{\\link[',
    switch(EXPR = type, v3 = 'SeuratObject', sce = 'S4Vectors'),
    ']{',
    switch(EXPR = type, v3 = 'Graph}s', sce = 'SelfHits}'),
    '}; by default, loads all ',
    switch(EXPR = axis, obsp = 'graphs', varp = 'networks')
  ))
}

rd_outgest_xlayers <- function(type = 'v3') {
  type <- type[1L]
  type <- match.arg(arg = type, choices = c('v3', 'sce'))
  return(switch(
    EXPR = type,
    v3 = paste(
      "A named character of X layers to add to the Seurat assay",
      "where the names are the names of Seurat slots and the values are the names",
      "of layers within \\code{X}; names should be one of:",
      "\\itemize{",
      "\\item \\dQuote{\\code{counts}} to add the layer as \\code{counts}",
      "\\item \\dQuote{\\code{data}} to add the layer as \\code{data}",
      "\\item \\dQuote{\\code{scale.data}} to add the layer as \\code{scale.data}",
      "}",
      "At least one of \\dQuote{\\code{counts}} or \\dQuote{\\code{data}} is required"
    ),
    sce = paste(
      'A character vector of X layers to add as assays in the main experiment;',
      'may optionally be named to set the name of the resulting assay',
      '(eg. \\code{X_layers = c(counts = "raw")} will load in X layer',
      '\\dQuote{\\code{raw}} as assay \\dQuote{\\code{counts}});',
      'by default, loads in all X layers'
    )
  ))
}
