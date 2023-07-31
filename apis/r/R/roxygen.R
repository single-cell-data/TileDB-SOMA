#' @title Functions for Reusable Documentation
#'
#' @description These functions are designed to be embedded in \pkg{roxygen2}
#' blocks using the \\Sexp macro; \pkg{roxygen2} templates have some issues
#' when used to document \code{R6} classes, so these function provide a more
#' stable shared documentation interface. Each function returns an
#' Rd (R documentation) formatted string. Some functions take parameters to
#' allow customization of generally-sharable documentation
#'
#' @noRd
#'
NULL

#' Document Virtual Classes
#'
#' @return An Rd-formatted string stating the the class being documented is
#' virtual and cannot be instantiated
#'
#' @noRd
#'
rd_return_virtual <- function() {
  return("This is a \\strong{virtual} class and cannot be directly instantiated")
}

#' Document Atomic Types
#'
#' @return An Rd-formatted itemized list of atomic types in \R
#'
#' @noRd
#'
rd_atomic <- function() {
  return(paste(
    '\\itemize{',
    paste0('\\item \\dQuote{\\code{', .SCALAR_TYPES(), '}}', collapse = '\n'),
    '}',
    sep = '\n'
  ))
}

#' Document Ephemeral Class Descriptions
#'
#' @param cls Type of ephemeral class to document; choose from:
#' \itemize{
#'  \item \dQuote{\code{collection}}
#'  \item \dQuote{\code{experiment}}
#'  \item \dQuote{\code{measurement}}
#' }
#' @param base Is the class being documented a base class or a final class
#'
#' @return An Rd-formatted description of the ephemeral class \code{cls};
#' includes a link to the documentation for the on-disk version of \code{cls}
#'
#' @noRd
#'
rd_ephemeral_cls <- function(
  cls = c('collection', 'experiment', 'measurement'),
  base = FALSE
) {
  stopifnot(is.character(cls), is_scalar_logical(base))
  cls <- match.arg(arg = cls)
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

#' Document Ephemeral Method Descriptions
#'
#' @return An Rd-formatted string stating that the method associated
#' with an ephemeral class is a stub that exists solely for compatibility
#' with on-disk versions
#'
#' @noRd
#'
rd_ephemeral_desc <- function() {
  return("Dummy method for ephemeral cobjects for compatibility with SOMA collections")
}

#' Document Ephemeral Method Errors
#'
#' @return An Rd-formatted string stating that the method associated with an
#' ephemeral class throws an error rather than returning a value
#'
#' @noRd
#'
rd_ephemeral_error <- function() {
  return("Throws an error as this method is not supported by ephemeral objects")
}

#' Document Ephemeral Fields
#'
#' @return An Rd-formatted string stating that the field associated
#' with an ephemeral class is a stub that exists solely for compatibility
#' with on-disk versions
#'
#' @noRd
#'
rd_ephemeral_field <- function() {
  return("Dummy field for ephemeral objects for compatibility with SOMA collections")
}

#' Document Ephemeral Method Parameters
#'
#' @return An Rd-formatted string stating that the parameters for methods
#' associated with an ephemeral class are ignored and exist solely for
#' compatibility with on-disk versions
#'
#' @noRd
#'
rd_ephemeral_param <- function() {
  return("Ignored for ephemeral objects")
}

#' Document SOMA Fields
#'
#' @param field Name of field to document; choose from:
#' \itemize{
#'  \item \dQuote{\code{X}}
#'  \item \dQuote{\code{ms}}
#'  \item \dQuote{\code{obs}}
#'  \item \dQuote{\code{obsm}}
#'  \item \dQuote{\code{obsp}}
#'  \item \dQuote{\code{var}}
#'  \item \dQuote{\code{varm}}
#'  \item \dQuote{\code{varp}}
#' }
#'
#' @return An Rd-formatted string describing the
#' structure of the SOMA field \code{field}
#'
#' @noRd
#'
rd_soma_field <- function(
  field = c('X', 'ms', 'obs', 'obsm', 'obsp', 'var', 'varm', 'varp')
) {
  stopifnot(is.character(field))
  field <- match.arg(arg = field)
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

#' Document \code{obs_index} and \code{var_index} for SOMA Outgestion
#'
#' @param type Type of outgestion being performed; choose from:
#' \itemize{
#'  \item \dQuote{\code{v3}} for Seurat v3 outgestion
#'  \item \dQuote{\code{sce}} for SingleCellExperiment outgestion
#' }
#' @param axis Axis data frame for index being documented; choose from:
#' \itemize{
#'  \item \dQuote{\code{obs}}
#'  \item \dQuote{\code{var}}
#' }
#'
#' @return An Rd-formatted string stating that \code{obs_index} or
#' \code{var_index} will be pulled from \code{obs} or \code{var} to use as
#' cell or feature names, and that by default the names will be pulled from
#' \dQuote{\code{soma_joinids}}
#'
#' @noRd
#'
rd_outgest_index <- function(type = c('v3', 'sce'), axis = c('obs', 'var')) {
  type <- match.arg(arg = type)
  axis <- match.arg(arg = axis)
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

#' Document Metadata Name Parameters for SOMA Outgestion
#'
#' @param type Type of outgestion being performed; choose from:
#' \itemize{
#'  \item \dQuote{\code{v3}} for Seurat v3 outgestion
#'  \item \dQuote{\code{sce}} for SingleCellExperiment outgestion
#' }
#' @param axis Axis data frame to be documented, choose from:
#' \itemize{
#'  \item \dQuote{\code{obs}}
#'  \item \dQuote{\code{var}}
#' }
#'
#' @return An Rd-formatted string stating that \code{obs_column_names} or
#' \code{var_column_names} will be pulled from \code{obs} or \code{var} and
#' where in the resulting objects the data frames will be stored. Also states
#' that by default, all columns from \code{obs} or \code{var} will be loaded
#'
#' @noRd
#'
rd_outgest_metadata_names <- function(type = c('v3', 'sce'), axis = c('obs', 'var')) {
  type <- match.arg(arg = type)
  axis <- match.arg(arg = axis)
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

#' Document *m Layers for SOMA Outgestion
#'
#' @param type Type of outgestion being performed; choose from:
#' \itemize{
#'  \item \dQuote{\code{v3}} for Seurat v3 outgestion
#'  \item \dQuote{\code{sce}} for SingleCellExperiment outgestion
#' }
#' @param axis Axis *m layers to be documented, choose from:
#' \itemize{
#'  \item \dQuote{\code{obsm}}
#'  \item \dQuote{\code{varm}}
#' }
#'
#' @return An Rd-formatted string for documenting \code{obsm_layers} or
#' \code{varm_layers}; includes examples of how to set the arguments and
#' details the default behavior
#'
#' @noRd
#'
rd_outgest_mlayers <- function(type = c('v3', 'sce'), axis = c('obsm', 'varm')) {
  type <- match.arg(arg = type)
  axis <- match.arg(arg = axis)
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

#' Document *p Layers for SOMA Outgestion
#'
#' @param type Type of outgestion being performed; choose from:
#' \itemize{
#'  \item \dQuote{\code{v3}} for Seurat v3 outgestion
#'  \item \dQuote{\code{sce}} for SingleCellExperiment outgestion
#' }
#' @param axis Axis *p layers to be documented, choose from:
#' \itemize{
#'  \item \dQuote{\code{obsp}}
#'  \item \dQuote{\code{varp}}
#' }
#'
#' @return An Rd-formatted string for documenting \code{obsp_layers} or
#' \code{varp_layers} and details the default behavior
#'
#' @noRd
#'
rd_outgest_players <- function(type = c('v3', 'sce'), axis = c('obsp', 'varp')) {
  type <- match.arg(arg = type)
  axis <- match.arg(arg = axis)
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

#' Document X Layers for SOMA Outgestion
#'
#' @param type Type of outgestion being performed; choose from:
#' \itemize{
#'  \item \dQuote{\code{v3}} for Seurat v3 outgestion
#'  \item \dQuote{\code{sce}} for SingleCellExperiment outgestion
#' }
#'
#' @return An Rd-formatted string for documenting \code{X_layers}; includes
#' examples of how \code{X_layers} should be formatted
#'
#' @noRd
#'
rd_outgest_xlayers <- function(type = c('v3', 'sce')) {
  type <- match.arg(arg = type)
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
