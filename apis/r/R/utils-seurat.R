#' Create \pkg{Seurat}-Style Names
#'
#' Convert AnnData-style names to \pkg{Seurat}-style names
#'
#' @param x A character vector of names
#' @param type Type of conversion to perform; choose from:
#' \itemize{
#'  \item \dQuote{\code{embeddings}}: convert AnnData-style \code{obsm} names
#'  \item \dQuote{\code{loadings}}: convert AnnData-style \code{varm} names
#' }
#'
#' @return \code{x} with names converted to Seurat-style names
#' based on \code{type}
#'
#' @keywords internal
#'
#' @noRd
#'
.anndata_to_seurat_reduc <- function(x, type = c('embeddings', 'loadings')) {
  if (is.null(x)) {
    return(NULL)
  }
  stopifnot(is.character(x), is.character(type))
  type <- type[1L]
  type <- match.arg(type)
  return(switch(
    EXPR = type,
    embeddings = tolower(gsub(pattern = '^X_', replacement = '', x = x)),
    loadings = {
      m <- regexpr(pattern = '[[:upper:]]+', text = x)
      x <- tolower(unlist(regmatches(x = x, m = m)))
      x[x == 'pc'] <- 'pca'
      x[x == 'ic'] <- 'ica'
      x
    }
  ))
}

#' Check \pkg{SeuratObject} Installation Status
#'
#' Check to see that valid version of \pkg{SeuratObject} is installed
#'
#' @inheritParams base::requireNamespace
#'
#' @return If \code{quietly}, then invisibly returns the installation status;
#' otherwise, errors if a valid version \pkg{SeuratObject} is unavailable or
#' invisibly returns \code{TRUE}
#'
#' @keywords internal
#'
#' @noRd
#'
.check_seurat_installed <- function(quietly = FALSE) {
  pkg <- 'SeuratObject'
  checks <- c(
    installed = requireNamespace(pkg, quietly = TRUE),
    version = tryCatch(
      expr = utils::packageVersion(pkg) >= .MINIMUM_SEURAT_VERSION(),
      error = function(...) {
        return(FALSE)
      }
    )
  )
  if (isTRUE(quietly)) {
    return(invisible(all(checks)))
  }
  if (!checks['installed']) {
    stop(sQuote(pkg), " must be installed", call. = FALSE)
  }
  if (!checks['version']) {
    stop(
      sQuote(pkg),
      " must be version ",
      .MINIMUM_SEURAT_VERSION('c'),
      " or higher",
      call. = FALSE
    )
  }
  return(invisible(TRUE))
}

#' Load \link[SeuratObject:Graph]{Nearest-Neighbor Graphs}
#'
#' Read in \link[SeuratObject:Graph]{nearest-neighbor graph} objects from a
#' \code{\link{SOMAExperimentAxisQuery}} object
#'
#' @param query A \code{\link{SOMAExperimentAxisQuery}} object
#' @param obsp_layers Names of arrays in \code{obsp} to load in as
#' \code{\link[SeuratObject]{Graph}s}; by default, loads all graphs
#' @template param-obs-index
#'
#' @return A named list of \code{\link[SeuratObject]{Graph}} objects
#'
#' @keywords internal
#'
#' @noRd
#'
.get_seurat_graphs <- function(query, obsp_layers = NULL, obs_index = NULL) {
  .check_seurat_installed()
  stopifnot(
    "'query' must be a SOMAExperimentAxisQuery object" = inherits(
      x = query,
      what = 'SOMAExperimentAxisQuery'
    ),
    "'obsp_layers' must be a character vector" = is_character_or_null(obsp_layers),
    "'obs_index' must be a single character value" = is.null(obs_index) ||
      (is_scalar_character(obs_index) && !is.na(obs_index))
  )
  ms_graphs <- tryCatch(expr = query$ms$obsp$names(), error = null)
  if (is.null(ms_graphs)) {
    stop("No nearest-neighbor graphs found", call. = FALSE)
  }
  obsp_layers <- obsp_layers %||% ms_graphs
  res <- stats::setNames(
    object = vector(mode = 'list', length = length(obsp_layers)),
    nm = obsp_layers
  )
  for (grph in obsp_layers) {
    mat <- tryCatch(
      expr = query$to_seurat_graph(obsp_layer = grph, obs_index = obs_index),
      error = function(e) {
        warning(conditionMessage(e), call. = FALSE, immediate. = TRUE)
        return(NULL)
      }
    )
    if (is.null(mat)) {
      next
    }
    res[[grph]] <- mat
  }
  return(Filter(f = Negate(is.null), x = res))
}

#' Load \link[SeuratObject:DimReduc]{Dimensional Reductions}
#'
#' Read in \link[SeuratObject:DimReduc]{dimensional reduction} objects from a
#' \code{\link{SOMAExperimentAxisQuery}} object
#'
#' @param query A \code{\link{SOMAExperimentAxisQuery}} object
#' @param obsm_layers Names of arrays in \code{obsm} to load in as the
#' cell embeddings; pass \code{FALSE} to suppress loading in any
#' dimensional reductions; by default, loads all dimensional
#' reduction information
#' @param varm_layers Named vector of arrays in \code{varm} to load in as
#' the feature loadings; names must be names of array in \code{obsm} (eg.
#' \code{varm_layers = c(X_pca = 'PCs')}); will try to determine
#' \code{varm_layers} from \code{obsm_layers}
#' @template param-obs-index
#' @template param-var-index
#'
#' @return A named list of \code{\link[SeuratObject]{DimReduc}} objects where
#' the names are the \pkg{Seurat}-style names
#'
#' @keywords internal
#'
#' @noRd
#'
.get_seurat_reductions <- function(
  query,
  obsm_layers = NULL,
  varm_layers = NULL,
  obs_index = NULL,
  var_index = NULL
) {
  .check_seurat_installed()
  stopifnot(
    "'query' must be a SOMAExperimentAxisQuery object" = inherits(
      x = query,
      what = 'SOMAExperimentAxisQuery'
    ),
    "'obsm_layers' must be a character vector" = is_character_or_null(obsm_layers),
    "'varm_layers' must be a named character vector" = is.null(varm_layers) ||
      (is.character(varm_layers) && is_named(varm_layers, allow_empty = FALSE)) ||
      is_scalar_logical(varm_layers),
    "'obs_index' must be a single character value" = is.null(obs_index) ||
      (is_scalar_character(obs_index) && !is.na(obs_index)),
    "'var_index' must be a single character value" = is.null(var_index) ||
      (is_scalar_character(var_index) && !is.na(var_index))
  )
  ms_embed <- tryCatch(expr = query$ms$obsm$names(), error = null)
  if (is.null(ms_embed)) {
    stop("No reductions found", call. = FALSE)
  }
  names(ms_embed) <- .anndata_to_seurat_reduc(ms_embed)
  obsm_layers <- obsm_layers %||% ms_embed
  res <- stats::setNames(
    object = vector(mode = 'list', length = length(obsm_layers)),
    nm = .anndata_to_seurat_reduc(obsm_layers)
  )
  # Match loadings to embeddings
  ms_load <- tryCatch(expr = query$ms$varm$names(), error = null)
  if (isTRUE(varm_layers)) {
    varm_layers <- NULL
  } else if (rlang::is_na(varm_layers)) {
    varm_layers <- FALSE
  }
  if (is.null(ms_load) && !isFALSE(varm_layers)) {
    warning("No loadings found", call. = FALSE, immediate. = TRUE)
    varm_layers <- FALSE
  }
  if (!isFALSE(varm_layers)) {
    names(ms_load) <- ms_embed[.anndata_to_seurat_reduc(ms_load, 'loadings')]
    varm_layers <- varm_layers %||% ms_load
    reduc_misisng <- setdiff(x = names(varm_layers), y = names(ms_load))
    if (length(reduc_misisng) == length(varm_layers)) {
      warning(
        "None of the reductions specified in 'varm_layers' can be found",
        call. = FALSE,
        immediate. = TRUE
      )
      varm_layers <- FALSE
    } else if (length(reduc_misisng)) {
      warning(
        paste(
          strwrap(paste(
            "The reductions for the following loadings cannot be found in 'varm':",
            sQuote(varm_layers[reduc_misisng]),
            collapse = ', '
          )),
          collapse = '\n'
        ),
        call. = FALSE,
        immediate. = TRUE
      )
      varm_layers <- varm_layers[!names(varm_layers) %in% reduc_misisng]
    }
  }
  # Read in reductions and add to `object`
  for (embed in obsm_layers) {
    if (embed %in% names(ms_embed)) {
      embed <- ms_embed[embed]
    }
    rname <- .anndata_to_seurat_reduc(embed)
    reduc <- tryCatch(
      expr = query$to_seurat_reduction(
        obsm_layer = embed,
        varm_layer = ifelse(
          embed %in% names(varm_layers),
          yes = varm_layers[embed],
          no = FALSE
        ),
        obs_index = obs_index,
        var_index = var_index
      ),
      error = function(e) {
        warning(conditionMessage(e), call. = FALSE, immediate. = TRUE)
        return(NULL)
      }
    )
    if (is.null(reduc)) {
      next
    }
    res[[rname]] <- reduc
  }
  return(Filter(f = Negate(is.null), x = res))
}

#' Minimum Version of \pkg{SeuratObject}
#'
#' Fetch the minimum required version of \pkg{SeuratObject}
#'
#' @param repr Representation of the version; choose from:
#' \itemize{
#'  \item \dQuote{\code{v}} to return a \code{\link[base]{package_version}}
#'  \item \dQuote{\code{c}} to return a \code{\link[base]{character}}
#' }
#'
#' @return The minimum required version of \pkg{SeuratObject}
#'
#' @keywords internal
#'
#' @noRd
#'
.MINIMUM_SEURAT_VERSION <- function(repr = c('v', 'c')) {
  repr <- repr[1L]
  repr <- match.arg(arg = repr)
  version <- '4.1.0'
  return(switch(EXPR = repr, v = package_version(version), version))
}
