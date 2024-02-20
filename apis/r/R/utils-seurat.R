
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

#' @importFrom methods getClassDef slotNames
#'
.load_seurat_command <- function(uns, ms_names) {
  key <- 'seurat_commands'
  check_package('jsonlite')
  check_package('SeuratObject', version = .MINIMUM_SEURAT_VERSION())
  stopifnot(
    "'uns' must be a SOMACollection" = inherits(uns, what = 'SOMACollection'),
    "'ms_names' must be a character vector with no empty strings" = is.character(ms_names) &&
      all(nzchar(ms_names))
  )
  if (!(key %in% uns$names() && inherits(logs <- uns$get(key), what = 'SOMACollection'))) {
    stop(errorCondition(
      "Cannot find a SOMACollection with command logs in 'uns'",
      class = c("noCommandLogsError", "missingCollectionError")
    ))
  }
  slots <- slotNames(getClassDef('SeuratCommand', package = 'SeuratObject'))
  hint <- uns_hint('1d')
  lognames <- logs$names()
  commands <- setNames(vector('list', length = length(lognames)), lognames)
  for (x in lognames) {
    spdl::info("Attempting to read command log {}", x)
    xdf <- logs$get(x)
    if (!inherits(xdf, 'SOMADataFrame')) {
      spdl::warn("Log {} is invalid: not a SOMADataFrame", x)
      next
    }
    xhint <- tryCatch(xdf$get_metadata(names(hint)), error = function(...) '')
    if (xhint != hint[[1L]]) {
      spdl::warn("Log {} is invalid: not a one-dimensional character data frame")
      next
    }
    spdl::info("Reading in and decoding command log")
    tbl <- xdf$read(column_names = 'values')$concat()
    enc <- as.data.frame(tbl)[['values']]
    cmdlist <- jsonlite::fromJSON(enc)
    if (!(is.null(cmdlist$assay.used) || cmdlist$assay.used %in% ms_names)) {
      spdl::info("Skipping command log {}: assay used not requested", x)
      next
    }
    spdl::info("Decoding command log parameters")
    for (param in names(cmdlist)) {
      cmdlist[[param]] <- if (param == 'time.stamp') {
        ts <- sapply(
          jsonlite::fromJSON(cmdlist[[param]]),
          FUN = function(dt) tryCatch(.decode_from_char(dt), error = function(...) dt),
          simplify = FALSE,
          USE.NAMES = TRUE
        )
        class(ts) <- c('POSIXlt', 'POSIXt')
        as.POSIXct(ts)
      } else if (is.character(cmdlist[[param]])) {
        .decode_from_char(cmdlist[[param]])
      } else {
        cmdlist[[param]]
      }
    }
    spdl::info("Assembling command log")
    params <- cmdlist[setdiff(names(cmdlist), slots)]
    cmdlist <- c(cmdlist[setdiff(names(cmdlist), names(params))], list(params = params))
    commands[[x]] <- do.call(new, c(cmdlist, Class = 'SeuratCommand'))
  }
  commands <- Filter(Negate(is.null), x = commands)
  spdl::info("Returning {} command log(s)", length(commands))
  idx <- order(sapply(commands, methods::slot, name = 'time.stamp'))
  return(commands[idx])
}

.MINIMUM_SEURAT_VERSION <- function(repr = c('v', 'c')) {
  repr <- repr[1L]
  repr <- match.arg(arg = repr)
  version <- '4.1.0'
  return(switch(EXPR = repr, v = package_version(version), version))
}
