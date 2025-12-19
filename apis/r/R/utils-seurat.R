.anndata_to_seurat_reduc <- function(x, type = c("embeddings", "loadings")) {
  if (is.null(x)) {
    return(NULL)
  }
  stopifnot(is.character(x), is.character(type))
  type <- type[1L]
  type <- match.arg(type)
  return(switch(
    EXPR = type,
    embeddings = tolower(gsub(pattern = "^X_", replacement = "", x = x)),
    loadings = {
      m <- regexpr(pattern = "[[:upper:]]+", text = x)
      x <- tolower(unlist(regmatches(x = x, m = m)))
      x[x == "pc"] <- "pca"
      x[x == "ic"] <- "ica"
      x
    }
  ))
}

.load_seurat_command <- function(uns, ms_names) {
  key <- "seurat_commands"
  check_package("jsonlite")
  check_package("SeuratObject", version = .MINIMUM_SEURAT_VERSION())
  stopifnot(
    "'uns' must be a SOMACollection" = inherits(uns, what = "SOMACollection"),
    "'ms_names' must be a character vector with no empty strings" = is.character(
      ms_names
    ) &&
      all(nzchar(ms_names))
  )
  if (
    !(key %in%
      uns$names() &&
      inherits(logs <- uns$get(key), what = "SOMACollection"))
  ) {
    stop(errorCondition(
      "Cannot find a SOMACollection with command logs in 'uns'",
      class = c("noCommandLogsError", "missingCollectionError")
    ))
  }
  slots <- methods::slotNames(methods::getClassDef(
    "SeuratCommand",
    package = "SeuratObject"
  ))
  hint <- uns_hint("1d")
  lognames <- logs$names()
  commands <- stats::setNames(
    vector("list", length = length(lognames)),
    lognames
  )
  for (x in lognames) {
    soma_info(sprintf("Attempting to read command log %s", x))
    xdf <- logs$get(x)
    if (!inherits(xdf, "SOMADataFrame")) {
      soma_warn(sprintf("Log %s is invalid: not a SOMADataFrame", x))
      next
    }
    xhint <- tryCatch(xdf$get_metadata(names(hint)), error = function(...) "")
    if (xhint != hint[[1L]]) {
      soma_warn(sprintf(
        "Log %s is invalid: not a one-dimensional character data frame"
      ))
      next
    }
    soma_info("Reading in and decoding command log")
    tbl <- xdf$read(column_names = "values")$concat()
    enc <- as.data.frame(tbl)[["values"]]
    cmdlist <- jsonlite::fromJSON(enc)
    if (!(is.null(cmdlist$assay.used) || cmdlist$assay.used %in% ms_names)) {
      soma_info(sprintf("Skipping command log %s: assay used not requested", x))
      next
    }
    soma_info("Decoding command log parameters")
    for (param in names(cmdlist)) {
      cmdlist[[param]] <- if (param == "time.stamp") {
        ts <- sapply(
          jsonlite::fromJSON(cmdlist[[param]]),
          FUN = function(dt) {
            tryCatch(.decode_from_char(dt), error = function(...) dt)
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
        class(ts) <- c("POSIXlt", "POSIXt")
        as.POSIXct(ts)
      } else if (is.character(cmdlist[[param]])) {
        .decode_from_char(cmdlist[[param]])
      } else {
        cmdlist[[param]]
      }
    }
    soma_info("Assembling command log")
    params <- cmdlist[setdiff(names(cmdlist), slots)]
    cmdlist <- c(
      cmdlist[setdiff(names(cmdlist), names(params))],
      list(params = params)
    )
    commands[[x]] <- do.call(methods::new, c(cmdlist, Class = "SeuratCommand"))
  }
  commands <- Filter(Negate(is.null), x = commands)
  soma_info(sprintf("Returning %s command log(s)", length(commands)))
  idx <- order(sapply(commands, methods::slot, name = "time.stamp"))
  return(commands[idx])
}

.assay_version_hint <- function(type = c('v3', 'v5')) {
  type <- match.arg(type)
  return(list(soma_ecosystem_seurat_assay_version = type))
}

.assay_obs_hint <- function(assay) {
  stopifnot(
    "'assay' must be a single, non-empty character value" = is.character(
      assay
    ) &&
      length(assay) == 1L &&
      nzchar(assay) &&
      !is.na(assay)
  )
  return(sprintf("soma_ecosystem_seurat_assay_cells_%s", assay))
}

.layer_hint <- function(lyr) {
  stopifnot(
    "'lyr' must be a non-empty character vector" = is.character(lyr) &&
      length(lyr) &&
      all(nzchar(lyr)) &&
      !any(is.na(lyr))
  )
  if (length(lyr) > 1L) {
    lyr <- paste0('[', paste(dQuote(lyr, FALSE), collapse = ','), ']')
  }
  return(list(soma_ecosystem_seurat_v5_default_layers = lyr))
}

.ragged_array_hint <- function() {
  list(soma_ecosystem_seurat_v5_ragged = 'ragged')
}

.MINIMUM_SEURAT_VERSION <- function(repr = c('v', 'c')) {
  repr <- repr[1L]
  repr <- match.arg(arg = repr)
  version <- "4.1.0"
  return(switch(EXPR = repr, v = package_version(version), version))
}

.write_seurat_assay <- function(
  x,
  uri = NULL,
  soma_parent,
  ...,
  ingest_mode = 'write',
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  soma_context = NULL,
  relative = TRUE
) {
  check_package('SeuratObject', version = .MINIMUM_SEURAT_VERSION())
  if (!inherits(x, what = 'Assay')) {
    check_package('SeuratObject', version = '5.0.1')
    if (!inherits(x, 'Assay5')) {
      stop("'x' must be a Seurat v3 or v5 assay object")
    }
  }
  stopifnot(
    "'uri' must be a single character value" = is.null(uri) ||
      is_scalar_character(uri),
    "'soma_parent' must be a SOMACollection" = inherits(
      x = soma_parent,
      what = 'SOMACollectionBase'
    ),
    "'relative' must be a single logical value" = is_scalar_logical(relative)
  )

  v5 <- inherits(x, what = 'Assay5')

  # Find `shape` if and only if we're called from `write_soma.Seurat()`
  parents <- unique(sys.parents())
  idx <- which(vapply_lgl(
    X = parents,
    FUN = function(i) {
      return(identical(sys.function(i), write_soma.Seurat))
    }
  ))
  shape <- if (length(idx) == 1L) {
    get("shape", envir = sys.frame(parents[idx]))
  } else {
    NULL
  }
  shape <- rev(shape)

  # Create a proper URI
  uri <- uri %||%
    sub(pattern = '_$', replacement = '', x = SeuratObject::Key(x))
  uri <- .check_soma_uri(
    uri = uri,
    soma_parent = soma_parent,
    relative = relative
  )

  # Create the measurement
  ms <- SOMAMeasurementCreate(
    uri = uri,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    soma_context = soma_context
  )
  ms$set_metadata(.assay_version_hint(ifelse(v5, yes = 'v5', no = 'v3')))
  X <- if (!'X' %in% ms$names()) {
    SOMACollectionCreate(
      uri = file_path(ms$uri, 'X'),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      soma_context = soma_context
    )
  } else if (isTRUE(relative)) {
    SOMACollectionOpen(uri = file_path(ms$uri, 'X'), mode = 'WRITE')
  } else {
    ms$X
  }
  withCallingHandlers(
    .register_soma_object(X, soma_parent = ms, key = 'X', relative = relative),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(X$close(), add = TRUE, after = FALSE)

  # Write `X` matrices
  if (v5) {
    # Pull presence matrices from the v5 assay
    cells_matrix <- methods::slot(x, name = 'cells')
    features_matrix <- methods::slot(x, name = 'features')

    # Write `X` matrices
    for (layer in SeuratObject::Layers(x)) {
      ldat <- SeuratObject::LayerData(x, layer = layer)
      if (!.s3_method_defined("write_soma", class(ldat))) {
        rlang::warn(
          message = sprintf(
            "Unknown matrix type %s (layer %s)",
            class(ldat)[1L],
            layer
          ),
          class = "unknownMatrixTypeWarning"
        )
        next
      }
      type <- .type_hint(if (is.matrix(ldat)) {
        "matrix"
      } else if (methods::is(ldat, "IterableMatrix")) {
        "dgCMatrix"
      } else {
        class(ldat)
      })
      if (all(features_matrix[, layer]) && all(cells_matrix[, layer])) {
        soma_info(sprintf("Adding '%s' matrix as '%s'", layer, layer))
        tryCatch(
          expr = {
            arr <- write_soma(
              x = ldat,
              uri = layer,
              soma_parent = X,
              sparse = TRUE,
              transpose = TRUE,
              ingest_mode = ingest_mode,
              shape = shape,
              key = layer,
              platform_config = platform_config,
              tiledbsoma_ctx = tiledbsoma_ctx,
              soma_context = soma_context
            )
            arr$set_metadata(type)
          },
          error = function(err) {
            if (layer %in% SeuratObject::DefaultLayer(x)) {
              stop(err)
            }
            err_to_warn(err)
          }
        )
        next
      }
      ldat <- Matrix::t(as(ldat, "TsparseMatrix"))
      idx <- which(cells_matrix[, layer])
      jdx <- which(features_matrix[, layer])
      coo <- data.frame(
        soma_dim_0 = bit64::as.integer64(idx[ldat@i + 1L] - 1L),
        soma_dim_1 = bit64::as.integer64(jdx[ldat@j + 1L] - 1L),
        soma_data = ldat@x
      )
      atype <- arrow::infer_type(coo$soma_data)
      rt <- r_type_from_arrow_type(atype)
      if (rt == 'integer' && .is_integerish(coo$soma_data)) {
        coo$soma_data <- as.integer(coo$soma_data)
      }
      shape <- c(max(coo$soma_dim_0), max(coo$soma_dim_1)) + 1L
      arr <- X$add_new_sparse_ndarray(
        key = layer,
        type = atype,
        shape = as.integer(shape)
      )
      arr$.write_coordinates(coo)
      arr$set_metadata(.ragged_array_hint())
      arr$set_metadata(type)
    }
  } else {
    for (slot in c("counts", "data", "scale.data")) {
      mat <- SeuratObject::GetAssayData(x, slot)
      if (SeuratObject::IsMatrixEmpty(mat)) {
        next
      }

      # Skip 'data' slot if it's identical to 'counts'
      if (
        slot == "data" &&
          identical(mat, SeuratObject::GetAssayData(x, "counts"))
      ) {
        soma_info("Skipping 'data' slot because it's identical to 'counts'")
        next
      }

      # Pad 'scale.data'
      if (!identical(x = dim(mat), y = dim(x))) {
        soma_info(sprintf(
          "Padding layer '%s' to match dimensions of assay",
          slot
        ))
        mat <- pad_matrix(
          x = mat,
          rowidx = match(x = rownames(mat), table = rownames(x)),
          colidx = match(x = colnames(mat), table = colnames(x)),
          shape = dim(x),
          sparse = TRUE,
          rownames = rownames(x),
          colnames = colnames(x)
        )
      }

      layer <- gsub(pattern = '\\.', replacement = '_', x = slot)
      soma_info(sprintf("Adding '%s' matrix as '%s'", slot, layer))
      tryCatch(
        expr = write_soma(
          x = mat,
          uri = layer,
          soma_parent = X,
          sparse = TRUE,
          transpose = TRUE,
          ingest_mode = ingest_mode,
          shape = shape,
          key = layer,
          platform_config = platform_config,
          tiledbsoma_ctx = tiledbsoma_ctx,
          soma_context = soma_context
        ),
        error = function(err) {
          if (slot == 'data') {
            stop(err)
          }
          err_to_warn(err)
        }
      )
    }
  }

  # Write feature-level metadata
  var_df <- .df_index(
    x = x[[]],
    alt = 'features',
    axis = 'var',
    prefix = 'seurat'
  )
  var_df[[attr(x = var_df, which = 'index')]] <- rownames(x)
  soma_info("Adding feature-level metadata")
  write_soma(
    x = var_df,
    uri = 'var',
    soma_parent = ms,
    key = 'var',
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    soma_context = soma_context
  )

  # Check for any potentially-missed data
  if (!class(x)[1L] %in% c('Assay', 'Assay5')) {
    warning(
      paste(
        strwrap(paste0(
          "Extended assays (eg. ",
          class(x)[1L],
          ") are not fully supported; core Assay data has been written but ",
          "additional slots have been skipped"
        )),
        collapse = '\n'
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  return(ms)
}
