#!/usr/bin/env Rscript

wd <- Sys.getenv("_CI_WORKING_DIRECTORY_", unset = getwd())
pkgtype <- Sys.getenv("_CI_PKG_TYPE_", unset = getOption("pkgType"))

if (isTRUE(as.logical(Sys.getenv("_CI_USE_BIOC_")))) {
  message("Setting up Bioconductor repos")
  setRepositories(ind = 1:3)
}

desc <- as.list(x = as.data.frame(x = read.dcf(
  file = file.path(wd, "DESCRIPTION"),
  fields = c(
    "Depends",
    "Imports",
    "LinkingTo",
    "Suggests",
    "Enhances",
    "Additional_repositories"
  )
)))

desc <- Filter(f = Negate(f = is.na), x = desc)

if (!is.null(desc$Additional_repositories)) {
  desc$Additional_repositories <- trimws(x = unlist(x = strsplit(
    x = desc$Additional_repositories,
    split = ','
  )))
}

repos <- c(desc$Additional_repositories, getOption("repos"))
desc$Additional_repositories <- NULL

message("Identifying dependencies to install")
desc <- lapply(
  X = desc,
  FUN = function(x) {
    x <- trimws(x = unlist(x = strsplit(x = x, split = ",")))
    x <- strsplit(x = x, split = " \\(>= ")
    for (i in seq_along(along.with = x)) {
      x[[i]][2L] <- x[[i]][2L]
    }
    x <- as.data.frame(x = t(x = data.frame(x)))
    row.names(x = x) <- NULL
    names(x = x) <- c("Package", "Version")
    x$Version <- sub(pattern = "\\)$", replacement = "", x = x$Version)
    if (length(idx <- which(x = x$Package == "R"))) {
      x <- x[-idx, , drop = FALSE]
    }
    if (!nrow(x = x)) {
      return(x)
    }
    x$Version[is.na(x = x$Version)] <- "0.0.0"
    x$Installed <- logical(length = nrow(x = x))
    pb <- txtProgressBar(max = nrow(x), style = 3L, file = stderr())
    on.exit(expr = close(con = pb), add = TRUE)
    for (i in seq_len(length.out = nrow(x = x))) {
      x$Installed[i] <- requireNamespace(x$Package[i], quietly = TRUE) &&
        utils::packageVersion(x$Package[i]) >= x$Version[i]
      setTxtProgressBar(pb = pb, value = i)
    }
    return(x)
  }
)

desc <- do.call(what = "rbind", args = desc)
row.names(x = desc) <- NULL

desc <- desc[!desc$Installed, , drop = FALSE]

if (nrow(x = desc)) {
  message(
    "Found ",
    nrow(x = desc),
    " missing ",
    ngettext(n = nrow(x = desc), msg1 = "dependency", msg2 = "dependencies")
  )
  message("Finding additional downstream dependencies")
  db <- utils::available.packages(repos = repos, type = pkgtype)
  deps <- unique(x = unlist(
    x = tools::package_dependencies(packages = desc$Package, db = db),
    use.names = FALSE
  ))
  deps <- setdiff(x = deps, y = c(utils::sessionInfo()$basePkgs, desc$Package))
  deps <- Filter(
    f = Negate(f = \(x) isTRUE(requireNamespace(x, quietly = TRUE))),
    x = deps
  )
  if (length(x = deps)) {
    message(
      "Found ",
      length(x = deps),
      " additional ",
      ngettext(n = length(x = deps), msg1 = "pacakge", msg2 = "packages")
    )
  }
  pkgs <- union(x = desc$Package, y = deps)
  message(
    "Installing ",
    length(x = pkgs),
    ngettext(n = length(x = pkgs), msg1 = " package", msg2 = " packages")
  )
  install.packages(desc$Package, type = pkgtype, repos = repos)
}
