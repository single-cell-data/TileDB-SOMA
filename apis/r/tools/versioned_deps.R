#!/usr/bin/env Rscript

if (!requireNamespace("rprojroot", quietly = TRUE)) {
  stop("Cannot find rprojroot", call. = FALSE)
}

root <- rprojroot::find_package_root_file()

deps <- read.dcf(file.path(root, "DESCRIPTION"))[1L, c("Depends", "Imports")]
deps <- Filter(
  f = nzchar,
  x = trimws(x = unlist(
    x = strsplit(
      x = gsub(pattern = ",", replacement = "\n", x = deps),
      split = "\n"
    ),
    use.names = FALSE
  ))
)

pattern <- sprintf(
  fmt = "\\(>= %s\\)$",
  .standard_regexps()$valid_package_version
)
versions <- vector(mode = "character", length = length(x = deps))
for (i in seq_along(along.with = deps)) {
  if (!grepl(pattern = pattern, x = deps[i])) {
    versions[i] <- NA_character_
    names(x = versions)[i] <- deps[i]
    next
  }
  version <- regmatches(
    x = deps[i],
    m = regexpr(pattern = pattern, text = deps[i])
  )
  versions[i] <- trimws(x = sub(
    pattern = "\\(>=",
    replacement = "",
    x = sub(pattern = "\\)$", replacement = "", x = version)
  ))
  names(x = versions)[i] <- trimws(x = sub(
    pattern = pattern,
    replacement = "",
    x = deps[i]
  ))
}
versions <- Filter(f = Negate(f = is.na), x = versions)

if (length(idx <- which(x = names(x = versions) == "R"))) {
  versions <- versions[-idx]
}
cat(names(versions), sep = "\n")
