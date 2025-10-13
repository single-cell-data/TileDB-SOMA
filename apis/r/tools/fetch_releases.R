#!/usr/bin/env Rscript

deps <- c("gh", "rprojroot")
for (i in deps) {
  if (!requireNamespace(i, quietly = TRUE)) {
    stop("Cannot find package ", sQuote(x = i, q = c("{", "}", "", "")))
  }
}

if (!gh::gh_token_exists()) {
  stop("No GitHub token found")
}

root <- rprojroot::find_package_root_file()
dir.create(
  path = file.path(root, "inst", "extdata"),
  showWarnings = FALSE,
  recursive = TRUE
)
repo <- Sys.getenv("GITHUB_REPOSITORY", "single-cell-data/TileDB-SOMA")
has_next <- sprintf(fmt = "rel=%s", dQuote(x = "next", q = FALSE))

resp <- list(link = has_next)
releases <- list()
message("Fetching releases for '", repo, "'")
while (grepl(pattern = has_next, x = resp$link)) {
  endpoint <- if (length(x = resp) == 1L) {
    sprintf(fmt = "/repos/%s/releases?per_page=100", repo)
  } else {
    link <- Filter(
      f = function(x) {
        return(grepl(pattern = has_next, x = x))
      },
      x = trimws(x = unlist(x = strsplit(x = resp$link, split = ",")))
    )
    sub(
      pattern = sprintf(fmt = ">; %s$", has_next),
      replacement = "",
      x = sub(pattern = "^<", replacement = "", x = link)
    )
  }
  response <- gh::gh(endpoint = endpoint)
  releases <- c(
    releases,
    Filter(
      f = Negate(f = is.null),
      lapply(
        X = response,
        FUN = function(x) {
          if (x$draft || x$prerelease) {
            return(NULL)
          }
          if (is.na(x = numeric_version(x = x$tag_name, strict = FALSE))) {
            return(NULL)
          }
          return(data.frame(
            Version = x$tag_name,
            Date = format(
              strptime(x = x$created_at, format = "%Y-%m-%dT%H:%M:%SZ"),
              format = "%Y-%m-%d"
            )
          ))
        }
      )
    )
  )
  resp <- attr(x = response, which = "response", exact = TRUE)
}

message("Assembling releases data frame")
releases <- do.call(what = rbind, args = releases)
file <- file.path(root, "inst", "extdata", "releases.dcf")
message("Writing releases to '", file,"'")
write.dcf(x = releases, file = file)
