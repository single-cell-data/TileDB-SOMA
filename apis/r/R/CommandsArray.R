#' Storage for SeuratCommand
#'
#' TileDB array with options for storing SeuratCommand objects.
#' @export

CommandsArray <- R6::R6Class(
  classname = "CommandsArray",
  inherit = AnnotationArray,

  public = list(
    #' @field verbose Print status messages
    verbose = TRUE,

    #' @description Store the Seurat Command history to TileDB
    #' @param x a named list of Seurat Command objects
    from_named_list_of_commands = function(x) {
      stopifnot(
        "CommandsArray input must be named list of Seurat Command" = is.list(x) && is_named_list(x)
      )
      for (command in x) {
        stopifnot(
          "CommandsArray input must be named list of Seurat Command" = inherits(command, "SeuratCommand")
        )
      }

      # Convert from list of objects to list of dataframes.
      command_dataframes <- lapply(x, as.data.frame.SeuratCommand)

      # Convert from list of dataframes to single dataframe.
      command_dataframe <- do.call("rbind", command_dataframes)

      # Add an index column to preserve original ordering. Else, if we index on
      # name, commands will be sorted by name when read back even if sorted
      # differently when written.
      command_dataframe$index <- seq_len(nrow(command_dataframe))

      if (!self$exists()) {
        private$create_empty_array(command_dataframe, "index")
      } else {
        if (self$verbose) {
          message(
            sprintf("Updating existing %s at '%s'", self$class(), self$uri)
          )
        }
      }
      private$ingest_data(command_dataframe)
    },

    #' @description Retrieve the Seurat Command history from TileDB
    #' @return A named list of Seurat Command objects
    to_named_list_of_commands = function() {
      if (self$verbose) message("Reading command history into memory")

      arr <- self$tiledb_array(return_as = "data.frame")[]
      df <- arr[]
      rows <- split(df, seq(nrow(df)))

      named_list_of_commands <- lapply(rows, function(row) {
        new("SeuratCommand",
          name        = row$name,
          time.stamp  = row$time.stamp,
          assay.used  = row$assay.used,
          call.string = row$call.string,
          params      = jsonlite::fromJSON(row$encoded_params)
        )
      })
      command_names <- lapply(rows, function(row) { row$name })
      names(named_list_of_commands) <- command_names

      named_list_of_commands
    }
  )
)

# Coerce a Seurat Command to a data.frame, using JSON serialization of the
# command's parameters
#' @importFrom methods slot
as.data.frame.SeuratCommand <- function(x, row.names = FALSE, optional = FALSE, ...) {
    data.frame(
        name = slot(x, "name"),
        time.stamp = slot(x, "time.stamp"),
        assay.used = slot(x, "assay.used") %||% NA_character_,
        call.string = paste0(slot(x, "call.string"), collapse = ""),
        encoded_params = as.character(jsonlite::toJSON(slot(x, "params"), auto_unbox = TRUE))
    )
}
