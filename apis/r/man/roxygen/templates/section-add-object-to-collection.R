#' @section Adding new objects to a collection:
#'
#' The [`<%= class %>`] class provides a number of type-specific methods for
#' adding new a object to the collection, such as `add_new_sparse_ndarray()` and
#' `add_new_dataframe()`. These methods will create the new object and add it as
#' member of the `<%= class %>`. The new object will always inherit the parent
#' context (see [`SOMATileDBContext`]) and, by default, its platform
#' configuration (see [`PlatformConfig`]). However, the user can override the
#' default platform configuration by passing a custom configuration to the
#' `platform_config` argument.
#'
#' @section Carrara (TileDB v3) behavior:
#'
#' When working with Carrara URIs (`tiledb://workspace/teamspace/...`), child
#' objects created at a URI nested under a parent collection are **automatically
#' added** as members of the parent. This means:
#'
#' \itemize{
#'   \item You do not need to call `add_new_collection()` after creating a child
#'     at a nested URI—the child is already a member.
#'   \item For backward compatibility, calling `add_new_collection()` on an
#'     already-registered child is a **no-op** and will not cause an error.
#'   \item The member name must match the relative URI segment (e.g., creating
#'     at `parent_uri/child` automatically adds the child with key `"child"`).
#' }
