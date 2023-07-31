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
