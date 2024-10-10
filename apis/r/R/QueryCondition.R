#  MIT License
#
#  Copyright (c) 2021-2024 TileDB Inc.
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.

# ================================================================
#' Create a 'tiledbsoma_query_condition' object from an expression
#'
#' The grammar for query conditions is constrained to the operators
#' (\code{">"}, \code{">="}, \code{"<"}, \code{"<="}, \code{"=="},
#' \code{"!="}, \code{"%in%"}, \code{"%nin%"}), and three boolean operators
#' (\code{"&&"}, also as \code{"&"}, (\code{"||"}, also as \code{"|"}, and
#' \code{"!"} for negation.  Note that we locally define \code{"%nin%"} as
#' \code{Negate()} call around \code{%in%)} which extends R a little for this
#' use case.
#'
#' Expressions, in the R language syntax, are parsed locally by this function.
#'
#' @param expr An expression that is understood by the TileDB grammar for
#' query conditions.
#'
#' @param schema The Arrow schema for the array for which a query
#' condition is being prepared. This is necessary to obtain type information
#' about left-hand sides of query expressions.
#'
#' @param strict A boolean toogle to, if set, errors if a non-existing
#' attribute is selected or filtered on, defaults to 'TRUE'; if 'FALSE' a
#' warning is shown but execution proceeds.
#'
#' @param somactx SOMAContext pointer.
#'
#' @return A `tiledbsoma_query_condition` object.
#' @export
parse_query_condition_new <- function(
  expr,
  schema,
  strict=TRUE,
  somactx
  ) {

  stopifnot("The schema argument must be an Arrow Schema" =
    is(schema, "ArrowObject") &&
    is(schema, "Schema"))

    # ----------------------------------------------------------------
    # Helpers for walking the parse tree

    # Operators
    `%!in%` <- Negate(`%in%`)
    .is_in_operator <- function(node) {
        return(tolower(as.character(node)) %in% c("%in%", "%nin%"))
    }
    .is_comparison_operator <- function(node) {
        return(tolower(as.character(node)) %in% c(">", ">=", "<", "<=", "==", "!=", "%in%", "%nin%"))
    }
    .is_boolean_operator <- function(node) {
        return(as.character(node) %in% c("&&", "||", "!", "&", "|"))
    }

    # Leaf nodes
    .is_ascii <- function(node) {
        return(grepl("^[[:alnum:]_]+$", node))
    }
    .is_integer <- function(node) {
        return(grepl("^[[:digit:]]+$", as.character(node)))
    }
    .is_double <- function(node) {
        return(grepl("^[[:digit:]\\.]+$", as.character(node)) && length(grepRaw(".", as.character(node), fixed = TRUE, all = TRUE)) == 1)
    }

    .error_function <- if (strict) stop else warning

    .map_op_to_character <- function(x) {
        return(switch(x, `>` = "GT", `>=` = "GE", `<` = "LT", `<=` = "LE", `==` = "EQ", `!=` = "NE"))
    }

    .map_bool_to_character <- function(x) {
        return(switch(x, `&&` = "AND", `&` = "AND", `||` = "OR", `|` = "OR", `!` = "NOT"))
    }

    # ----------------------------------------------------------------
    # Map the R parse tree (from base-r `substitute`) to a TileDB core QueryCondition

    .parse_tree_to_qc <- function(node, debug=FALSE) {
        if (is.symbol(node)) {
            stop("Unexpected symbol in expression: ", format(node))

        } else if (.is_boolean_operator(node[1])) {
            spdl::debug("[parseqc] boolop [{}] [{}] [{}]",
                as.character(node[2]),
                as.character(node[1]),
                as.character(node[3]))

            return(tiledbsoma_query_condition_combine(
                .parse_tree_to_qc(node[[2]]),
                .parse_tree_to_qc(node[[3]]),
                .map_bool_to_character(as.character(node[1])),
                somactx))

        } else if (.is_in_operator(node[1])) {
            spdl::debug("[parseqc] inop [{}] [{}] [{}]",
                as.character(node[2]),
                as.character(node[1]),
                as.character(node[3]))

            attr_name <- as.character(node[2])
            r_op_name <- tolower(as.character(node[1]))
            tdb_op_name <- if (r_op_name == "%in%") "IN" else "NOT_IN"

            # XXX EXTRACT HELPER
            arrow_field <- schema[[attr_name]]
            if (is.null(arrow_field)) {
                .error_function("No attribute '", attr_name, "' is present.", call. = FALSE)
            }
            arrow_type_name <- arrow_field$type$name
            is_enum <- is(arrow_field$type, "DictionaryType")

            values <- eval(parse(text=as.character(node[3])))
            if (arrow_type_name == "int32" && !is_enum) {
                values <- as.integer(values)
            }

            return(tiledbsoma_query_condition_in_nin(attr_name, tdb_op_name, values, somactx))

        } else if (.is_comparison_operator(node[1])) {
            spdl::debug("[parseqc] cmpop [{}] [{}] [{}]",
                as.character(node[2]),
                as.character(node[1]),
                as.character(node[3]))

            op_name <- as.character(node[1])
            attr_name <- as.character(node[2])
            rhs_text <- as.character(node[3])

            arrow_field <- schema[[attr_name]]
            if (is.null(arrow_field)) {
                .error_function("No attribute '", attr_name, "' is present.", call. = FALSE)
            }
            arrow_type_name <- arrow_field$type$name

            # Take care of factor (aka "enum" case) and set the data type to ASCII
            if (arrow_type_name == "dictionary") {
                arrow_type_name <- "utf8"
            }

            # General case of extracting appropriate value given type info
            return(tiledbsoma_query_condition_from_triple(
                attr_name = attr_name,
                value = switch(
                    arrow_type_name,
                    ascii = rhs_text,
                    utf8 = rhs_text,
                    bool = as.logical(rhs_text),
                    ## XXX DATETIME_MS = as.POSIXct(rhs_text),
                    ## XXX DATETIME_DAY = as.Date(rhs_text),
                    as.numeric(rhs_text)),
                arrow_type_name = arrow_type_name,
                op_name = .map_op_to_character(op_name),
                qc = tiledbsoma_empty_query_condition(somactx)))

        } else {
            stop("Unexpected token in expression: ", format(node))
        }
    }

    # Use base-r `substitute` to map the user-provided expression to a parse tree
    parse_tree <- substitute(expr)

    # Map the parse tree to TileDB core QueryCondition
    return(.parse_tree_to_qc(parse_tree, debug))
}

# ================================================================
#' An S4 class for a TileDB QueryCondition object
#'
#' @slot ptr An external pointer to the underlying implementation
#' @slot init A logical variable tracking if the query condition object has been
#' initialized
#' @exportClass tiledbsoma_query_condition
setClass(
    "tiledbsoma_query_condition",
    slots = list(ptr = "externalptr", init = "logical"))

# ================================================================
#' Creates a 'tiledbsoma_query_condition' object
#'
#' @param ctx (optional) A TileDB Ctx object; if not supplied the default
#' context object is retrieved
#' @return A 'tiledbsoma_query_condition' object
#' @export
tiledbsoma_empty_query_condition <- function(somactx) {
    stopifnot("The argument must be a ctx object" = is(ctx, "externalptr"))
    ptr <- libtiledbsoma_empty_query_condition(somactx)
    query_condition <- new("tiledbsoma_query_condition", ptr = ptr, init = FALSE)
    invisible(query_condition)
}

# ================================================================
#' Initialize a 'tiledbsoma_query_condition' object
#'
#' Initializes (and possibly allocates) a query condition object using a triplet of
#' attribute name, comparison value, and operator.  Six types of conditions are supported,
#' they all take a single scalar comparison argument and attribute to compare against.
#' At present only integer or numeric attribute comparisons are implemented.
#' @param attr_name A character value with the scheme attribute name
#' @param value A scalar value that the attribute is compared against
#' @param arrow_type_name A character value with the TileDB data type of the attribute column, for
#' example 'float' or 'int32'
#' @param op_name A character value with the comparison operation. This must be one of
#' 'LT', 'LE', 'GT', 'GE', 'EQ', 'NE'.
#' @param qc A 'tiledbsoma_query_condition' object to be initialized by this call.
#' @return The initialized 'tiledbsoma_query_condition' object
#' @export
tiledbsoma_query_condition_from_triple <- function(
    attr_name,
    value,
    arrow_type_name,
    op_name,
    qc) {

    stopifnot(
        "Argument 'qc' with query condition object required" = inherits(qc, "tiledbsoma_query_condition"),
        "Argument 'attr_name' must be character" = is.character(attr_name),
        "Argument 'value' must be of length one" = (
            is.vector(value) ||
            bit64::is.integer64(value) ||
            inherits(value, "POSIXt") ||
            inherits(value, "Date")) && all.equal(length(value),1),
        "Argument 'arrow_type_name' must be character" = is.character(arrow_type_name),
        "Argument 'op_name' must be character" = is.character(op_name))

    op_name <- match.arg(op_name, c("LT", "LE", "GT", "GE", "EQ", "NE"))
    # If arrow_type_name is int64 or uint64 but the class of value does not yet inherit from
    # integer64, cast.
    if (grepl("int64", arrow_type_name) && !inherits(value, "integer64")) {
        value <- bit64::as.integer64(value)
    }
    libtiledbsoma_query_condition_from_triple(qc@ptr, attr_name, value, arrow_type_name, op_name)
    qc@init <- TRUE
    invisible(qc)
}

# ================================================================
#' Combine two 'tiledbsoma_query_condition' objects
#'
#' Combines two query condition objects using a relatiional operator.
#'
#' @param lhs A 'tiledbsoma_query_condition' object on the left-hand side of the relation
#' @param rhs A 'tiledbsoma_query_condition' object on the right-hand side of the relation
#' @param op_name A character value with the relation, which must be one of 'AND', 'OR' or 'NOT'.
#' @param somactx SOMAContext pointer.
#' @return The combined 'tiledbsoma_query_condition' object
#' @export
tiledbsoma_query_condition_combine <- function(lhs, rhs, op_name, somactx) {
    stopifnot(
        "Argument 'lhs' must be a query condition object" = is(lhs, "tiledbsoma_query_condition"),
        "Argument 'rhs' must be a query condition object" = is(rhs, "tiledbsoma_query_condition"),
        "Argument 'op_name' must be a character" = is.character(op_name))
    op_name <- match.arg(op_name, c("AND", "OR", "NOT"))
    qc <- tiledbsoma_empty_query_condition(somactx)
    qc@ptr <- libtiledbsoma_query_condition_combine(lhs@ptr, rhs@ptr, op_name)
    qc@init <- TRUE
    invisible(qc)
}

# ================================================================
#' Create a query condition for vector 'IN' and 'NOT_IN' operations
#'
#' Uses \sQuote{IN} and \sQuote{NOT_IN} operators on given attribute
#'
#' @param attr_name A character value with the schema attribute name.
#'
#' @param op_name A character value with the chosen set operation. This must be one of
#' \sQuote{IN} or \sQuote{NOT_IN}.
#'
#' @param values A vector wiith the given values. Supported types are integer, double,
#' integer64, and character.
#'
#' @param somactx SOMAContext pointer.
#'
#' @return A query-condition object is returned
#' @export
tiledbsoma_query_condition_in_nin <- function(
  attr_name,
  op_name = "IN",
  values,
  somactx) {
    stopifnot("Argument 'attr_name' must be character" = is.character(attr_name),
              "Argument 'values' must be int, double, int64 or char" =
                  (is.numeric(values) || bit64::is.integer64(values) || is.character(values)),
              "Argument 'op_name' must be one of 'IN' or 'NOT_IN'" = op_name %in% c("IN", "NOT_IN"))

    qc <- tiledbsoma_empty_query_condition(somactx)
    qc@ptr <- libtiledbsoma_query_condition_in_nin(somactx, attr_name, op_name, values)
    qc@init <- TRUE
    invisible(qc)
}
