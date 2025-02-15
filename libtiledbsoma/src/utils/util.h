/**
 * @file   util.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the utility functions
 */

#ifndef UTIL_H
#define UTIL_H

#include <regex>
#include <span>
#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include "../soma/soma_column.h"
#include "common.h"
#include "nanoarrow/nanoarrow.hpp"

namespace tiledbsoma::util {

using VarlenBufferPair = std::pair<std::string, std::vector<uint64_t>>;

template <typename T>
VarlenBufferPair to_varlen_buffers(std::vector<T> data, bool arrow = true);

template <class T>
std::vector<T> to_vector(const std::span<T>& s) {
    return std::vector<T>(s.begin(), s.end());
}

/**
 * @brief Check if the provided URI is a TileDB Cloud URI.
 *
 * @param uri URI to check
 * @return true URI is a TileBD Cloud URI
 * @return false URI is not a TileBD Cloud URI
 */
bool is_tiledb_uri(std::string_view uri);

/**
 * @brief Remove all trailing '/' from URI.
 *
 * @param uri URI
 * @return std::string URI without trailing '/'
 */
std::string rstrip_uri(std::string_view uri);

/**
 * @brief Take a bitmap and return a vector containing the uint8_t
 * representation. If the bitmap is null, then return a nullopt.
 *
 * @param bitmap Pointer to the start of the bitmap
 * @param length Total number of elements
 * @param offset Optionally offset the data
 * @return std::optional<std::vector<uint8_t>>
 */
std::optional<std::vector<uint8_t>> bitmap_to_uint8(
    const uint8_t* bitmap, size_t length, size_t offset = 0);

std::shared_ptr<SOMAColumn> find_column_by_name(
    std::span<const std::shared_ptr<SOMAColumn>> columns,
    std::string_view name);

/**
 * @brief Construct the enumeration label given a dictionary's index and value
 * ArrowSchemas. The name is extracted from the index schema. The datatype is
 * extracted from the value schema. If the format string is u or z (regular
 * string or regular binary), then convert that to U or Z (large string or large
 * binary).
 *
 * @param schema ArrowSchema representing the dictionary values
 * @return enmr dtype as Arrow format string
 */
std::string get_enmr_label(
    ArrowSchema* index_schema, ArrowSchema* value_schema);

}  // namespace tiledbsoma::util

#endif
