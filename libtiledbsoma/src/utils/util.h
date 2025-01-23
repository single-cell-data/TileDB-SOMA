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
 * @brief Take an arrow schema and array containing bool
 * data in bits and return a vector containing the uint8_t
 * representation
 *
 * @param schema the ArrowSchema which must be format 'b'
 * @param array the ArrowArray holding Boolean data
 * @return std::vector<uint8_t>
 */
std::vector<uint8_t> cast_bit_to_uint8(ArrowSchema* schema, ArrowArray* array);

std::shared_ptr<SOMAColumn> find_column_by_name(
    std::span<const std::shared_ptr<SOMAColumn>> columns,
    std::string_view name);

}  // namespace tiledbsoma::util

#endif
