/**
 * @file   logger_public.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * This file defines simple logging functions that can be exposed (by expedient)
 * to the public API. Their implementations are in `logger.cc`. See the
 * documentation in `logger.h` for the full story.
 */

#pragma once
#ifndef TILEDB_LOGGER_PUBLIC_H
#define TILEDB_LOGGER_PUBLIC_H

#ifndef TILEDBSOMA_EXPORT
// TILEDBSOMA_EXPORT is defined by the auto-generated tiledbsoma_export.h
// which is included from the top-level tiledbsoma header. We don't include
// it here to simplify the include paths (and avoid copying all headers to
// a single directory like TileDB does). In case the symbol is not defined,
// define it here as empty.
#define TILEDBSOMA_EXPORT
#endif

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

namespace tiledbsoma {

/** Set log level for global logger and optionally set a logfile. */
TILEDBSOMA_EXPORT void LOG_CONFIG(
    const std::string& level, const std::string& logfile = "");

/** Set log level for global logger. */
void LOG_SET_LEVEL(const std::string& level);

/** Set log file for global logger. */
void LOG_SET_FILE(const std::string& logfile);

/** Check if global logger is logging debug messages. */
bool LOG_DEBUG_ENABLED();

/** Logs a trace message. */
TILEDBSOMA_EXPORT void LOG_TRACE(const std::string& msg);

/** Logs a debug message. */
TILEDBSOMA_EXPORT void LOG_DEBUG(const std::string& msg);

/** Logs an info message. */
TILEDBSOMA_EXPORT void LOG_INFO(const std::string& msg);

/** Logs a warning. */
TILEDBSOMA_EXPORT void LOG_WARN(const std::string& msg);

/** Logs an error. */
void LOG_ERROR(const std::string& msg);

/** Logs a critical error and exits with a non-zero status. */
void LOG_FATAL(const std::string& msg);

/** Convert TileDB timestamp (in ms) to human readable timestamp. */
std::string asc_timestamp(uint64_t timestamp_ms);

}  // namespace tiledbsoma

#endif  // TILEDB_LOGGER_PUBLIC_H
