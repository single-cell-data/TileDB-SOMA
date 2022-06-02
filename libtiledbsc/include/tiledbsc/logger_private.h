/**
 * @file   logger_private.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2017-2022 TileDB, Inc.
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
 * This file defines simple logging functions that *cannot* be exposed to the
 * public API. Their implementations are in `logger.cc`. See the documentation
 * in `logger.h` for the full story.
 */

#pragma once
#ifndef TILEDB_LOGGER_PRIVATE_H
#define TILEDB_LOGGER_PRIVATE_H

#include "logger.h"

namespace tiledbsc {

/** Set log level for global logger and optionally set a logfile. */
void LOG_CONFIG(const std::string& level, const std::string& logfile = "");

/** Set log level for global logger. */
void LOG_SET_LEVEL(const std::string& level);

/** Set log file for global logger. */
void LOG_SET_FILE(const std::string& logfile);

/** Check if global logger is logging debug messages. */
bool LOG_DEBUG_ENABLED();

/** Logs a trace message. */
void LOG_TRACE(const std::string& msg);

/** Logs a formatted trace message. */
template <typename Arg1, typename... Args>
void LOG_TRACE(const char* fmt, const Arg1& arg1, const Args&... args) {
    global_logger().trace(fmt, arg1, args...);
}

/** Logs a debug message. */
void LOG_DEBUG(const std::string& msg);

/** Logs a formatted debug message. */
template <typename Arg1, typename... Args>
void LOG_DEBUG(const char* fmt, const Arg1& arg1, const Args&... args) {
    global_logger().debug(fmt, arg1, args...);
}

/** Logs an info message. */
void LOG_INFO(const std::string& msg);

/** Logs a formatted info message. */
template <typename Arg1, typename... Args>
void LOG_INFO(const char* fmt, const Arg1& arg1, const Args&... args) {
    global_logger().info(fmt, arg1, args...);
}

/** Logs a warning. */
void LOG_WARN(const std::string& msg);

/** Logs a formatted warning message. */
template <typename Arg1, typename... Args>
void LOG_WARN(const char* fmt, const Arg1& arg1, const Args&... args) {
    global_logger().warn(fmt, arg1, args...);
}

/** Logs an error. */
void LOG_ERROR(const std::string& msg);

/** Logs a formatted error message. */
template <typename Arg1, typename... Args>
void LOG_ERROR(const char* fmt, const Arg1& arg1, const Args&... args) {
    global_logger().error(fmt, arg1, args...);
}

/** Logs a critical error and exits with a non-zero status. */
void LOG_FATAL(const std::string& msg);

/** Logs a formatted critical error and exits with a non-zero status. */
template <typename Arg1, typename... Args>
void LOG_FATAL(const char* fmt, const Arg1& arg1, const Args&... args) {
    global_logger().critical(fmt, arg1, args...);
    exit(1);
}

/** Convert TileDB timestamp (in ms) to human readable timestamp. */
std::string asc_timestamp(uint64_t timestamp_ms);

}  // namespace tiledbsc

#endif
