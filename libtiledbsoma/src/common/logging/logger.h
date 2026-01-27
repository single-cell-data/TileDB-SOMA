/**
 * @file   logger.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines simple logging functions that can be exposed (by expedient)
 * to the public API. Their implementations are in `impl/logger.cc`. See the
 * documentation in `impl/logger.h` for the full story.
 */

#ifndef COMMON_LOGGER_H
#define COMMON_LOGGER_H

#include <memory>
#include <stdexcept>

namespace tiledbsoma::common::logging {

/** Set log level for global logger and optionally set a logfile. */
void LOG_CONFIG(std::string_view level, std::string_view logfile = "");

/** Set log level for global logger. */
void LOG_SET_LEVEL(std::string_view level);

/** Set log file for global logger. */
void LOG_SET_FILE(std::string_view logfile);

/** Check if global logger is logging debug messages. */
bool LOG_DEBUG_ENABLED();

/** Logs a trace message. */
void LOG_TRACE(std::string_view msg);

/** Logs a debug message. */
void LOG_DEBUG(std::string_view msg);

/** Logs an info message. */
void LOG_INFO(std::string_view msg);

/** Logs a warning. */
void LOG_WARN(std::string_view msg);

/** Logs an error. */
void LOG_ERROR(std::string_view msg);

/** Logs a critical error and exits with a non-zero status. */
void LOG_FATAL(std::string_view msg);

/** Convert TileDB timestamp (in ms) to human readable timestamp. */
std::string asc_timestamp(uint64_t timestamp_ms);

}  // namespace tiledbsoma::common::logging

#endif  // COMMON_LOGGER_H
