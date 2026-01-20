/**
 * @file   logger.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines class Logger, declared in `logger.h`, and the public logging
 * functions, declared in `common/logging/logger.h`.
 */

#include "logger.h"

#include <spdlog/cfg/env.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace tiledbsoma::common::logging {

namespace impl {

bool sv_compare(std::string_view lhs, std::string_view rhs) {
    return std::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), [](const char& l, const char& r) {
        return ::tolower(l) == ::tolower(r);
    });
}

// Set the default logging format
// %^ : start color range
// [Year-month-day 24hr-min-second.microsecond]
// [logger]
// [Process: id]
// [Thread: id]
// [log level]
// text to log...
// %$ : end color range
inline constexpr std::string_view LOG_PATTERN{"%^[%Y-%m-%d %H:%M:%S.%e] [%n] [Process: %P] [Thread: %t] [%l] %v%$"};
inline constexpr std::string_view CONSOLE_LOGGER{"tiledbsoma"};
inline constexpr std::string_view FILE_LOGGER{"tiledbsoma-file"};

/* ********************************* */
/*     CONSTRUCTORS & DESTRUCTORS    */
/* ********************************* */

Logger::Logger() {
    logger_ = spdlog::get(CONSOLE_LOGGER.data());
    if (logger_ == nullptr) {
        logger_ = spdlog::stdout_color_mt(CONSOLE_LOGGER.data());
        logger_->set_pattern(LOG_PATTERN.data());
#if !defined(_WIN32)
        // change color of critical messages
        auto console_sink = static_cast<spdlog::sinks::stdout_color_sink_mt*>(logger_->sinks().back().get());
        console_sink->set_color(spdlog::level::critical, console_sink->red_bold);
#endif
    }
    set_level("INFO");
    // Examples:
    // SPDLOG_LEVEL=trace name-of-program
    // SPDLOG_LEVEL=tiledbsoma=trace name-of-program
    spdlog::cfg::load_env_levels();
}

Logger::~Logger() {
    spdlog::drop(CONSOLE_LOGGER.data());
    if (spdlog::get(FILE_LOGGER.data()) != nullptr) {
        spdlog::drop(FILE_LOGGER.data());
    }
}

void Logger::set_level(std::string_view level) {
    if (sv_compare(level, "fatal") || level[0] == 'f') {
        level_ = spdlog::level::critical;
    } else if (sv_compare(level, "error")) {
        level_ = spdlog::level::err;
    } else if (sv_compare(level, "warn")) {
        level_ = spdlog::level::warn;
    } else if (sv_compare(level, "info")) {
        level_ = spdlog::level::info;
    } else if (sv_compare(level, "debug")) {
        level_ = spdlog::level::debug;
    } else if (sv_compare(level, "trace")) {
        level_ = spdlog::level::trace;
    } else {
        level_ = spdlog::level::critical;
    }
    logger_->set_level(level_);
}

void Logger::set_logfile(std::string_view filename) {
    if (!logfile_.empty()) {
        // LOG_WARN("Already logging messages to {}", logfile_);
        return;
    }

    logfile_ = filename;

    try {
        auto file_logger = spdlog::basic_logger_mt(FILE_LOGGER.data(), filename.data());
        file_logger->set_pattern(LOG_PATTERN.data());
        file_logger->set_level(level_);
    } catch (spdlog::spdlog_ex& e) {
        // log message and exit if file logger cannot be created
        logger_->error(e.what());
    }

    // add sink to existing logger
    // (https://github.com/gabime/spdlog/wiki/4.-Sinks)
    auto file_sink = spdlog::get(FILE_LOGGER.data())->sinks().back();
    logger_->sinks().push_back(file_sink);
    logger_->flush_on(spdlog::level::info);
}

bool Logger::debug_enabled() {
    return (level_ == spdlog::level::debug) || (level_ == spdlog::level::trace);
}

}  // namespace impl
/* ********************************* */
/*              GLOBAL               */
/* ********************************* */

impl::Logger& global_logger() {
    static impl::Logger l;
    return l;
}

/** Set log level for global logger. */
void LOG_CONFIG(std::string_view level, std::string_view logfile) {
    if (!level.empty()) {
        global_logger().set_level(level);
    }
    if (!logfile.empty()) {
        global_logger().set_logfile(logfile);
    }
}

/** Set log level for global logger. */
void LOG_SET_LEVEL(std::string_view level) {
    global_logger().set_level(level);
}

/** Set log file for global logger. */
void LOG_SET_FILE(std::string_view logfile) {
    global_logger().set_logfile(logfile);
}

/** Check if global logger is logging debug messages. */
bool LOG_DEBUG_ENABLED() {
    return global_logger().debug_enabled();
}

/** Logs a trace message. */
void LOG_TRACE(std::string_view msg) {
    global_logger().trace(msg.data());
}

/** Logs a debug message. */
void LOG_DEBUG(std::string_view msg) {
    global_logger().debug(msg.data());
}

/** Logs an info message. */
void LOG_INFO(std::string_view msg) {
    global_logger().info(msg.data());
}

/** Logs a warning. */
void LOG_WARN(std::string_view msg) {
    global_logger().warn(msg.data());
}

/** Logs an error. */
void LOG_ERROR(std::string_view msg) {
    global_logger().error(msg.data());
}

/** Logs a critical error and exits with a non-zero status. */
void LOG_FATAL(std::string_view msg) {
    global_logger().critical(msg.data());
    exit(1);
}

/** Convert TileDB timestamp (in ms) to human readable timestamp. */
std::string asc_timestamp(uint64_t timestamp_ms) {
    auto time_sec = static_cast<time_t>(timestamp_ms) / 1000;
    std::string time_str = asctime(gmtime(&time_sec));
    time_str.pop_back();  // remove newline
    time_str += " UTC";
    return time_str;
}

}  // namespace tiledbsoma::common::logging
