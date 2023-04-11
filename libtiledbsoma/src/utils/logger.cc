/**
 * @file   logger.cc
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
 * This file defines class Logger, declared in logger.h, and the public logging
 * functions, declared in logger_public.h.
 */

#include "logger.h"

#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace tiledbsoma {

// Set the default logging format
// %^ : start color range
// [Year-month-day 24hr-min-second.microsecond]
// [logger]
// [Process: id]
// [Thread: id]
// [log level]
// text to log...
// %$ : end color range
const std::string LOG_PATTERN =
    "%^[%Y-%m-%d %H:%M:%S.%e] [%n] [Process: %P] [Thread: %t] [%l] %v%$";
const std::string CONSOLE_LOGGER = "tiledbsoma";
const std::string FILE_LOGGER = "tiledbsoma-file";

/* ********************************* */
/*     CONSTRUCTORS & DESTRUCTORS    */
/* ********************************* */

Logger::Logger() {
    logger_ = spdlog::get(CONSOLE_LOGGER);
    if (logger_ == nullptr) {
#if !defined(R_BUILD)
        logger_ = spdlog::stdout_color_mt(CONSOLE_LOGGER);
#else
        logger_ = spdlog::r_sink_mt(CONSOLE_LOGGER);
        spdlog::set_default_logger(logger_);
#endif
#if !defined(R_BUILD) && !defined(_WIN32)
        // change color of critical messages
        auto console_sink = static_cast<spdlog::sinks::stdout_color_sink_mt*>(
            logger_->sinks().back().get());
        console_sink->set_color(
            spdlog::level::critical, console_sink->red_bold);
#endif
    }
    logger_->set_pattern(LOG_PATTERN);
    set_level("INFO");
}

Logger::~Logger() {
    spdlog::drop(CONSOLE_LOGGER);
    if (spdlog::get(FILE_LOGGER) != nullptr) {
        spdlog::drop(FILE_LOGGER);
    }
}

void Logger::trace(const char* msg) {
    logger_->trace(msg);
}

void Logger::debug(const char* msg) {
    logger_->debug(msg);
}

void Logger::info(const char* msg) {
    logger_->info(msg);
}

void Logger::warn(const char* msg) {
    logger_->warn(msg);
}

void Logger::error(const char* msg) {
    logger_->error(msg);
}

void Logger::critical(const char* msg) {
    logger_->critical(msg);
}

void Logger::set_level(const std::string& level_in) {
    // convert level to lower case
    std::string level = level_in;
    std::for_each(
        level.begin(), level.end(), [](char& c) { c = ::tolower(c); });

    if (level == "fatal" || level[0] == 'f') {
        level_ = spdlog::level::critical;
    } else if (level == "error" || level[0] == 'e') {
        level_ = spdlog::level::err;
    } else if (level == "warn" || level[0] == 'w') {
        level_ = spdlog::level::warn;
    } else if (level == "info" || level[0] == 'i') {
        level_ = spdlog::level::info;
    } else if (level == "debug" || level[0] == 'd') {
        level_ = spdlog::level::debug;
    } else if (level == "trace" || level[0] == 't') {
        level_ = spdlog::level::trace;
    } else {
        set_level("WARN");
        //        LOG_WARN("Illegal log level = {}, using log level FATAL",
        //        level);
        level_ = spdlog::level::critical;
    }
    logger_->set_level(level_);
}

void Logger::set_logfile(const std::string& filename) {
    if (!logfile_.empty()) {
        // LOG_WARN("Already logging messages to {}", logfile_);
        return;
    }

    logfile_ = filename;

    try {
        auto file_logger = spdlog::basic_logger_mt(FILE_LOGGER, filename);
        file_logger->set_pattern(LOG_PATTERN);
        file_logger->set_level(level_);
    } catch (spdlog::spdlog_ex& e) {
        // log message and exit if file logger cannot be created
        LOG_FATAL(e.what());
    }

    // add sink to existing logger
    // (https://github.com/gabime/spdlog/wiki/4.-Sinks)
    auto file_sink = spdlog::get(FILE_LOGGER)->sinks().back();
    logger_->sinks().push_back(file_sink);
    logger_->flush_on(spdlog::level::info);
}

bool Logger::debug_enabled() {
    return (level_ == spdlog::level::debug) || (level_ == spdlog::level::trace);
}

/* ********************************* */
/*              GLOBAL               */
/* ********************************* */

Logger& global_logger() {
    static Logger l;
    return l;
}

/** Set log level for global logger. */
void LOG_CONFIG(const std::string& level, const std::string& logfile) {
    if (!level.empty()) {
        global_logger().set_level(level);
    }
    if (!logfile.empty()) {
        global_logger().set_logfile(logfile);
    }
}

/** Set log level for global logger. */
void LOG_SET_LEVEL(const std::string& level) {
    global_logger().set_level(level);
}

/** Set log file for global logger. */
void LOG_SET_FILE(const std::string& logfile) {
    global_logger().set_logfile(logfile);
}

/** Check if global logger is logging debug messages. */
bool LOG_DEBUG_ENABLED() {
    return global_logger().debug_enabled();
}

/** Logs a trace message. */
void LOG_TRACE(const std::string& msg) {
    global_logger().trace(msg.c_str());
}

/** Logs a debug message. */
void LOG_DEBUG(const std::string& msg) {
    global_logger().debug(msg.c_str());
}

/** Logs an info message. */
void LOG_INFO(const std::string& msg) {
    global_logger().info(msg.c_str());
}

/** Logs a warning. */
void LOG_WARN(const std::string& msg) {
    global_logger().warn(msg.c_str());
}

/** Logs an error. */
void LOG_ERROR(const std::string& msg) {
    global_logger().error(msg.c_str());
}

/** Logs a critical error and exits with a non-zero status. */
void LOG_FATAL(const std::string& msg) {
    global_logger().critical(msg.c_str());
#if !defined(R_BUILD)
    exit(1);
#endif
}

/** Convert TileDB timestamp (in ms) to human readable timestamp. */
std::string asc_timestamp(uint64_t timestamp_ms) {
    auto time_sec = static_cast<time_t>(timestamp_ms) / 1000;
    std::string time_str = asctime(gmtime(&time_sec));
    time_str.pop_back();  // remove newline
    time_str += " UTC";
    return time_str;
}

}  // namespace tiledbsoma
