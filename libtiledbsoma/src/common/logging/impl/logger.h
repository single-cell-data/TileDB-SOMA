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
 * This file defines class Logger, which is implemented as a wrapper around
 * `spdlog`. By policy `spdlog` must remain encapsulated as an implementation
 * and not be exposed as a dependency of the TileDB library. Accordingly, this
 * header should not be included as a header in any other header file. For
 * inclusion in a header (notably for use within the definition of
 * template-dependent functions), include the header `common/logging/logger.h`.
 *
 * The reason for this restriction is a technical limitation in template
 * instantiation. Part of the interface to `spdlog` consists of template
 * functions with variadic template arguments. Instantiation of such function
 * does not instantiate a variadic function (for exmaple `printf`) but rather a
 * function with a fixed number of arguments that depend upon the argument list.
 * Such variadic template argument lists cannot be forwarded across the
 * boundaries of compilation units, so exposing variadic template arguments
 * necessarily exposes the dependency upon `spdlog`. Thus this file `logger.h`,
 * which does have such arguments, must remain entirely within the library, but
 * `common/logging/logger.h`, which does not have such arguments, may be exposed without
 * creating an additional external dependency.
 */

#ifndef COMMON_LOGGER_IMPL_H
#define COMMON_LOGGER_IMPL_H

#include <spdlog/spdlog.h>

namespace tiledbsoma::common::logging::impl {
class Logger {
   public:
    Logger();
    ~Logger();

    /**
     * A formatted trace statment.
     *
     * @param fmt A fmtlib format string, see http://fmtlib.net/latest/ for
     *     details.
     * @param args optional positional arguments to format.
     */
    template <typename... Args>
    void trace(const char* fmt, const Args&... args) {
        logger_->trace(fmt, args...);
    }

    /**
     * A formatted debug statment.
     *
     * @param fmt A fmtlib format string, see http://fmtlib.net/latest/ for
     *     details.
     * @param args optional positional arguments to format.
     */
    template <typename... Args>
    void debug(const char* fmt, const Args&... args) {
        logger_->debug(fmt, args...);
    }

    /**
     * A formatted info statment.
     *
     * @param fmt A fmtlib format string, see http://fmtlib.net/latest/ for
     *     details.
     * @param args optional positional arguments to format.
     */
    template <typename... Args>
    void info(const char* fmt, const Args&... args) {
        logger_->info(fmt, args...);
    }

    /**
     * A formatted warn statment.
     *
     * @param fmt A fmtlib format string, see http://fmtlib.net/latest/ for
     *     details.
     * @param args optional positional arguments to format.
     */
    template <typename... Args>
    void warn(const char* fmt, const Args&... args) {
        logger_->warn(fmt, args...);
    }

    /** A formatted error statement.
     *
     * @param fmt A fmtlib format string, see http://fmtlib.net/latest/ for
     * details.
     * @param args optional positional arguments to format.
     */
    template <typename... Args>
    void error(const char* fmt, const Args&... args) {
        logger_->error(fmt, args...);
    }

    /**
     * A formatted critical statment.
     *
     * @param fmt A fmtlib format string, see http://fmtlib.net/latest/ for
     *     details.
     * @param args optional positional arguments to format.
     */
    template <typename... Args>
    void critical(const char* fmt, const Args&... args) {
        logger_->critical(fmt, args...);
    }

    /** Verbosity level. */
    enum class Level : char {
        FATAL,
        ERR,
        WARN,
        INFO,
        DBG,
        TRACE,
    };

    /**
     * Set the logger level.
     *
     * @param level log level string (FATAL|ERROR|WARN|INFO|DEBUG|TRACE)
     */
    void set_level(std::string_view level);

    /**
     * Set the logger output file.
     *
     * @param filename
     */
    void set_logfile(std::string_view filename);

    /**
     * Return true if debug messages are enabled.
     */
    bool debug_enabled();

   private:
    /* ********************************* */
    /*         PRIVATE ATTRIBUTES        */
    /* ********************************* */

    /** The logger object. */
    std::shared_ptr<spdlog::logger> logger_;
    spdlog::level::level_enum level_;
    std::string logfile_;
};
};  // namespace tiledbsoma::common::logging::impl

#endif