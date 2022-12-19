// clang-format off
//                                               Emacs make this -*- mode: C++; -*-
//
// R and spdlog -- simplified interface borrowing from RcppSpdlog
//
// RcppSpdlog brings its own copies of spdlog and fmt
// This header defines a matching sink but assumming it is part
// of larger build that already has spdlog and fmt included

#pragma once

// define R_R_H and USING_R so that REprintf and R_FlushConsole are used
// by spdlog, and also declare them hereby
#define R_R_H
#define USING_R
extern "C" {
    void Rprintf(const char *, ...);
    void REprintf(const char *, ...);
    void R_FlushConsole(void);
}

// this define is important to not include another logger pulling in stdout
#define SPDLOG_DISABLE_DEFAULT_LOGGER 1

// include main header for spdlog -- others headers needed for extra features
#include "spdlog/spdlog.h"
#include "spdlog/sinks/base_sink.h"

#include <mutex>

namespace spdlog {
    namespace sinks {

        template<typename Mutex>
        class r_sink : public base_sink<Mutex> {

        protected:
            void sink_it_(const spdlog::details::log_msg& msg) override {
                // log_msg is a struct containing the log entry info like level, timestamp,
                // thread id etc.; msg.raw contains pre formatted log

                // If needed (very likely but not mandatory), the sink formats the message before
                // sending it to its final destination:
                spdlog::memory_buf_t formatted;
                spdlog::sinks::base_sink<Mutex>::formatter_->format(msg, formatted);
                Rprintf(fmt::to_string(formatted).c_str());
            }

            void flush_() override {
                R_FlushConsole();
            }
        };

        using r_sink_mt = r_sink<std::mutex>;
        using r_sink_st = r_sink<details::null_mutex>;

    } // namespace sinks

    // factory functions
    template<typename Factory = spdlog::synchronous_factory>
    inline std::shared_ptr<logger> r_sink_mt(const std::string &logger_name) {
        return Factory::template create<sinks::r_sink_mt>(logger_name);
    }

    template<typename Factory = spdlog::synchronous_factory>
    inline std::shared_ptr<logger> r_sink_st(const std::string &logger_name) {
        return Factory::template create<sinks::r_sink_st>(logger_name);
    }

} // namespace spdlog
