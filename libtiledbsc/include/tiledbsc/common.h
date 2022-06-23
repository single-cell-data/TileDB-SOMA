#ifndef TILEDBSC_COMMON_H
#define TILEDBSC_COMMON_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <string>

namespace tiledbsc {

class TileDBSCError : public std::runtime_error {
   public:
    explicit TileDBSCError(const char* m)
        : std::runtime_error(m){};
    explicit TileDBSCError(std::string m)
        : std::runtime_error(m.c_str()){};

   public:
    virtual const char* what() const noexcept override {
        return std::runtime_error::what();
    };
};

};  // namespace tiledbsc

#endif  // TILEDBSC_COMMON_H
