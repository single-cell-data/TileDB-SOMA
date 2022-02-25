#include <string>

#include <tiledbsc/common.h>

namespace tiledbsc {

TileDBSCError::TileDBSCError(const char* m) : std::runtime_error(m) {};
TileDBSCError::TileDBSCError(std::string m) : std::runtime_error(m.c_str()) {};

const char* TileDBSCError::what() const noexcept {
  return std::runtime_error::what();
}

}; // end namespace tiledbsc