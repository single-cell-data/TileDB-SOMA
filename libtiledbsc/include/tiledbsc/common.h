#include <stdexcept>

#include "tiledbsc_export.h"

namespace tiledbsc {

class TILEDBSC_EXPORT TileDBSCError : std::runtime_error {
public:
  explicit TileDBSCError(const char *m);
  explicit TileDBSCError(std::string m);

public:
  virtual const char *what() const noexcept override;
};

}; // namespace tiledbsc