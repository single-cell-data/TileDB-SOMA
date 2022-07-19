#ifndef TABLE_BUFFER_H
#define TABLE_BUFFER_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <span>
#include <tiledb/tiledb>

#include "tiledbsc/column_buffer.h"
#include "tiledbsc/common.h"
#include "tiledbsc/logger_public.h"

namespace tiledbsc {

using namespace tiledb;

using ColumnBuffers = std::map<std::string, std::shared_ptr<ColumnBuffer>>;

/**
 * @brief Class to store data for a TileDB dimension or attribute.
 *
 */
class TableBuffer {
   public:
    //===================================================================
    //= public non-static
    //===================================================================
    TableBuffer() = delete;

    TableBuffer(ColumnBuffers&& columns)
        : columns_(std::move(columns)){};

    ColumnBuffers& get_columns() {
        return columns_;
    }

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    ColumnBuffers columns_;
};

}  // namespace tiledbsc
#endif
