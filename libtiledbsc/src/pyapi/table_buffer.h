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

/**
 * @brief Class to store ColumnBuffers.
 *
 * TODO: Currently used for testing only, may be removed.
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
