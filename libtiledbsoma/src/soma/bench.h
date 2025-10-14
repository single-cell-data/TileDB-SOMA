/**
 * @file   column_buffer.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This declares the column buffer API
 */

#ifndef BENCH_H
#define BENCH_H

#include <span>
#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include "../utils/arrow_adapter.h"

namespace tiledbsoma {

using namespace tiledb;

class MemoryBench {
   public:
    //===================================================================
    //= public static
    //===================================================================

    static void release_schema(ArrowSchema* schema);
    static void release_vector_array(ArrowArray* array);
    static void release_pointer_array(ArrowArray* array);
    static void release_column_buffer_array(ArrowArray* array);

    static ArrowTable allocate_vector(uint64_t size);

    static ArrowTable allocate_pointer(uint64_t size);

    static ArrowTable allocate_column_buffer(uint64_t size);
};

}  // namespace tiledbsoma
#endif
