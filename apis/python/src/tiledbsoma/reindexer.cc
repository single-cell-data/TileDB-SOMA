/**
 * @file   indexer.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2024 TileDB, Inc.
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
 * This file defines the Reindexer bindings.
 */

#include <tiledbsoma/reindexer/reindexer.h>
// #include <tiledbsoma/utils/carrow.h>
#include "common.h"

#define DENUM(x) .value(#x, TILEDB_##x)
namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

/***
 * Handle general lookup for Re-indexer
 * @param indexer reference to the indexer
 * @param lookups input values to be looked up
 * @return looked up values
 */
py::array_t<int64_t> get_indexer_general_aux(
    IntIndexer& indexer, py::array_t<int64_t> lookups) {
    auto input_buffer = lookups.request();
    int64_t* input_ptr = static_cast<int64_t*>(input_buffer.ptr);
    size_t size = input_buffer.shape[0];
    auto results = py::array_t<int64_t>(size);
    auto results_buffer = results.request();
    int64_t* results_ptr = static_cast<int64_t*>(results_buffer.ptr);
    indexer.lookup(input_ptr, results_ptr, size);
    return results;
}
py::array_t<int64_t> get_indexer_general(
    IntIndexer& indexer, py::array_t<int64_t> lookups) {
    try {
        return get_indexer_general_aux(indexer, lookups);
    } catch (const std::exception& e) {
        throw TileDBSOMAError(e.what());
    }
}

/***
 * Helper function to provide data and schema for an arrow object
 * @param object python object
 * @param arrow_array extracted array data
 * @param arrow_schema extracted array schema
 */
//
void extract_py_array_schema(
    const pybind11::handle object,
    ArrowArray& arrow_array,
    ArrowSchema& arrow_schema) {
    uintptr_t arrow_schema_ptr = (uintptr_t)(&arrow_schema);
    uintptr_t arrow_array_ptr = (uintptr_t)(&arrow_array);

    // Using array._export_to_c to get arrow array and schema.
    object.attr("_export_to_c")(arrow_array_ptr, arrow_schema_ptr);
}

/***
 * Handle pyarrow-based lookup for Re-indexer
 * @param indexer reference to the indexer
 * @py_arrow_array pyarrow inputs to be looked up
 * @return looked up values
 */
py::array_t<int64_t> get_indexer_py_arrow_aux(
    IntIndexer& indexer, py::object py_arrow_array) {
    // Check if it is not a pyarrow array or pyarrow chunked array
    if (!py::hasattr(py_arrow_array, "_export_to_c") &&
        !py::hasattr(py_arrow_array, "chunks") &&
        !py::hasattr(py_arrow_array, "combine_chunks")) {
        // Handle the general case (no py arrow objects)
        return get_indexer_general(indexer, py_arrow_array);
    }

    py::list array_chunks;
    if (py::hasattr(py_arrow_array, "chunks")) {
        array_chunks = py_arrow_array.attr("chunks").cast<py::list>();
    } else {
        array_chunks.append(py_arrow_array);
    }

    // Calculate the total size of the input chunked array.
    int total_size = 0;
    for (const pybind11::handle array : array_chunks) {
        ArrowSchema arrow_schema;
        ArrowArray arrow_array;
        extract_py_array_schema(array, arrow_array, arrow_schema);
        total_size += arrow_array.length;

        arrow_schema.release(&arrow_schema);
        arrow_array.release(&arrow_array);
    }

    // Allocate the output
    auto results = py::array_t<int64_t>(total_size);
    auto results_buffer = results.request();
    int64_t* results_ptr = static_cast<int64_t*>(results_buffer.ptr);

    // Write output (one chunk at a time)
    int write_offset = 0;
    for (const pybind11::handle array : array_chunks) {
        ArrowSchema arrow_schema;
        ArrowArray arrow_array;
        extract_py_array_schema(array, arrow_array, arrow_schema);
        auto input_ptr = (int64_t*)arrow_array.buffers[1];
        indexer.lookup(
            input_ptr, results_ptr + write_offset, arrow_array.length);
        write_offset += arrow_array.length;

        arrow_schema.release(&arrow_schema);
        arrow_array.release(&arrow_array);
    }
    return results;
}

py::array_t<int64_t> get_indexer_py_arrow(
    IntIndexer& indexer, py::object py_arrow_array) {
    try {
        return get_indexer_py_arrow_aux(indexer, py_arrow_array);
    } catch (const std::exception& e) {
        throw TileDBSOMAError(e.what());
    }
}

void load_reindexer(py::module& m) {
    // Efficient C++ re-indexing (aka hashing unique key values to an index
    // between 0 and number of keys - 1) based on khash
    py::class_<IntIndexer>(m, "IntIndexer")
        .def(py::init<>())
        .def(py::init<std::shared_ptr<SOMAContext>>())
        .def(
            "map_locations",
            [](IntIndexer& indexer, py::array_t<int64_t> keys) {
                try {
                    indexer.map_locations(keys.data(), keys.size());
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })
        // Perform lookup for a large input array of keys and writes the
        // looked up values into previously allocated array (works for the
        // cases in which python and R pre-allocate the array)
        .def("get_indexer_general", get_indexer_general)
        // If the input is not arrow (does not have _export_to_c attribute),
        // it will be handled using a general input method.
        .def("get_indexer_pyarrow", get_indexer_py_arrow);
}

}  // namespace libtiledbsomacpp
