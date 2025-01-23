/**
 * @file   indexer.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the Reindexer bindings.
 */

#include <tiledbsoma/reindexer/reindexer.h>
#include "common.h"

#define DENUM(x) .value(#x, TILEDB_##x)
namespace libtiledbsomacpp {

namespace py = pybind11;
using npy_api = pybind11::detail::npy_api;
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

    py::gil_scoped_release release;
    indexer.lookup(input_ptr, results_ptr, size);
    py::gil_scoped_acquire acquire;

    return results;
}
py::array_t<int64_t> get_indexer_general(
    IntIndexer& indexer, py::array_t<int64_t> lookups) {
    if (lookups.ndim() != 1) {
        throw std::invalid_argument(
            "IntIndexer only supports arrays of dimension 1");
    }
    if (!lookups.dtype().is(py::dtype::of<int64_t>())) {
        throw py::type_error("IntIndexer only supports array of type int64");
    }

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

        bool type_ok = (strcmp(arrow_schema.format, "l") == 0);

        arrow_schema.release(&arrow_schema);
        arrow_array.release(&arrow_array);

        if (!type_ok)
            throw TileDBSOMAError(
                "IntIndexer only supports array of type int64");
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
        auto input_ptr = (int64_t*)arrow_array.buffers[1] + arrow_array.offset;

        py::gil_scoped_release release;
        indexer.lookup(
            input_ptr, results_ptr + write_offset, arrow_array.length);
        py::gil_scoped_acquire acquire;

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
        .def(
            py::init<std::shared_ptr<SOMAContext>>(),
            py::arg("context").noconvert())
        .def(
            "map_locations",
            [](IntIndexer& indexer, py::array keys) {
                if (keys.ndim() != 1) {
                    throw std::invalid_argument(
                        "IntIndexer only supports arrays of dimension 1");
                }
                if (!keys.dtype().is(py::dtype::of<int64_t>())) {
                    throw py::type_error(
                        "IntIndexer only supports array of type int64");
                }

                auto keys_int64 = py::cast<py::array_t<int64_t>>(keys);
                try {
                    auto keys_p = keys_int64.data();
                    auto keys_sz = keys_int64.size();

                    py::gil_scoped_release release;
                    indexer.map_locations(keys_p, keys_sz);

                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        // Perform lookup for a large input array of keys and writes the
        // looked up values into previously allocated array (works for the
        // cases in which python and R pre-allocate the array)
        .def("get_indexer_general", get_indexer_general, py::arg().noconvert())

        // If the input is not arrow (does not have _export_to_c attribute),
        // it will be handled using a general input method.
        .def(
            "get_indexer_pyarrow", get_indexer_py_arrow, py::arg().noconvert());
}

}  // namespace libtiledbsomacpp
