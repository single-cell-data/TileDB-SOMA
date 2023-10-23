/**
 * @file   soma_sparse_ndarray.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2023 TileDB, Inc.
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
 * This file defines the SOMASparseNDArray bindings.
 */

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsoma/tiledbsoma>

#include "common.h"

using namespace tiledbsoma;

namespace py = pybind11;

namespace tiledbsoma {
void load_soma_sparse_ndarray(py::module &m) {
    py::class_<SOMASparseNDArray>(m, "SOMASparseNDArray")

    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::map<std::string, std::string>, std::vector<std::string>, ResultOrder, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMASparseNDArray::open))
    .def_static("exists", &SOMASparseNDArray::exists)
    .def("reopen", py::overload_cast<OpenMode, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMASparseNDArray::open))
    .def("close", &SOMASparseNDArray::close)
    .def_property_readonly("closed", [](
        SOMASparseNDArray& soma_sparse_ndarr) -> bool { 
        return soma_sparse_ndarr.is_open();
    })
    .def("type", &SOMASparseNDArray::type)
    .def("ctx", &SOMASparseNDArray::ctx)
    .def("is_sparse", &SOMASparseNDArray::is_sparse)
    .def("uri", &SOMASparseNDArray::uri)
    .def("schema", &SOMASparseNDArray::schema)
    .def_property_readonly("schema", [](SOMASparseNDArray& soma_sparse_ndarr) -> py::object {
        auto pa = py::module::import("pyarrow");
        auto pa_schema_import = pa.attr("Schema").attr("_import_from_c");
        return pa_schema_import(py::capsule(soma_sparse_ndarr.schema().get()));
    })
    .def_property_readonly("timestamp", &SOMASparseNDArray::timestamp)
    .def_property_readonly("index_column_names", &SOMASparseNDArray::index_column_names)
        .def("non_empty_domain", [](SOMASparseNDArray& soma_sparse_ndarr, std::string name, py::dtype dtype){
        switch (np_to_tdb_dtype(dtype)) {
        case TILEDB_UINT64:
            return py::cast(soma_sparse_ndarr.non_empty_domain<uint64_t>(name));
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS:
        case TILEDB_INT64:
            return py::cast(soma_sparse_ndarr.non_empty_domain<int64_t>(name));
        case TILEDB_UINT32:
            return py::cast(soma_sparse_ndarr.non_empty_domain<uint32_t>(name));
        case TILEDB_INT32:
            return py::cast(soma_sparse_ndarr.non_empty_domain<int32_t>(name));
        case TILEDB_UINT16:
            return py::cast(soma_sparse_ndarr.non_empty_domain<uint16_t>(name));
        case TILEDB_INT16:
            return py::cast(soma_sparse_ndarr.non_empty_domain<int16_t>(name));
        case TILEDB_UINT8:
            return py::cast(soma_sparse_ndarr.non_empty_domain<uint8_t>(name));
        case TILEDB_INT8:
            return py::cast(soma_sparse_ndarr.non_empty_domain<int8_t>(name));
        case TILEDB_FLOAT64:
            return py::cast(soma_sparse_ndarr.non_empty_domain<double>(name));
        case TILEDB_FLOAT32:
            return py::cast(soma_sparse_ndarr.non_empty_domain<float>(name));
        case TILEDB_STRING_UTF8:
        case TILEDB_STRING_ASCII: 
            return py::cast(soma_sparse_ndarr.non_empty_domain_var(name));
        default:
            throw TileDBSOMAError("Unsupported dtype for nonempty domain.");
        }
    })
    .def("domain", [](SOMASparseNDArray& soma_sparse_ndarr, std::string name, py::dtype dtype) {
        switch (np_to_tdb_dtype(dtype)) {
        case TILEDB_UINT64:
            return py::cast(soma_sparse_ndarr.domain<uint64_t>(name));
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS:
        case TILEDB_INT64:
            return py::cast(soma_sparse_ndarr.domain<int64_t>(name));
        case TILEDB_UINT32:
            return py::cast(soma_sparse_ndarr.domain<uint32_t>(name));
        case TILEDB_INT32:
            return py::cast(soma_sparse_ndarr.domain<int32_t>(name));
        case TILEDB_UINT16:
            return py::cast(soma_sparse_ndarr.domain<uint16_t>(name));
        case TILEDB_INT16:
            return py::cast(soma_sparse_ndarr.domain<int16_t>(name));
        case TILEDB_UINT8:
            return py::cast(soma_sparse_ndarr.domain<uint8_t>(name));
        case TILEDB_INT8:
            return py::cast(soma_sparse_ndarr.domain<int8_t>(name));
        case TILEDB_FLOAT64:
            return py::cast(soma_sparse_ndarr.domain<double>(name));
        case TILEDB_FLOAT32:
            return py::cast(soma_sparse_ndarr.domain<float>(name));
        case TILEDB_STRING_UTF8:
        case TILEDB_STRING_ASCII: {
            std::pair<std::string, std::string> str_domain;
            return py::cast(std::make_pair("", ""));
        }
        default:
            throw TileDBSOMAError("Unsupported dtype for Dimension's domain");
        }
    })
    .def("shape", &SOMASparseNDArray::shape)
    .def("ndim", &SOMASparseNDArray::ndim)
    .def("nnz", &SOMASparseNDArray::nnz)
    .def("read_next", &SOMASparseNDArray::read_next)
    .def("write", &SOMASparseNDArray::write)
    .def("set_metadata", &SOMASparseNDArray::set_metadata)
    .def("delete_metadata", &SOMASparseNDArray::delete_metadata)
    .def_property_readonly("meta", [](SOMASparseNDArray& soma_sparse_ndarr) -> py::dict {
        py::dict results;
            
        for (auto const& [key, val] : soma_sparse_ndarr.get_metadata()){
            tiledb_datatype_t tdb_type = std::get<MetadataInfo::dtype>(val);
            uint32_t value_num = std::get<MetadataInfo::num>(val);
            const void *value = std::get<MetadataInfo::value>(val);

            if(tdb_type == TILEDB_STRING_UTF8){
                results[py::str(key)] = py::str(std::string((const char*)value, value_num));
            }else if(tdb_type == TILEDB_STRING_ASCII){
                results[py::str(key)] = py::bytes(std::string((const char*)value, value_num));
            }else{
                py::dtype value_type = tdb_to_np_dtype(tdb_type, 1);
                results[py::str(key)] = py::array(value_type, value_num, value);
            }
        }
        return results;
    })
    .def("get_metadata", 
        py::overload_cast<const std::string&>(&SOMASparseNDArray::get_metadata))
    .def("get_metadata", py::overload_cast<>(&SOMASparseNDArray::get_metadata))
    .def("has_metadata", &SOMASparseNDArray::has_metadata)
    .def("metadata_num", &SOMASparseNDArray::metadata_num);
}
}