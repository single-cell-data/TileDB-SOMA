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
void load_soma_dense_ndarray(py::module &m) {
    py::class_<SOMADenseNDArray>(m, "SOMADenseNDArray")

    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::map<std::string, std::string>, std::vector<std::string>, ResultOrder, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADenseNDArray::open))
    .def_static("exists", &SOMADenseNDArray::exists)
    .def("reopen", py::overload_cast<OpenMode, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADenseNDArray::open))
    .def("close", &SOMADenseNDArray::close)
    .def("is_open", &SOMADenseNDArray::is_open)
    .def("type", &SOMADenseNDArray::type)
    .def("ctx", &SOMADenseNDArray::ctx)
    .def("is_sparse", &SOMADenseNDArray::is_sparse)
    .def("uri", &SOMADenseNDArray::uri)
    .def("schema", &SOMADenseNDArray::schema)
    .def_property_readonly("schema", [](SOMADenseNDArray& soma_dense_ndarr) -> py::object {
        auto pa = py::module::import("pyarrow");
        auto pa_schema_import = pa.attr("Schema").attr("_import_from_c");
        return pa_schema_import(py::capsule(soma_dense_ndarr.schema().get()));
    })
    .def_property_readonly("timestamp", &SOMADenseNDArray::timestamp)
    .def_property_readonly("index_column_names", &SOMADenseNDArray::index_column_names)
        .def("non_empty_domain", [](SOMADenseNDArray& soma_dense_ndarr, std::string name, py::dtype dtype){
        switch (np_to_tdb_dtype(dtype)) {
        case TILEDB_UINT64:
            return py::cast(soma_dense_ndarr.non_empty_domain<uint64_t>(name));
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
            return py::cast(soma_dense_ndarr.non_empty_domain<int64_t>(name));
        case TILEDB_UINT32:
            return py::cast(soma_dense_ndarr.non_empty_domain<uint32_t>(name));
        case TILEDB_INT32:
            return py::cast(soma_dense_ndarr.non_empty_domain<int32_t>(name));
        case TILEDB_UINT16:
            return py::cast(soma_dense_ndarr.non_empty_domain<uint16_t>(name));
        case TILEDB_INT16:
            return py::cast(soma_dense_ndarr.non_empty_domain<int16_t>(name));
        case TILEDB_UINT8:
            return py::cast(soma_dense_ndarr.non_empty_domain<uint8_t>(name));
        case TILEDB_INT8:
            return py::cast(soma_dense_ndarr.non_empty_domain<int8_t>(name));
        case TILEDB_FLOAT64:
            return py::cast(soma_dense_ndarr.non_empty_domain<double>(name));
        case TILEDB_FLOAT32:
            return py::cast(soma_dense_ndarr.non_empty_domain<float>(name));
        case TILEDB_STRING_UTF8:
        case TILEDB_STRING_ASCII: 
            return py::cast(soma_dense_ndarr.non_empty_domain_var(name));
        default:
            throw TileDBSOMAError("Unsupported dtype for nonempty domain.");
        }
    })
    .def("domain", [](SOMADenseNDArray& soma_dense_ndarr, std::string name, py::dtype dtype) {
        switch (np_to_tdb_dtype(dtype)) {
        case TILEDB_UINT64:
            return py::cast(soma_dense_ndarr.domain<uint64_t>(name));
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
            return py::cast(soma_dense_ndarr.domain<int64_t>(name));
        case TILEDB_UINT32:
            return py::cast(soma_dense_ndarr.domain<uint32_t>(name));
        case TILEDB_INT32:
            return py::cast(soma_dense_ndarr.domain<int32_t>(name));
        case TILEDB_UINT16:
            return py::cast(soma_dense_ndarr.domain<uint16_t>(name));
        case TILEDB_INT16:
            return py::cast(soma_dense_ndarr.domain<int16_t>(name));
        case TILEDB_UINT8:
            return py::cast(soma_dense_ndarr.domain<uint8_t>(name));
        case TILEDB_INT8:
            return py::cast(soma_dense_ndarr.domain<int8_t>(name));
        case TILEDB_FLOAT64:
            return py::cast(soma_dense_ndarr.domain<double>(name));
        case TILEDB_FLOAT32:
            return py::cast(soma_dense_ndarr.domain<float>(name));
        case TILEDB_STRING_UTF8:
        case TILEDB_STRING_ASCII: {
            std::pair<std::string, std::string> str_domain;
            return py::cast(std::make_pair("", ""));
        }
        default:
            throw TileDBSOMAError("Unsupported dtype for Dimension's domain");
        }
    })
    .def("shape", &SOMADenseNDArray::shape)
    .def("ndim", &SOMADenseNDArray::ndim)
    .def("read_next", &SOMADenseNDArray::read_next)
    .def("write", &SOMADenseNDArray::write)
    .def("set_metadata", &SOMADenseNDArray::set_metadata)
    .def("delete_metadata", &SOMADenseNDArray::delete_metadata)
    .def_property_readonly("meta", [](SOMADenseNDArray& soma_dense_ndarr) -> py::dict {
        py::dict results;
            
        for (auto const& [key, val] : soma_dense_ndarr.get_metadata()){
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
        py::overload_cast<const std::string&>(&SOMADenseNDArray::get_metadata))
    .def("get_metadata", py::overload_cast<>(&SOMADenseNDArray::get_metadata))
    .def("has_metadata", &SOMADenseNDArray::has_metadata)
    .def("metadata_num", &SOMADenseNDArray::metadata_num);
}
}