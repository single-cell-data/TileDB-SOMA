#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsc/tiledbsc>

#include "arrow_adapter.h"

#define DENUM(x) .value(#x, TILEDB_##x)

using namespace tiledbsc;

namespace py = pybind11;
using namespace py::literals;

/**
 * @brief Convert ColumnBuffer to Arrow array.
 *
 * @param cb ColumnBuffer
 * @return py::object Arrow array
 */
py::object to_array(ColumnBuffer& cb) {
    auto pa = py::module::import("pyarrow");
    auto pa_array_import = pa.attr("Array").attr("_import_from_c");

    auto [array, schema] = ArrowAdapter::to_arrow(cb);
    return pa_array_import(py::capsule(array.get()), py::capsule(schema.get()));
}

/**
 * @brief Convert ColumnBuffers to Arrow table.
 *
 * @param cbs ColumnBuffers
 * @return py::object
 */
py::object to_table(ColumnBuffers& cbs) {
    auto pa = py::module::import("pyarrow");
    auto pa_table_from_arrays = pa.attr("Table").attr("from_arrays");

    py::list names;
    py::list arrays;

    for (auto& [name, column] : cbs) {
        names.append(name);
        arrays.append(to_array(*column));
    }

    return pa_table_from_arrays(arrays, names);
}

// TODO: Convert Arrow array to ColumnBuffer
/*
    _export_to_c(...) method of pyarrow.lib.StringArray instance
    Array._export_to_c(self, out_ptr, out_schema_ptr=0)

    Export to a C ArrowArray struct, given its pointer.

    If a C ArrowSchema struct pointer is also given, the array type
    is exported to it at the same time.

    Parameters
    ----------
    out_ptr: int
        The raw pointer to a C ArrowArray struct.
    out_schema_ptr: int (optional)
        The raw pointer to a C ArrowSchema struct.

    Be careful: if you don't pass the ArrowArray struct to a consumer,
    array memory will leak.  This is a low-level function intended for
    expert users.
*/

/**
 * @brief pybind11 bindings
 *
 */
PYBIND11_MODULE(libtiledbsc, m) {
    m.doc() = "TileDB-SOMA acceleration library";

    py::class_<SOMA>(m, "SOMA")
        .def(
            py::init(
                [](std::string_view uri,
                   std::optional<std::map<std::string, std::string>> config) {
                    if (config.has_value()) {
                        auto cfg = Config(*config);
                        return SOMA::open(uri, cfg);
                    } else {
                        return SOMA::open(uri);
                    }
                }),
            "uri"_a,
            "config"_a = py::none())
        .def("list_arrays", &SOMA::list_arrays)
        // SOMAQuery (0 = return value) will keep SOMA alive (1 = this)
        .def("query", &SOMA::query, py::keep_alive<0, 1>());

    py::class_<SOMAQuery>(m, "SOMAQuery")
        .def("next_results", [](SOMAQuery& sq) -> std::optional<py::object> {
            auto buffers = sq.next_results();
            if (buffers.has_value()) {
                return to_table(*buffers);
            }
            return std::nullopt;
        });

    //===============================================================
    // Code below is provided for testing
    //===============================================================
    m.def(
        "to_arrow",
        [](std::map<std::string, std::shared_ptr<ColumnBuffer>> tb) {
            return to_table(tb);
        });

    py::class_<ColumnBuffer, std::shared_ptr<ColumnBuffer>>(m, "ColumnBuffer")
        .def(
            py::init([](std::string& name,
                        tiledb_datatype_t type,
                        int num_cells,
                        py::array data,
                        std::optional<py::array_t<uint64_t>> offsets,
                        std::optional<py::array_t<uint8_t>> validity) {
                auto o = offsets ?
                             std::span<uint64_t>(
                                 (uint64_t*)offsets->data(), offsets->size()) :
                             std::span<uint64_t>{};
                auto v = validity ?
                             std::span<uint8_t>(
                                 (uint8_t*)validity->data(), validity->size()) :
                             std::span<uint8_t>{};

                return ColumnBuffer::create(
                    name,
                    type,
                    num_cells,
                    std::span<std::byte>(
                        (std::byte*)data.data(), data.nbytes()),
                    o,
                    v);
            }),
            py::arg("name"),
            py::arg("type"),
            py::arg("num_cells"),
            py::arg("data"),
            py::arg("offsets") = std::nullopt,
            py::arg("validity") = std::nullopt)

        .def("to_arrow", [](ColumnBuffer& cb) { return to_array(cb); })

        // WARNING: these functions copy!
        .def(
            "data",
            [](ColumnBuffer& buf) -> py::array {
                switch (buf.type()) {
                    case TILEDB_INT32:
                        return py::array_t<int32_t>(
                            buf.data<int32_t>().size(),
                            buf.data<int32_t>().data());
                    case TILEDB_INT64:
                        return py::array_t<int64_t>(
                            buf.data<int64_t>().size(),
                            buf.data<int64_t>().data());
                    case TILEDB_FLOAT32:
                        return py::array_t<float>(
                            buf.data<float>().size(), buf.data<float>().data());
                    case TILEDB_FLOAT64:
                        return py::array_t<double>(
                            buf.data<double>().size(),
                            buf.data<double>().data());
                    case TILEDB_STRING_ASCII:
                        return py::array_t<char>(
                            buf.data<char>().size(), buf.data<char>().data());
                    default:
                        throw TileDBSCError(
                            "[ColumnBuffer] Unsupported type: " + buf.type());
                }
            })
        .def(
            "offsets",
            [](ColumnBuffer& buf) -> std::optional<py::array> {
                if (!buf.is_var()) {
                    return std::nullopt;
                }

                auto o = buf.offsets();
                return py::array_t<uint64_t>(o.size(), o.data());
            })
        .def("validity", [](ColumnBuffer& buf) -> std::optional<py::array> {
            if (!buf.is_nullable()) {
                return std::nullopt;
            }

            auto v = buf.validity();
            return py::array_t<uint8_t>(v.size(), v.data());
        });

    // clang-format off
    py::enum_<tiledb_datatype_t>(m, "DataType", py::module_local())
    DENUM(INT32)
    DENUM(INT64)
    DENUM(FLOAT32)
    DENUM(FLOAT64)
    DENUM(CHAR)
    DENUM(INT8)
    DENUM(UINT8)
    DENUM(INT16)
    DENUM(UINT16)
    DENUM(UINT32)
    DENUM(UINT64)
    DENUM(STRING_ASCII)
    DENUM(STRING_UTF8)
    DENUM(STRING_UTF16)
    DENUM(STRING_UTF32)
    DENUM(STRING_UCS2)
    DENUM(STRING_UCS4)
    DENUM(ANY)
    DENUM(DATETIME_YEAR)
    DENUM(DATETIME_WEEK)
    DENUM(DATETIME_DAY)
    DENUM(DATETIME_HR)
    DENUM(DATETIME_MIN)
    DENUM(DATETIME_SEC)
    DENUM(DATETIME_MS)
    DENUM(DATETIME_US)
    DENUM(DATETIME_NS)
    DENUM(DATETIME_PS)
    DENUM(DATETIME_FS)
    DENUM(DATETIME_AS)
    DENUM(TIME_HR)
    DENUM(TIME_MIN)
    DENUM(TIME_SEC)
    DENUM(TIME_MS)
    DENUM(TIME_US)
    DENUM(TIME_NS)
    DENUM(TIME_PS)
    DENUM(TIME_FS)
    DENUM(TIME_AS);
    // clang-format on
}
