#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsc/tiledbsc>

#include "tiledbsc/tiledb_arrow.h"

#define DENUM(x) .value(#x, TILEDB_##x)

using namespace tiledbsc;

namespace py = pybind11;

PYBIND11_MODULE(pytiledbsc, m) {
    m.doc() = "TileDB-SingleCell python library";

    // TODO: ColumnBuffer is useful for testing, but may be removed later
    py::class_<ColumnBuffer, std::shared_ptr<ColumnBuffer>>(m, "ColumnBuffer")
        .def(
            py::init([](std::string& name,
                        tiledb_datatype_t type,
                        py::int_ num_cells,
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

        .def(
            "to_arrow",
            [](ColumnBuffer& buf) {
                auto [array, schema] = ArrowAdapter::to_arrow(buf);
                return py::make_tuple(
                    py::capsule(array.get()), py::capsule(schema.get()));
            })

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

    py::class_<SOMA, std::shared_ptr<SOMA>>(m, "SOMA")
        .def(py::init([](std::string_view uri) { return SOMA::open(uri); }))
        .def("list_arrays", &SOMA::list_arrays);

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
