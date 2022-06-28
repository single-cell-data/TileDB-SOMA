#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsc/tiledbsc>

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
                        py::array data) {
                return ColumnBuffer::create(
                    name,
                    type,
                    num_cells,
                    std::span<std::byte>(
                        (std::byte*)data.data(), data.nbytes()));
            }),
            py::arg("name"),
            py::arg("type"),
            py::arg("num_cells"),
            py::arg("data"))

        // WARNING: these functions copy!
        .def("data", &ColumnBuffer::py_array);

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
