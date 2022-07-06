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
}
