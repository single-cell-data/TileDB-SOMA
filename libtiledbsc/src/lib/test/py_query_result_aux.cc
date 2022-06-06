#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include <tiledb/tiledb.h>
#include <tiledb/tiledb>

#include <tiledbsc/managed_query.h>
#include <tiledbsc/query_result.h>
#include <tiledbsc/sc_arrowio.h>
#include <tiledbsc/tiledbsc.h>

#define DENUM(x) .value(#x, TILEDB_##x)

namespace tiledbsc_aux {

using namespace std;
using namespace tiledbsc;

namespace py = pybind11;

PYBIND11_MODULE(py_query_result_aux, m) {
    py::class_<tiledbsc::arrow::ArrowPair>(m, "ArrowPair")
        .def(
            "schema",
            [](tiledbsc::arrow::ArrowPair& p) -> py::int_ {
                return (ptrdiff_t)p.schema;
            })
        .def("array", [](tiledbsc::arrow::ArrowPair& p) -> py::int_ {
            return (ptrdiff_t)p.array;
        });

    py::class_<tiledbsc::QueryResult>(m, "QueryResult")
        .def(py::init<std::map<string, BufferSet>>())
        .def(
            "to_arrow", &QueryResult::to_arrow, py::arg("name") = std::nullopt);

    py::enum_<tiledb_datatype_t>(m, "DataType", py::module_local()) DENUM(
        INT32) DENUM(INT64) DENUM(FLOAT32) DENUM(FLOAT64) DENUM(CHAR)
        DENUM(INT8) DENUM(UINT8) DENUM(INT16) DENUM(UINT16) DENUM(UINT32) DENUM(
            UINT64) DENUM(STRING_ASCII) DENUM(STRING_UTF8) DENUM(STRING_UTF16)
            DENUM(STRING_UTF32) DENUM(STRING_UCS2) DENUM(STRING_UCS4) DENUM(ANY)
                DENUM(DATETIME_YEAR) DENUM(DATETIME_WEEK) DENUM(DATETIME_DAY)
                    DENUM(DATETIME_HR) DENUM(DATETIME_MIN) DENUM(DATETIME_SEC)
                        DENUM(DATETIME_MS) DENUM(DATETIME_US) DENUM(DATETIME_NS)
                            DENUM(DATETIME_PS) DENUM(DATETIME_FS)
                                DENUM(DATETIME_AS) DENUM(TIME_HR)
                                    DENUM(TIME_MIN) DENUM(TIME_SEC)
                                        DENUM(TIME_MS) DENUM(TIME_US)
                                            DENUM(TIME_NS) DENUM(TIME_PS)
                                                DENUM(TIME_FS) DENUM(TIME_AS);

    py::class_<tiledbsc::BufferSet>(m, "BufferSet")
        .def(py::init(
            [](std::string& name,
               tiledb_datatype_t data_type,
               py::int_ nbytes,
               py::array data,
               std::optional<py::array_t<uint64_t>> offsets = std::nullopt,
               std::optional<py::array_t<byte>> validity = std::nullopt) {
                auto offsets_v = offsets ?
                                     std::make_optional<std::vector<uint64_t>>(
                                         offsets.value().data(),
                                         offsets.value().data() +
                                             offsets.value().size()) :
                                     std::nullopt;
                auto validity_v = validity ?
                                      std::make_optional<std::vector<byte>>(
                                          validity.value().data(),
                                          validity.value().data() +
                                              validity.value().size()) :
                                      std::nullopt;

                tiledbsc::BufferSet bs{
                    name,
                    data_type,
                    nbytes,
                    std::vector<byte>(
                        (byte*)data.data(), (byte*)data.data() + data.nbytes()),
                    offsets_v,
                    validity_v,
                };
                return std::make_unique<tiledbsc::BufferSet>(bs);
            }))
        /*
         * WARNING: these functions copy!
         */
        .def(
            "data",
            [](BufferSet& bfs) -> py::array_t<byte> {
                return py::array_t<byte>(
                    bfs.data<byte>().size(), bfs.data<byte>().data());
            })
        .def(
            "offsets",
            [](BufferSet& bfs) -> std::optional<py::array_t<uint64_t>> {
                if (!bfs.isvar())
                    return std::nullopt;

                auto o = bfs.offsets();
                return py::array_t<uint64_t>(o.size(), o.data());
            })
        .def(
            "validity", [](BufferSet& bfs) -> std::optional<py::array_t<byte>> {
                if (!bfs.isnullable())
                    return std::nullopt;

                auto v = bfs.validity();
                return py::array_t<byte>(v.size(), v.data());
            });
}

};  // namespace tiledbsc_aux