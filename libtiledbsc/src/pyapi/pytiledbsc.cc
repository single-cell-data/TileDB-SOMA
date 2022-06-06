#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include <tiledbsc/managed_query.h>
#include <tiledbsc/tiledbsc.h>

namespace {}

namespace tiledbsc {

using namespace std;

namespace py = pybind11;

PYBIND11_MODULE(pytiledbsc, m) {
    py::class_<QueryResult>(m, "QueryResult");
    //.def(py::init<>());
}

};  // namespace tiledbsc