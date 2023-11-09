#include <tiledbsoma/tiledbsoma>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

using namespace tiledbsoma;

namespace tiledbsoma {

namespace py = pybind11;
using namespace py::literals;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void load_soma_array(py::module &);
void load_soma_dataframe(py::module &);
void load_query_condition(py::module &);

PYBIND11_MODULE(pytiledbsoma, m) {
  py::enum_<OpenMode>(m, "OpenMode")
      .value("read", OpenMode::read)
      .value("write", OpenMode::write);

  py::enum_<ResultOrder>(m, "ResultOrder")
      .value("automatic", ResultOrder::automatic)
      .value("rowmajor", ResultOrder::rowmajor)
      .value("colmajor", ResultOrder::colmajor);

  m.doc() = "SOMA acceleration library";

  m.def("version", []() { return tiledbsoma::version::as_string(); });

  m.def(
      "config_logging",
      [](const std::string& level, const std::string& logfile) {
          LOG_CONFIG(level, logfile);
      },
      "level"_a,
      "logfile"_a = "");

  m.def("info", &LOG_INFO, "message"_a = "");
  m.def("debug", &LOG_DEBUG, "message"_a = "");

  m.def(
      "tiledbsoma_stats_enable",
      []() { tiledbsoma::stats::enable(); },
      "Enable TileDB internal statistics. Lifecycle: experimental.");
  m.def(
      "tiledbsoma_stats_disable",
      []() { tiledbsoma::stats::disable(); },
      "Disable TileDB internal statistics. Lifecycle: experimental.");
  m.def(
      "tiledbsoma_stats_reset",
      []() { tiledbsoma::stats::reset(); },
      "Reset all TileDB internal statistics to 0. Lifecycle: experimental.");
  m.def(
      "tiledbsoma_stats_dump",
      []() {
          py::print(tiledbsoma::version::as_string());
          std::string stats = tiledbsoma::stats::dump();
          py::print(stats);
      },
      "Print TileDB internal statistics. Lifecycle: experimental.");

  load_soma_array(m);
  load_soma_dataframe(m);
  load_query_condition(m);
}

}; 