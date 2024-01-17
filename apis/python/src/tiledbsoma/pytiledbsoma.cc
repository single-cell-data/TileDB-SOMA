#include <tiledbsoma/tiledbsoma>
#include <tiledbsoma/reindexer/reindexer.h>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include "common.h"

namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void load_soma_array(py::module &);
void load_soma_object(py::module &);
void load_soma_dataframe(py::module &);
void load_query_condition(py::module &);

PYBIND11_MODULE(pytiledbsoma, m) {
    py::register_exception<TileDBSOMAError>(m, "SOMAError");

    /* We need to make sure C++ TileDBSOMAError is translated to a correctly-typed 
    * Python error
    */
    py::register_exception_translator([](std::exception_ptr p) {
    auto tiledb_soma_error =
        (py::object)py::module::import("tiledbsoma").attr("SOMAError");

    try {
        if (p)
        std::rethrow_exception(p);
    } catch (const TileDBSOMAError &e) {
        PyErr_SetString(tiledb_soma_error.ptr(), e.what());
    } catch (const TileDBSOMAPyError &e) {
        PyErr_SetString(tiledb_soma_error.ptr(), e.what());
    } catch (py::builtin_exception &e) {
        throw;
    };
    });

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

    // Efficient C++ re-indexing (aka hashing unique key values to an index
    // between 0 and number of keys - 1) based on khash
    py::class_<IntIndexer>(m, "IntIndexer")
        .def(py::init<>())
        .def(py::init<std::vector<int64_t>&, int>())
        .def(
            "map_locations",
            [](IntIndexer& indexer,
                py::array_t<int64_t> keys,
                int num_threads) {
                auto buffer = keys.request();
                int64_t* data = static_cast<int64_t*>(buffer.ptr);
                size_t length = buffer.shape[0];
                indexer.map_locations(keys.data(), keys.size(), num_threads);
            })
        .def(
            "map_locations",
            [](IntIndexer& indexer,
                std::vector<int64_t> keys,
                int num_threads) {
                indexer.map_locations(keys.data(), keys.size(), num_threads);
            })
        // Perform lookup for a large input array of keys and return the looked
        // up value array (passing ownership from C++ to python)
        .def(
            "get_indexer",
            [](IntIndexer& indexer, py::array_t<int64_t> lookups) {
                auto input_buffer = lookups.request();
                int64_t* input_ptr = static_cast<int64_t*>(input_buffer.ptr);
                size_t size = input_buffer.shape[0];
                auto results = py::array_t<int64_t>(size);
                auto results_buffer = results.request();
                size_t results_size = results_buffer.shape[0];

                int64_t* results_ptr = static_cast<int64_t*>(
                    results_buffer.ptr);

                indexer.lookup(input_ptr, results_ptr, size);
                return results;
            })
        // Perform lookup for a large input array of keys and writes the looked
        // up values into previously allocated array (works for the cases in
        // which python and R pre-allocate the array)
        .def(
            "get_indexer",
            [](IntIndexer& indexer,
                py::array_t<int64_t> lookups,
                py::array_t<int64_t>& results) {
                auto input_buffer = lookups.request();
                int64_t* input_ptr = static_cast<int64_t*>(input_buffer.ptr);
                size_t size = input_buffer.shape[0];

                auto results_buffer = results.request();
                int64_t* results_ptr = static_cast<int64_t*>(
                    results_buffer.ptr);
                size_t results_size = input_buffer.shape[0];
                indexer.lookup(input_ptr, input_ptr, size);
            });

    load_soma_array(m);
    load_soma_object(m);
    load_soma_dataframe(m);
    load_query_condition(m);
}

}; 
