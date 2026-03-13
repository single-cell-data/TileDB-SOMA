/**
 * @file   soma_collection.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the SOMACollection bindings.
 */

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsoma/tiledbsoma>

#include "common.h"

namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

void load_soma_collection(py::module& m) {
    py::class_<SOMACollectionBase, SOMAGroup, SOMAObject, py::smart_holder>(m, "SOMACollectionBase")
        // .def(
        //     "__iter__",
        //     [](SOMACollectionBase& collection) { return py::make_iterator(collection.begin(), collection.end()); },
        //     py::keep_alive<0, 1>())
        .def(
            "add",
            [](std::shared_ptr<SOMACollectionBase> collection,
               const std::string& uri,
               URIType uri_type,
               const std::string& name,
               const std::string& soma_type) { collection->set(uri, uri_type, name, soma_type); },
            "uri"_a,
            "uri_type"_a,
            "name"_a,
            "soma_type"_a)
        .def(
            "add",
            [](std::shared_ptr<SOMACollectionBase> collection,
               const std::string& uri,
               URIType uri_type,
               const std::string& name,
               const std::string& soma_type,
               std::shared_ptr<SOMAObject> member,
               bool managed) { collection->set(uri, uri_type, name, soma_type, member, managed); },
            "uri"_a,
            "uri_type"_a,
            "name"_a,
            "soma_type"_a,
            "member"_a,
            "managed"_a)
        .def("get", &SOMACollectionBase::get);

    py::class_<SOMACollection, SOMACollectionBase, py::smart_holder>(m, "SOMACollection")
        .def_static(
            "open",
            py::overload_cast<
                std::string_view,
                OpenMode,
                std::shared_ptr<SOMAContext>,
                std::optional<std::pair<uint64_t, uint64_t>>>(&SOMACollection::open),
            "uri"_a,
            py::kw_only(),
            "mode"_a,
            "context"_a,
            "timestamp"_a = py::none(),
            py::call_guard<py::gil_scoped_release>())
        .def_static(
            "create",
            [](std::shared_ptr<SOMAContext> ctx, std::string_view uri, std::optional<TimestampRange> timestamp) {
                try {
                    SOMACollection::create(uri, ctx, timestamp);
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
            },
            py::kw_only(),
            "ctx"_a,
            "uri"_a,
            "timestamp"_a = py::none());

    py::class_<SOMAExperiment, SOMACollectionBase, py::smart_holder>(m, "SOMAExperiment")
        .def_static(
            "open",
            &SOMAExperiment::open,
            "uri"_a,
            py::kw_only(),
            "mode"_a,
            "context"_a,
            "timestamp"_a = py::none(),
            py::call_guard<py::gil_scoped_release>())
        .def_static(
            "create",
            [](std::shared_ptr<SOMAContext> ctx, std::string_view uri, std::optional<TimestampRange> timestamp) {
                try {
                    SOMAExperiment::create(uri, ctx, timestamp);
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
            },
            py::kw_only(),
            "ctx"_a,
            "uri"_a,
            "timestamp"_a = py::none())
        .def_property_readonly("obs", &SOMAExperiment::obs)
        .def_property_readonly("ms", &SOMAExperiment::ms);

    py::class_<SOMAMeasurement, SOMACollectionBase, py::smart_holder>(m, "SOMAMeasurement")
        .def_static(
            "open",
            &SOMAMeasurement::open,
            "uri"_a,
            py::kw_only(),
            "mode"_a,
            "context"_a,
            "timestamp"_a = py::none(),
            py::call_guard<py::gil_scoped_release>())
        .def_static(
            "create",
            [](std::shared_ptr<SOMAContext> ctx, std::string_view uri, std::optional<TimestampRange> timestamp) {
                try {
                    SOMAMeasurement::create(uri, ctx, timestamp);
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
            },
            py::kw_only(),
            "ctx"_a,
            "uri"_a,
            "timestamp"_a = py::none())
        .def_property_readonly("var", &SOMAMeasurement::var)
        .def_property_readonly("X", &SOMAMeasurement::X)
        .def_property_readonly("obsm", &SOMAMeasurement::obsm)
        .def_property_readonly("obsp", &SOMAMeasurement::obsp)
        .def_property_readonly("varm", &SOMAMeasurement::varm)
        .def_property_readonly("varp", &SOMAMeasurement::varp);

    py::class_<SOMAScene, SOMACollectionBase, py::smart_holder>(m, "SOMAScene")
        .def_static(
            "create",
            [](std::shared_ptr<SOMAContext> ctx,
               std::string_view uri,
               const std::optional<std::vector<std::string>>& axis_names,
               const std::optional<std::vector<std::optional<std::string>>>& axis_units,
               std::optional<TimestampRange> timestamp) {
                if (axis_units.has_value() && !axis_names.has_value()) {
                    throw TileDBSOMAError("Cannot provide axis units without axis names.");
                }
                std::optional<SOMACoordinateSpace> coord_space{std::nullopt};
                if (axis_names.has_value()) {
                    if (axis_units.has_value()) {
                        coord_space = SOMACoordinateSpace(axis_names.value(), axis_units.value());
                    } else {
                        coord_space = SOMACoordinateSpace(axis_names.value());
                    }
                }
                try {
                    SOMAScene::create(uri, ctx, coord_space, timestamp);
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
            },
            py::kw_only(),
            "ctx"_a,
            "uri"_a,
            "axis_names"_a,
            "axis_units"_a,
            "timestamp"_a = py::none())
        .def_static(
            "open",
            &SOMAScene::open,
            "uri"_a,
            py::kw_only(),
            "mode"_a,
            "context"_a,
            "timestamp"_a = py::none(),
            py::call_guard<py::gil_scoped_release>())
        .def_property_readonly("img", &SOMAScene::img)
        .def_property_readonly("obsl", &SOMAScene::obsl)
        .def_property_readonly("varl", &SOMAScene::varl);

    py::class_<SOMAMultiscaleImage, SOMACollectionBase, py::smart_holder>(m, "SOMAMultiscaleImage")
        .def_static(
            "create",
            [](std::shared_ptr<SOMAContext> ctx,
               std::string_view uri,
               const std::vector<std::string>& axis_names,
               const std::vector<std::optional<std::string>>& axis_units,
               std::optional<TimestampRange> timestamp) {
                SOMACoordinateSpace coord_space{axis_names, axis_units};
                try {
                    SOMAMultiscaleImage::create(uri, ctx, coord_space, timestamp);
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
            },
            py::kw_only(),
            "ctx"_a,
            "uri"_a,
            "axis_names"_a,
            "axis_units"_a,
            "timestamp"_a = py::none())
        .def_static(
            "open",
            &SOMAMultiscaleImage::open,
            "uri"_a,
            py::kw_only(),
            "mode"_a,
            "context"_a,
            "timestamp"_a = py::none(),
            py::call_guard<py::gil_scoped_release>());
}
}  // namespace libtiledbsomacpp
