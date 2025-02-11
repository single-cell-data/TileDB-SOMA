/**
 * @file   soma_object.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the SOMAObject bindings.
 */

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <memory>
#include <tiledbsoma/tiledbsoma>

#include "common.h"

namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

void load_soma_object(py::module& m) {
    py::class_<SOMAObject>(m, "SOMAObject")

        .def_static(
            "open",
            [](std::string_view uri,
               OpenMode mode,
               std::shared_ptr<SOMAContext> context,
               std::optional<std::pair<uint64_t, uint64_t>> timestamp,
               std::optional<std::string> clib_type) -> py::object {
                auto soma_obj = ([&]() {
                    py::gil_scoped_release release;
                    return SOMAObject::open(
                        uri, mode, context, timestamp, clib_type);
                })();

                std::optional<std::string> soma_obj_type;
                try {
                    soma_obj_type = soma_obj->type();
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }

                if (!soma_obj_type) {
                    assert(
                        false &&
                        "Unreachable code: The missing soma_object_type case "
                        "is already handled. This indicates an "
                        "unexpected failure to catch exceptions by "
                        "SOMAObject::open");
                }

                std::transform(
                    soma_obj_type->begin(),
                    soma_obj_type->end(),
                    soma_obj_type->begin(),
                    [](unsigned char c) { return std::tolower(c); });

                if (soma_obj_type == "somadataframe")
                    return py::cast(dynamic_cast<SOMADataFrame&>(*soma_obj));
                else if (soma_obj_type == "somapointclouddataframe")
                    return py::cast(
                        dynamic_cast<SOMAPointCloudDataFrame&>(*soma_obj));
                else if (soma_obj_type == "somageometrydataframe")
                    return py::cast(
                        dynamic_cast<SOMAGeometryDataFrame&>(*soma_obj));
                else if (soma_obj_type == "somasparsendarray")
                    return py::cast(
                        dynamic_cast<SOMASparseNDArray&>(*soma_obj));
                else if (soma_obj_type == "somadensendarray")
                    return py::cast(dynamic_cast<SOMADenseNDArray&>(*soma_obj));
                else if (soma_obj_type == "somacollection")
                    return py::cast(dynamic_cast<SOMACollection&>(*soma_obj));
                else if (soma_obj_type == "somaexperiment")
                    return py::cast(dynamic_cast<SOMAExperiment&>(*soma_obj));
                else if (soma_obj_type == "somameasurement")
                    return py::cast(dynamic_cast<SOMAMeasurement&>(*soma_obj));
                else if (soma_obj_type == "somascene")
                    return py::cast(dynamic_cast<SOMAScene&>(*soma_obj));
                else if (soma_obj_type == "somamultiscaleimage")
                    return py::cast(
                        dynamic_cast<SOMAMultiscaleImage&>(*soma_obj));

                assert(
                    false &&
                    "Unreachable code: All possible SOMA object types are "
                    "already handled. This indicates a logic error or an "
                    "unexpected failure to catch exceptions by "
                    "SOMAObject::open");

                return py::none();  // Unreached, but appeases a compiler
                                    // warning
            },
            "uri"_a,
            "mode"_a,
            "context"_a,
            py::kw_only(),
            "timestamp"_a = py::none(),
            "clib_type"_a = py::none())
        .def_property_readonly("type", &SOMAObject::type);
};
}  // namespace libtiledbsomacpp
