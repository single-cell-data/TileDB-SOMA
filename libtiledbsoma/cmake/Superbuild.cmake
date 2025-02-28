#
# Superbuild.cmake
#
# Licensed under the MIT License.
# Copyright (c) TileDB, Inc.
#

include(ExternalProject)

############################################################
# Common variables
############################################################

# Build paths for external projects
set(EP_BASE "${CMAKE_CURRENT_BINARY_DIR}/externals")
set(EP_SOURCE_DIR "${EP_BASE}/src")
set(EP_INSTALL_PREFIX "${EP_BASE}/install")

# A variable that will hold extra variables to pass to the regular
# non-superbuild build process as CMake arguments.
set(FORWARD_EP_CMAKE_ARGS)

# Variable that will hold a list of all the external projects added
# as a part of the superbuild.
set(EXTERNAL_PROJECTS)

# Forward any additional CMake args to the non-superbuild.
set(INHERITED_CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
  -DCMAKE_BUILD_TYPE=$<CONFIG>
  -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
  -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
  -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DEP_BASE=${EP_BASE}
  -DFORCE_BUILD_TILEDB=${FORCE_BUILD_TILEDB}
  -DTILEDB_S3=${TILEDB_S3}
  -DTILEDB_AZURE=${TILEDB_AZURe}
  -DTILEDB_GCS=${TILEDB_GCS}
  -DTILEDB_HDFS=${TILEDB_HDFS}
  -DTILEDB_SERIALIZATION=${TILEDB_SERIALIZATION}
  -DTILEDB_WERROR=${TILEDB_WERROR}
  -DTILEDB_VERBOSE=${TILEDB_VERBOSE}
  -DTILEDB_SANITIZER=${TILEDB_SANITIZER}
  -DTileDB_DIR=${TileDB_DIR}
  -DENABLE_ARROW_EXPORT=${ENABLE_ARROW_EXPORT}
  -DOVERRIDE_INSTALL_PREFIX=${OVERRIDE_INSTALL_PREFIX}
  -DTILEDBSOMA_BUILD_STATIC=${TILEDBSOMA_BUILD_STATIC}
  -DTILEDBSOMA_BUILD_CLI=${TILEDBSOMA_BUILD_CLI}
  -DTILEDBSOMA_ENABLE_TESTING=${TILEDBSOMA_ENABLE_TESTING}
  -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
  -DSPDLOG_LINK_SHARED=${SPDLOG_LINK_SHARED}
)

############################################################
# Set up external projects for dependencies
############################################################

# These includes modify the EXTERNAL_PROJECTS variable.

#TBD include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindCLI11_EP.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindTileDB_EP.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindSpdlog_EP.cmake)

############################################################
# 'make format' target
############################################################

set(SCRIPTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../scripts")

find_package(ClangTools)
if (${CLANG_FORMAT_FOUND})
  # Runs clang-format and updates files in place.
  add_custom_target(format ${SCRIPTS_DIR}/run-clang-format.sh ${CMAKE_CURRENT_SOURCE_DIR}/src ${CLANG_FORMAT_BIN} 1
    `find ${CMAKE_CURRENT_SOURCE_DIR} -name \\*.cc -or -name \\*.c -or -name \\*.h`)

  # Runs clang-format and exits with a non-zero exit code# if any files need to
  # be reformatted
  add_custom_target(check-format ${SCRIPTS_DIR}/run-clang-format.sh ${CMAKE_CURRENT_SOURCE_DIR}/src ${CLANG_FORMAT_BIN} 0
    `find ${CMAKE_CURRENT_SOURCE_DIR} -name \\*.cc -or -name \\*.c -or -name \\*.h`)
endif()

############################################################
# Set up the regular build (i.e. non-superbuild).
############################################################

ExternalProject_Add(libtiledbsoma
  SOURCE_DIR ${PROJECT_SOURCE_DIR}
  CMAKE_ARGS
    -DSUPERBUILD=OFF
    ${INHERITED_CMAKE_ARGS}
    ${FORWARD_EP_CMAKE_ARGS}
  INSTALL_COMMAND ""
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libtiledbsoma
  DEPENDS ${EXTERNAL_PROJECTS}
)

# make install-libtiledbsoma
add_custom_target(install-libtiledbsoma
  COMMAND ${CMAKE_COMMAND} --build . --target install --config $<CONFIG>
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/libtiledbsoma
)

# make check
add_custom_target(check
  COMMAND ${CMAKE_COMMAND} --build . --target check --config $<CONFIG>
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/libtiledbsoma
)


# TODO: add command to build add links to apis/python/src/tiledbsoma/
# add_custom_target(link_target ALL
#                  COMMAND ${CMAKE_COMMAND} -E create_symlink ${target} ${link})
