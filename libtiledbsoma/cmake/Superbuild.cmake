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
# Filter out vcpkg-installed paths from CMAKE_PREFIX_PATH to avoid "extra path" warnings
# vcpkg will add these paths automatically when its toolchain file loads in the nested configure
# But we need to keep EP_INSTALL_PREFIX so TileDB can be found
if(CMAKE_PREFIX_PATH)
  set(CMAKE_PREFIX_PATH_FILTERED)
  foreach(path ${CMAKE_PREFIX_PATH})
    # Skip vcpkg-installed paths - vcpkg will add them automatically
    if(NOT path MATCHES "vcpkg_installed")
      list(APPEND CMAKE_PREFIX_PATH_FILTERED "${path}")
    endif()
  endforeach()
  set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH_FILTERED}")
endif()
# Always add EP_INSTALL_PREFIX so TileDB (built by superbuild) can be found
list(APPEND CMAKE_PREFIX_PATH "${EP_INSTALL_PREFIX}")
# Join the filtered list with semicolons for passing to CMake command
list(JOIN CMAKE_PREFIX_PATH ";" CMAKE_PREFIX_PATH_STRING)

# Forward CMAKE_TOOLCHAIN_FILE and CMAKE_OSX_ARCHITECTURES if set (needed for vcpkg in nested configure)
set(INHERITED_TOOLCHAIN_ARG)
if(CMAKE_TOOLCHAIN_FILE)
  list(APPEND INHERITED_TOOLCHAIN_ARG -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE})
endif()
if(CMAKE_OSX_ARCHITECTURES)
  list(APPEND INHERITED_TOOLCHAIN_ARG -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES})
endif()

set(INHERITED_CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH_STRING}
  -DCMAKE_BUILD_TYPE=$<CONFIG>
  ${INHERITED_TOOLCHAIN_ARG}
  -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
  -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
  -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DEP_INSTALL_PREFIX=${EP_INSTALL_PREFIX}
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

# Only include spdlog superbuild if vcpkg is not being used
# When using vcpkg, spdlog comes from vcpkg; superbuild only builds TileDB
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindSpdlog_EP.cmake)
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
