message(STATUS "Starting TileDB-SOMA build.")

set(TILEDBSOMA_INSTALL_TARGETS "")

# ###########################################################
# Find required dependencies
# See ../cmake/Modules/*.cmake for provenance/version info
# ###########################################################
find_package(TileDB_EP REQUIRED)
find_package(Spdlog_EP REQUIRED)

# ###########################################################
# Get source commit hash
# ###########################################################
find_package(Git REQUIRED)

execute_process(
  COMMAND "${GIT_EXECUTABLE}" describe --exact-match --tags HEAD
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  RESULT_VARIABLE res
  OUTPUT_VARIABLE BUILD_COMMIT_HASH
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# If we didn't find a tag name let's grab the SHA
if(res)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" describe --dirty=-modified --always
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RESULT_VARIABLE res
    OUTPUT_VARIABLE BUILD_COMMIT_HASH
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

set_property(GLOBAL APPEND
  PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${CMAKE_SOURCE_DIR}/.git/index")

message(STATUS "Building with commit hash ${BUILD_COMMIT_HASH}")

# ###########################################################
# Common object library
# ###########################################################
add_library(TILEDB_SOMA_OBJECTS OBJECT
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/managed_query.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_array.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_group.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_collection.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_experiment.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_measurement.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_dataframe.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_dense_ndarray.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_sparse_ndarray.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/array_buffers.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/column_buffer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/arrow_adapter.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/logger.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/stats.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/util.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/version.cc

  # TODO: uncomment when thread_pool is enabled
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/src/thread_pool/thread_pool.cc
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/src/thread_pool/status.cc
)

message(STATUS "Building TileDB without deprecation warnings")
target_compile_definitions(TILEDB_SOMA_OBJECTS PRIVATE
  -DBUILD_COMMIT_HASH="${BUILD_COMMIT_HASH}"
  -DTILEDB_DEPRECATED=
)

target_compile_options(
  TILEDB_SOMA_OBJECTS PRIVATE
)

set_property(TARGET TILEDB_SOMA_OBJECTS PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(TILEDB_SOMA_OBJECTS
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/external/include
  $<TARGET_PROPERTY:TileDB::tiledb_shared,INTERFACE_INCLUDE_DIRECTORIES>
  $<TARGET_PROPERTY:spdlog::spdlog,INTERFACE_INCLUDE_DIRECTORIES>
  ${pybind11_INCLUDE_DIRS}
)

# ###########################################################
# Compile options/definitions
# ###########################################################
if(SANITIZER)
  string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LOWER)

  if(NOT CMAKE_BUILD_TYPE_LOWER MATCHES "debug")
    message(FATAL_ERROR "Sanitizers only enabled for Debug build")
  endif()

  string(TOLOWER ${SANITIZER} SANITIZER)

  if(NOT SANITIZER MATCHES "^(address|memory|leak|thread|undefined)$")
    message(FATAL_ERROR "Unknown clang sanitizer: ${SANITIZER})")
  else()
    message(STATUS "The TileDB-SOMA library is compiled with sanitizer ${SANITIZER} enabled")
  endif()

  target_compile_options(TILEDB_SOMA_OBJECTS
    PRIVATE
    -g -fno-omit-frame-pointer -fno-optimize-sibling-calls -fsanitize=${SANITIZER}
  )
endif()

# ###########################################################
# tiledbsoma library target
# ###########################################################
if(TILEDBSOMA_BUILD_STATIC)
  add_library(tiledbsoma_static STATIC
    $<TARGET_OBJECTS:TILEDB_SOMA_OBJECTS>
  )
  list(APPEND TILEDBSOMA_INSTALL_TARGETS tiledbsoma_static)

  if(WIN32)
    # On Windows we must name the static library something else to avoid
    # name clash with the DLL's "import library" .lib file.
    set_target_properties(tiledbsoma_static
      PROPERTIES
      OUTPUT_NAME "tiledbsomastatic"
    )
  else()
    set_target_properties(tiledbsoma_static
      PROPERTIES
      OUTPUT_NAME "tiledbsoma"
    )
  endif()
else()
  add_library(tiledbsoma SHARED
    $<TARGET_OBJECTS:TILEDB_SOMA_OBJECTS>
  )

  list(APPEND TILEDBSOMA_INSTALL_TARGETS tiledbsoma)

  target_link_libraries(tiledbsoma
    PUBLIC
    TileDB::tiledb_shared
    spdlog::spdlog
  )
endif()

# Sanitizer linker flags
if(SANITIZER)
  if(TILEDBSOMA_BUILD_STATIC)
    target_link_libraries(tiledbsoma_static
      INTERFACE
      -fsanitize=${SANITIZER}
    )
  else()
    target_link_libraries(tiledbsoma
      INTERFACE
      -fsanitize=${SANITIZER}
    )
  endif()
endif()

# Install header files 
# target_sources FILE_SET is preferred with cmake>=3.23
# TODO Uncomment after finishing Python and R bindings
# set(TILEDB_SOMA_PUBLIC_HEADERS
#   ${CMAKE_CURRENT_SOURCE_DIR}/cpp_api/tiledbsoma
#   ${CMAKE_CURRENT_SOURCE_DIR}/cpp_api/soma_collection.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/cpp_api/soma_dataframe.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/cpp_api/soma_dense_ndarray.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/cpp_api/soma_experiment.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/cpp_api/soma_measurement.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/cpp_api/soma_object.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/cpp_api/soma_sparse_ndarray.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/cpp_api/logger_public.h
# )

install(FILES 
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/enums.h
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/logger_public.h
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/managed_query.h
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/array_buffers.h  
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/column_buffer.h  
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_array.h 
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_group.h
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_collection.h 
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_dataframe.h 
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_dense_ndarray.h
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_sparse_ndarray.h 
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_experiment.h 
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_measurement.h 
  ${CMAKE_CURRENT_SOURCE_DIR}/soma/soma_object.h 
  DESTINATION "include/tiledbsoma/soma"
)

install(FILES 
  ${CMAKE_CURRENT_SOURCE_DIR}/tiledbsoma/tiledbsoma
  DESTINATION "include/tiledbsoma"
)

install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/arrow_adapter.h  
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/carrow.h  
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/common.h  
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/stats.h  
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/util.h  
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/version.h
  DESTINATION "include/tiledbsoma/utils"
)

install(FILES 
  ${CMAKE_CURRENT_SOURCE_DIR}/external/include/span/span.hpp
  DESTINATION "include/tiledbsoma/soma/span"
)


# ###########################################################
# API symbol exports
# ###########################################################
include(GenerateExportHeader)

# Generates the file 'tiledbsoma_export.h' suitable for the current compiler.
generate_export_header(TILEDB_SOMA_OBJECTS
  BASE_NAME tiledbsoma
)

# Set variables in the parent scope so the tests can reference it.
set(TILEDB_SOMA_EXPORT_HEADER_NAME "tiledbsoma_export.h" PARENT_SCOPE)
set(TILEDB_SOMA_EXPORT_HEADER_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(TILEDB_SOMA_EXPORT_HEADER_DIR "${CMAKE_CURRENT_BINARY_DIR}" PARENT_SCOPE)
set(TILEDB_SOMA_EXPORT_HEADER "${CMAKE_CURRENT_BINARY_DIR}/tiledbsoma_export.h")
set(TILEDB_SOMA_EXPORT_HEADER "${TILEDB_SOMA_EXPORT_HEADER}" PARENT_SCOPE)

# Set related compiler settings
target_compile_definitions(TILEDB_SOMA_OBJECTS PRIVATE -DTILEDB_SOMA_OBJECTS_EXPORTS)
target_include_directories(TILEDB_SOMA_OBJECTS PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# Add the generated header to the public headers list
list(APPEND TILEDB_SOMA_PUBLIC_HEADERS
  "${TILEDB_SOMA_EXPORT_HEADER}"
)

# Add public headers, which are the ones that get installed.
if(NOT TILEDBSOMA_BUILD_STATIC)
  set_target_properties(tiledbsoma
    PROPERTIES
    PUBLIC_HEADER "${TILEDB_SOMA_PUBLIC_HEADERS}"
  )
endif()

# ###########################################################
# CLI executable target
# ###########################################################
if(TILEDBSOMA_BUILD_CLI)
  add_executable(tiledbsoma-cli
    ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli.cc
    $<TARGET_OBJECTS:TILEDB_SOMA_OBJECTS>
  )

  list(APPEND TILEDBSOMA_INSTALL_TARGETS tiledbsoma-cli)

  set_target_properties(tiledbsoma-cli PROPERTIES OUTPUT_NAME tdbsoma)

  target_link_libraries(tiledbsoma-cli
    PUBLIC

    # CLI11::CLI11
    # spdlog::spdlog
    tiledbsoma
    TileDB::tiledb_shared
  )

  # Sanitizer linker flags
  if(SANITIZER)
    target_link_libraries(tiledbsoma-cli
      INTERFACE
      -fsanitize=${SANITIZER}
    )
  endif()

  if(NOT APPLE AND NOT WIN32)
    target_link_libraries(tiledbsoma-cli PRIVATE pthread)
  endif()

  target_include_directories(tiledbsoma-cli
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/external/include
    ${TILEDB_SOMA_EXPORT_HEADER_DIR}
    ${pybind11_INCLUDE_DIRS}
  )
endif()

# ###########################################################
# Installation
# ###########################################################
# Get library directory for multiarch linux distros
include(GNUInstallDirs)

# Include module with function 'configure_package_config_file'
include(CMakePackageConfigHelpers)

# Set directory where TileDBConfig.cmake will be installed
set(CONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/tiledbsoma")

# Set name for export target file (will be installed to CONFIG_INSTALL_DIR)
set(TARGETS_EXPORT_NAME "TileDBSomaTargets")

# Path to generated cmake file
set(PROJECT_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/TileDBSomaConfig.cmake")

# Generate 'TileDBConfig.cmake'
# This process requires these variables to be defined at this point:
# * TARGETS_EXPORT_NAME
# * PROJECT_NAME
# * TILEDB_STATIC_DEP_STRING
configure_package_config_file(
  "${TILEDBSOMA_CMAKE_INPUTS_DIR}/Config.cmake.in"
  "${PROJECT_CONFIG}"
  INSTALL_DESTINATION "${CONFIG_INSTALL_DIR}"
)

# Install config file to <prefix>/lib/cmake/tiledbsoma/TileDBSomaConfig.cmake
install(
  FILES "${PROJECT_CONFIG}"
  DESTINATION "${CONFIG_INSTALL_DIR}"
)

# Set rpath so the TileDB-SOMA dynamic dependencies can be located.
if(NOT WIN32 AND NOT TILEDBSOMA_BUILD_STATIC)
  set_target_properties(tiledbsoma
    PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
  )
endif()

if(TILEDBSOMA_BUILD_CLI)
  set_target_properties(tiledbsoma-cli PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
endif()

install(
  TARGETS ${TILEDBSOMA_INSTALL_TARGETS}
  EXPORT ${TARGETS_EXPORT_NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tiledbsoma
)

# Install targets file to <prefix>/lib/cmake/tiledbsoma/TileDBSOMATargets.cmake
install(
  EXPORT "${TARGETS_EXPORT_NAME}"
  NAMESPACE "${PROJECT_NAME}::"
  DESTINATION "${CONFIG_INSTALL_DIR}"
)

# PKG Config file
configure_file(
  "${TILEDBSOMA_CMAKE_INPUTS_DIR}/tiledbsoma.pc.in"
  ${CMAKE_CURRENT_BINARY_DIR}/tiledbsoma.pc
  @ONLY
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/tiledbsoma.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
