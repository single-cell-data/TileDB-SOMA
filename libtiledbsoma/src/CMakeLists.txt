message(STATUS "Starting TileDB-SOMA build.")

set(TILEDBSOMA_INSTALL_TARGETS "")

############################################################
# Find required dependencies
############################################################

find_package(TileDB_EP REQUIRED)
find_package(Spdlog_EP REQUIRED)

############################################################
# Get source commit hash
############################################################

find_package(Git REQUIRED)

execute_process(
  COMMAND "${GIT_EXECUTABLE}" describe --exact-match --tags HEAD
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  RESULT_VARIABLE res
  OUTPUT_VARIABLE BUILD_COMMIT_HASH
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# If we didn't find a tag name let's grab the SHA
if (res)
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

############################################################
# Common object library
############################################################

add_library(TILEDB_SOMA_OBJECTS OBJECT
  ${CMAKE_CURRENT_SOURCE_DIR}/column_buffer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/logger.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/managed_query.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/soma_reader.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/util.cc
# TODO: uncomment when thread_pool is enabled
#  ${CMAKE_CURRENT_SOURCE_DIR}/../external/src/thread_pool/thread_pool.cc
#  ${CMAKE_CURRENT_SOURCE_DIR}/../external/src/thread_pool/status.cc
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
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/include
    $<TARGET_PROPERTY:TileDB::tiledb_shared,INTERFACE_INCLUDE_DIRECTORIES>
    $<TARGET_PROPERTY:spdlog::spdlog,INTERFACE_INCLUDE_DIRECTORIES>
    ${pybind11_INCLUDE_DIRS}
)

############################################################
# Compile options/definitions
############################################################

if (SANITIZER)
  string( TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LOWER )
  if (NOT CMAKE_BUILD_TYPE_LOWER MATCHES "debug")
    message(FATAL_ERROR "Sanitizers only enabled for Debug build")
  endif()
  string(TOLOWER ${SANITIZER} SANITIZER)
  if (NOT SANITIZER MATCHES "^(address|memory|leak|thread|undefined)$")
    message(FATAL_ERROR "Unknown clang sanitizer: ${SANITIZER})")
  else()
    message(STATUS "The TileDB-SOMA library is compiled with sanitizer ${SANITIZER} enabled")
  endif()
  target_compile_options(TILEDB_SOMA_OBJECTS
          PRIVATE
          -g -fno-omit-frame-pointer -fno-optimize-sibling-calls -fsanitize=${SANITIZER}
          )
endif()

############################################################
# tiledbsoma library target
############################################################

add_library(tiledbsoma STATIC
  $<TARGET_OBJECTS:TILEDB_SOMA_OBJECTS>
)

list(APPEND TILEDBSOMA_INSTALL_TARGETS tiledbsoma)

target_link_libraries(tiledbsoma
  PUBLIC
    TileDB::tiledb_shared
    spdlog::spdlog
)

# Sanitizer linker flags
if (SANITIZER)
  target_link_libraries(tiledbsoma
          INTERFACE
          -fsanitize=${SANITIZER}
          )
endif()

# Install header files and preserve directory structure
# target_sources FILE_SET is preferred with cmake>=3.23
install(DIRECTORY 
  ${CMAKE_CURRENT_SOURCE_DIR}/../external/include/span
  ${CMAKE_CURRENT_SOURCE_DIR}/../external/include/thread_pool
  ${CMAKE_CURRENT_SOURCE_DIR}/../include/tiledbsoma 
  TYPE INCLUDE
)

############################################################
# API symbol exports
############################################################

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

# Add public headers, which are the ones that get installed.
set_target_properties(tiledbsoma
  PROPERTIES
    PUBLIC_HEADER "${TILEDB_SOMA_EXPORT_HEADER}"
)

############################################################
# Python API bindings
############################################################

if (TILEDBSOMA_BUILD_PYTHON)
  add_subdirectory(pyapi)
endif()

############################################################
# CLI executable target
############################################################

if (TILEDBSOMA_BUILD_CLI)
  add_executable(tiledbsoma-cli
    ${CMAKE_CURRENT_SOURCE_DIR}/cli.cc
    $<TARGET_OBJECTS:TILEDB_SOMA_OBJECTS>
  )

  list(APPEND TILEDBSOMA_INSTALL_TARGETS tiledbsoma-cli)

  set_target_properties(tiledbsoma-cli PROPERTIES OUTPUT_NAME tdbsoma)

  target_link_libraries(tiledbsoma-cli
    PUBLIC
      #CLI11::CLI11
      #spdlog::spdlog
      tiledbsoma
      TileDB::tiledb_shared
  )

  # Sanitizer linker flags
  if (SANITIZER)
    target_link_libraries(tiledbsoma-cli
            INTERFACE
            -fsanitize=${SANITIZER}
            )
  endif()

  if (NOT APPLE AND NOT WIN32)
    target_link_libraries(tiledbsoma-cli PRIVATE pthread)
  endif()

  target_include_directories(tiledbsoma-cli
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/../include
      ${CMAKE_CURRENT_SOURCE_DIR}/../external/include
      ${TILEDB_SOMA_EXPORT_HEADER_DIR}
      ${pybind11_INCLUDE_DIRS}
  )

  set_target_properties(tiledbsoma-cli PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

endif()

############################################################
# Installation
############################################################

include(GNUInstallDirs)

install(
  TARGETS ${TILEDBSOMA_INSTALL_TARGETS}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tiledbsoma
)
