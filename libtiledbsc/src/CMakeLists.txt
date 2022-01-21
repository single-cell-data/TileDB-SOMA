message(STATUS "Starting TileDB-SC build.")

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


target_compile_definitions(TILEDB_SC_OBJECTS PRIVATE
  -DBUILD_COMMIT_HASH="${BUILD_COMMIT_HASH}"
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
    message(STATUS "The TileDB-SC library is compiled with sanitizer ${SANITIZER} enabled")
  endif()
  target_compile_options(TILEDB_SC_OBJECTS
          PRIVATE
          -g -fno-omit-frame-pointer -fno-optimize-sibling-calls -fsanitize=${SANITIZER}
          )
endif()

############################################################
# tiledbsc library target
############################################################

add_library(tiledbsc SHARED
)

target_link_libraries(tiledbsc
  PUBLIC
    HTSlib::HTSlib
    TileDB::tiledb_shared
    spdlog::spdlog
)

# Sanitizer linker flags
if (SANITIZER)
  target_link_libraries(tiledbsc
          INTERFACE
          -fsanitize=${SANITIZER}
          )
endif()

# List of API headers (to be installed)
set(TILEDB_SC_PUBLIC_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/cpp_api/tiledbsc.h
)

if (ENABLE_ARROW_EXPORT)
  list(APPEND TILEDB_SC_PUBLIC_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/c_api/arrow.h
  )
endif()

############################################################
# API symbol exports (and public headers for install)
############################################################

include(GenerateExportHeader)

# Generates the file 'tiledbsc_export.h' suitable for the current compiler.
generate_export_header(TILEDB_SC_OBJECTS
  BASE_NAME tiledbsc
)

# Set variables in the parent scope so the tests can reference it.
set(TILEDB_SC_EXPORT_HEADER_NAME "tiledbsc_export.h" PARENT_SCOPE)
set(TILEDB_SC_EXPORT_HEADER_DIR "${CMAKE_CURRENT_BINARY_DIR}" PARENT_SCOPE)
set(TILEDB_SC_EXPORT_HEADER "${CMAKE_CURRENT_BINARY_DIR}/tiledbsc_export.h")
set(TILEDB_SC_EXPORT_HEADER "${TILEDB_SC_EXPORT_HEADER}" PARENT_SCOPE)

# Set related compiler settings
target_compile_definitions(TILEDB_SC_OBJECTS PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)
target_include_directories(TILEDB_SC_OBJECTS PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# Add the generated header to the public headers list
list(APPEND TILEDB_SC_PUBLIC_HEADERS
  "${TILEDB_SC_EXPORT_HEADER}"
)

# Set the public headers, which are the ones that get installed.
set_target_properties(tiledbsc
  PROPERTIES
    PUBLIC_HEADER "${TILEDB_SC_PUBLIC_HEADERS}"
)

############################################################
# Installation
############################################################

include(GNUInstallDirs)

# Set rpath to be relative to the .so.
if (APPLE)
  set_target_properties(tiledbsc PROPERTIES INSTALL_RPATH "@loader_path/")
else()
  set_target_properties(tiledbsc PROPERTIES INSTALL_RPATH "$ORIGIN/")
endif()

set_property(
  TARGET tiledbsc-bin
  PROPERTY INSTALL_RPATH
    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
    "${CMAKE_INSTALL_PREFIX}/lib"
)

install(
  TARGETS tiledbsc tiledbsc-bin
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tiledbsc
)
