############################################################
# Definitions
############################################################

# Set root source directory as compiler definition for use in tests
get_filename_component(TILEDBSOMA_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../" ABSOLUTE)
add_compile_definitions(TILEDBSOMA_SOURCE_ROOT="${TILEDBSOMA_SOURCE_ROOT}")

############################################################
# Dependencies
############################################################

find_package(TileDB_EP REQUIRED)
find_package(Spdlog_EP REQUIRED)

############################################################
# SOMA unit test
############################################################

find_package(Catch_EP REQUIRED)

add_executable(unit_soma EXCLUDE_FROM_ALL
    $<TARGET_OBJECTS:TILEDB_SOMA_OBJECTS>
    unit_column_buffer.cc
    unit_managed_query.cc
    unit_soma_reader.cc
# TODO: uncomment when thread_pool is enabled
#    unit_thread_pool.cc
)

target_link_libraries(unit_soma
  PRIVATE
    Catch2::Catch2WithMain
    TileDB::tiledb_shared
)

target_include_directories(unit_soma
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/include
    ${pybind11_INCLUDE_DIRS}
    $<TARGET_PROPERTY:spdlog::spdlog,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(unit_soma PRIVATE CATCH_CONFIG_MAIN)

# Allow deprecated function for writing to an array with a timestamp
target_compile_options(unit_soma PRIVATE -Wno-deprecated-declarations)

if (NOT APPLE AND NOT WIN32)
    target_link_libraries(unit_soma PRIVATE pthread)
endif()

add_test(
    NAME "unit_soma"
    COMMAND $<TARGET_FILE:unit_soma> "--durations=yes"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_target(build_tests)
add_dependencies(build_tests
  unit_soma
)

############################################################
# make check
############################################################

# Add custom target 'check' in the root of the tree
add_custom_target(
  check COMMAND ${CMAKE_CTEST_COMMAND} -V -C ${CMAKE_BUILD_TYPE}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../
  DEPENDS
    build_tests
)
