############################################################
# Definitions
############################################################

# Set root source directory as compiler definition for use in tests
get_filename_component(TILEDBSOMA_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../" ABSOLUTE)
add_compile_definitions(TILEDBSOMA_SOURCE_ROOT="${TILEDBSOMA_SOURCE_ROOT}")


############################################################
# SOMA unit test
############################################################

find_package(Catch_EP REQUIRED)

add_executable(unit_soma
    $<TARGET_OBJECTS:TILEDB_SOMA_OBJECTS>
    $<TARGET_OBJECTS:TILEDBSOMA_NANOARROW_OBJECT>
    common.cc
    common.h
    unit_arrow_adapter.cc
    unit_column_buffer.cc
    unit_managed_query.cc
    unit_soma_array.cc
    unit_soma_group.cc
    unit_soma_dataframe.cc
    unit_soma_dense_ndarray.cc
    unit_soma_sparse_ndarray.cc
    unit_soma_collection.cc
    unit_soma_scene.cc
    unit_soma_point_cloud_dataframe.cc
    unit_soma_multiscale_image.cc
    test_indexer.cc
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
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/external/include
    ${pybind11_INCLUDE_DIRS}
    ${TILEDB_SOMA_EXPORT_HEADER_DIR}
    $<TARGET_PROPERTY:spdlog::spdlog,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(unit_soma
  PRIVATE
  CATCH_CONFIG_MAIN
)
target_compile_options(unit_soma
  PRIVATE
    ${TILEDBSOMA_COMPILE_OPTIONS}
    ${TILEDBSOMA_WERROR_OPTION}
)

if (NOT MSVC)
  # Allow deprecated function for writing to an array with a timestamp
  target_compile_options(unit_soma PRIVATE -Wno-deprecated-declarations)
endif()

if (NOT APPLE AND NOT WIN32)
    target_link_libraries(unit_soma PRIVATE pthread)
endif()

add_test(
    NAME "unit_soma"
    COMMAND $<TARGET_FILE:unit_soma> "--durations=yes"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

if (WIN32)
  add_custom_command(TARGET unit_soma POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:unit_soma> $<TARGET_FILE_DIR:unit_soma>
    COMMAND_EXPAND_LISTS
  )
endif()

add_custom_target(build_tests)
add_dependencies(build_tests
  unit_soma
)

############################################################
# make check
############################################################

# Add custom target 'check' in the root of the tree
add_custom_target(
  check COMMAND ${CMAKE_CTEST_COMMAND} -V -C $<CONFIG>
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../
  DEPENDS
    build_tests
)
