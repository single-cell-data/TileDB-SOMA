############################################################
# common definitions for tests
############################################################

# Set root source directory as compiler definition for use in tests
get_filename_component(TILEDBSC_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../../" ABSOLUTE)
add_compile_definitions(TILEDBSC_SOURCE_ROOT="${TILEDBSC_SOURCE_ROOT}")

############################################################
# util object library
############################################################

add_library(util OBJECT
    util.cc
)

target_link_libraries(util
    PRIVATE
    TileDB::tiledb_shared
)

target_include_directories(util
  PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/../
    ../../include)

target_compile_definitions(util PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)

if (TILEDBSC_TESTING)
    ## Test for util
    add_executable(unit_util)
    #    $<TARGET_OBJECTS:util>
    #)

    target_sources(unit_util PUBLIC
        test/unit_util.cc
    )

    target_link_libraries(unit_util
      PRIVATE
        Catch2::Catch2WithMain
        TileDB::tiledb_shared
        util
    )

    target_compile_definitions(unit_util PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)
    target_compile_definitions(unit_util PRIVATE CATCH_CONFIG_MAIN)

    add_test(
        NAME "unit_util"
        COMMAND $<TARGET_FILE:unit_util>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

############################################################
# buffer set object library
############################################################

add_library(buffer_set OBJECT
    buffer_set.cc
)
add_library(TileDB::buffer_set ALIAS buffer_set)

target_link_libraries(buffer_set
  PUBLIC
    TileDB::tiledb_shared
    util
    )

target_include_directories(buffer_set
    PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/../
    ../../include)

target_compile_definitions(buffer_set PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)

# TODO make this conditional
if (TILEDBSC_TESTING)
    find_package(Catch_EP REQUIRED)

    add_executable(unit_buffer_set)

    # Sources for tests
    target_sources(unit_buffer_set
        PUBLIC
        test/unit_buffer_set.cc
    )

    target_link_libraries(unit_buffer_set
      PRIVATE
        Catch2::Catch2WithMain
        util
        buffer_set
    )

    target_compile_definitions(unit_buffer_set PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)
    target_compile_definitions(unit_buffer_set PRIVATE CATCH_CONFIG_MAIN)

    add_test(
        NAME "unit_buffer_set"
        COMMAND $<TARGET_FILE:unit_buffer_set>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

############################################################
# scarrow object library
############################################################

add_library(sc_arrow OBJECT
    sc_arrow_export.cc
)
add_library(TileDB::sc_arrow ALIAS sc_arrow)

target_link_libraries(sc_arrow
  PUBLIC
    TileDB::tiledb_shared
    buffer_set
    )


target_include_directories(sc_arrow
    PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/../
    ../../include)

target_compile_definitions(sc_arrow PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)

# TODO make this conditional
if (TILEDBSC_TESTING)
    find_package(Catch_EP REQUIRED)

    add_executable(unit_sc_arrow)
    target_link_libraries(unit_sc_arrow
        PUBLIC
        Catch2::Catch2WithMain
        TileDB::tiledb_shared
        util
        sc_arrow
        buffer_set
        query_result
    )
    #target_link_libraries(unit_sc_arrow
    #    PRIVATE
    #)

    # Sources for tests
    target_sources(unit_sc_arrow PUBLIC
        test/unit_sc_arrow.cc
    )

    target_compile_definitions(unit_sc_arrow PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)
    target_compile_definitions(unit_sc_arrow PRIVATE CATCH_CONFIG_MAIN)

    add_test(
        NAME "unit_sc_arrow"
        COMMAND $<TARGET_FILE:unit_sc_arrow>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    ############################################################
    # pybind11 test harness
    #include(py_query_result_aux.cmake)
    ############################################################

    add_dependencies(unit_sc_arrow tiledbsc)
endif()

############################################################
# query result object library
############################################################

add_library(query_result OBJECT
    query_result.cc
)
add_library(TileDB::query_result ALIAS query_result)

target_link_libraries(query_result
    PUBLIC
    TileDB::tiledb_shared
    buffer_set
    sc_arrow
    )

target_include_directories(query_result
    PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/../
    ../../include)

target_compile_definitions(query_result PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)

# TODO make this conditional
if (TILEDBSC_TESTING)
    find_package(Catch_EP REQUIRED)

    add_executable(unit_query_result)
    target_link_libraries(unit_query_result
        PUBLIC
        Catch2::Catch2WithMain
        TileDB::tiledb_shared
        util
        buffer_set
        query_result
        sc_arrow
    )

    # Sources for tests
    target_sources(unit_query_result PUBLIC
        test/unit_query_result.cc
    )

    target_compile_definitions(unit_query_result PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)
    target_compile_definitions(unit_query_result PRIVATE CATCH_CONFIG_MAIN)

    add_test(
        NAME "unit_query_result"
        COMMAND $<TARGET_FILE:unit_query_result>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    ############################################################
    # pybind11 test harness
    include(py_query_result_aux.cmake)
    ############################################################

    add_dependencies(unit_query_result tiledbsc)
endif()


############################################################
# managed query object library
############################################################

add_library(managed_query OBJECT
    managed_query.cc
)
add_library(TileDB::managed_query ALIAS managed_query)

target_link_libraries(managed_query
    TileDB::tiledb_shared
    #$<TARGET_OBJECTS:buffer_set>
    #$<TARGET_OBJECTS:query_result>
    #$<TARGET_OBJECTS:sc_arrow>
    buffer_set
    sc_arrow
    query_result
    )

target_include_directories(managed_query
    PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/../
    ../../include)

target_compile_definitions(managed_query PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)

if (TILEDBSC_TESTING)
    find_package(Catch_EP REQUIRED)

    add_executable(unit_managed_query EXCLUDE_FROM_ALL)
    target_link_libraries(unit_managed_query
        PUBLIC
        Catch2::Catch2WithMain
        TileDB::tiledb_shared
        util
        buffer_set
        sc_arrow
        query_result
        managed_query
    )

    # Sources for code under test
    target_sources(unit_managed_query PUBLIC ${SOURCES})

    # Sources for tests
    target_sources(unit_managed_query PUBLIC
        test/unit_managed_query.cc
    )

    target_compile_definitions(unit_managed_query PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)
    target_compile_definitions(unit_managed_query PRIVATE CATCH_CONFIG_MAIN)

    add_test(
        NAME "unit_managed_query"
        COMMAND $<TARGET_FILE:unit_managed_query>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_dependencies(unit_managed_query tiledbsc)
endif()


############################################################
# array to array query object library
############################################################

add_library(ij_query OBJECT
    ij_query.cc
)
add_library(TileDB::ij_query ALIAS ij_query)

target_link_libraries(ij_query
    TileDB::tiledb_shared
    buffer_set
    managed_query
    )

target_include_directories(ij_query
    PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/../
    ../../include)


target_compile_definitions(ij_query PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)

# TODO make this conditional
if (TILEDBSC_TESTING)
    find_package(Catch_EP REQUIRED)

    add_executable(unit_ij_query)
    target_link_libraries(unit_ij_query
        PUBLIC
        Catch2::Catch2WithMain
        TileDB::tiledb_shared
        ij_query
        buffer_set
        managed_query
    )

    # Sources for code under test
    target_sources(unit_ij_query PUBLIC ${SOURCES})

    # Sources for tests
    target_sources(unit_ij_query PUBLIC
        test/unit_ij_query.cc
    )

    target_compile_definitions(unit_ij_query PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)
    target_compile_definitions(unit_ij_query PRIVATE CATCH_CONFIG_MAIN)

    add_test(
        NAME "unit_ij_query"
        COMMAND $<TARGET_FILE:unit_ij_query>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

endif()