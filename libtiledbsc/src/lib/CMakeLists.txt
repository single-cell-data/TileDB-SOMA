############################################################
# common definitions for tests
############################################################

# Set root source directory as compiler definition for use in tests
get_filename_component(TILEDBSC_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../../" ABSOLUTE)
add_compile_definitions(TILEDBSC_SOURCE_ROOT="${TILEDBSC_SOURCE_ROOT}")

############################################################
# buffer set object library
############################################################

add_library(buffer_set OBJECT
    buffer_set.cc
)
add_library(TileDB::buffer_set ALIAS buffer_set)

target_link_libraries(buffer_set TileDB::tiledb_shared)
target_include_directories(buffer_set
    PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/../
    ../../include)

target_compile_definitions(buffer_set PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)

# TODO make this conditional
if (FALSE)
    find_package(Catch_EP REQUIRED)

    add_executable(unit_buffer_set EXCLUDE_FROM_ALL)
    target_link_libraries(unit_buffer_set
        PUBLIC
        Catch2::Catch2WithMain
        TileDB::tiledb_shared
        buffer_set
    )

    # Sources for code under test
    target_sources(unit_buffer_set PUBLIC ${SOURCES})

    # Sources for tests
    target_sources(unit_buffer_set PUBLIC
        test/unit_buffer_set.cc
    )

    target_compile_definitions(unit_buffer_set PRIVATE CATCH_CONFIG_MAIN)

    add_test(
        NAME "unit_buffer_set"
        COMMAND $<TARGET_FILE:unit_buffer_set>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
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
    $<TARGET_OBJECTS:buffer_set>)

target_include_directories(managed_query
    PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/../
    ../../include)

target_compile_definitions(managed_query PRIVATE -DTILEDB_SC_OBJECTS_EXPORTS)

# TODO make this conditional
if (FALSE)
    find_package(Catch_EP REQUIRED)

    add_executable(unit_managed_query EXCLUDE_FROM_ALL)
    target_link_libraries(unit_managed_query
        PUBLIC
        Catch2::Catch2WithMain
        TileDB::tiledb_shared
        managed_query
    )

    # Sources for code under test
    target_sources(unit_managed_query PUBLIC ${SOURCES})

    # Sources for tests
    target_sources(unit_managed_query PUBLIC
        test/unit_managed_query.cc
    )

    target_compile_definitions(unit_managed_query PRIVATE CATCH_CONFIG_MAIN)

    add_test(
        NAME "unit_managed_query"
        COMMAND $<TARGET_FILE:unit_managed_query>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()