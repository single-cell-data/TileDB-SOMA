############################################################
# Definitions
############################################################

# Set root source directory as compiler definition for use in tests
get_filename_component(TILEDBSC_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../" ABSOLUTE)
add_compile_definitions(TILEDBSC_SOURCE_ROOT="${TILEDBSC_SOURCE_ROOT}")

############################################################
# Dependencies
############################################################

find_package(TileDB_EP REQUIRED)

############################################################
# SC unit test
############################################################

# TODO make this conditional
if (TRUE)
    find_package(Catch_EP REQUIRED)

    add_executable(unit_sc EXCLUDE_FROM_ALL)
    target_link_libraries(unit_sc
        PUBLIC
        Catch2::Catch2WithMain
        TileDB::tiledb_shared
    )

    target_link_libraries(unit_sc
        PRIVATE
        $<TARGET_OBJECTS:buffer_set>
        $<TARGET_OBJECTS:managed_query>
    )

    target_include_directories(unit_sc
       PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/../include/tiledbsc
        ../include)

    # Sources for tests
    target_sources(unit_sc PUBLIC
        unit_sc.cc
    )

    target_compile_definitions(unit_sc PRIVATE CATCH_CONFIG_MAIN)

    add_test(
        NAME "unit_sc"
        COMMAND $<TARGET_FILE:unit_sc>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

