#ifndef TILEDBSC_H
#define TILEDBSC_H

// local includes
#include "tiledbsc_export.h"

// external includes
#include <map>
#include <optional>
#include <string>

// TODO fixes build error on VS2019 due to "missing" include in <tiledb/type.h>
#include <stdexcept>
#include <tiledb/tiledb>

#include <tiledbsc/common.h>
#include <tiledbsc/managed_query.h>


namespace tiledbsc {

using namespace tiledb;

// Placeholders
class SCRanges {};

// Forward declaration
class ArrowTable;

/**
 * Struct representing a set of results from a query.
 *
 * @details
 * The object contains a set of buffers resulting from a generic
 * TileDB query against an array.
 */
struct SCResult {
    /**
     * @brief Returns the BufferGroup mapping dimension and attribute
     *        results to a BufferSet (buffer composition depends on
     *        attribute type).
     *
     * @throws TileDBSCError if the array is already open or other error occurred.
     */
    std::shared_ptr<QueryResult> query_result();


    /*** Arrow API ***/
    /**
     * @brief Returns the BufferGroup mapping dimension and attribute
     *        results to a BufferSet (buffer composition depends on
     *        attribute type).
     *
     * @throws TileDBSCError if the array is already open or other error occurred.
     */
    void to_arrow(void* schema, void* array);


    /**
     * @brief Returns the BufferGroup mapping dimension and attribute
     *        results to a BufferSet (buffer composition depends on
     *        attribute type).
     *
     * @throws TileDBSCError if the array is already open or other error occurred.
     */
    void release_callback(void* user_data);


    private:
    /* ********************************* */
    /*         PRIVATE ATTRIBUTES        */
    /* ********************************* */

    ResultBuffers buffers_;
};


}; // end namespace tiledbsc

#endif // TILEDBSC_H