#ifndef TILEDBSC_H
#define TILEDBSC_H

// external includes
#include <map>
#include <optional>
#include <string>

// TODO fixes build error on VS2019 due to "missing" include in <tiledb/type.h>
#include <stdexcept>
#include <tiledb/tiledb>

#include <tiledbsc/common.h>
#include <tiledbsc/managed_query.h>

// local includes
#include "tiledbsc_export.h"

using namespace tiledb;

namespace tiledbsc {

// Placeholders
class SCRanges {};
using SCConfig = std::map<std::string, std::string>;

// Forward declaration
class SCGroup;

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

    std::shared_ptr<BufferGroup> buffers_;
};

/**
 * Class representing a single, sliceable (TileDB) array in an SCGroup.
 *
 */
class TILEDBSC_EXPORT SCArray {
    public:
        SCArray(std::shared_ptr<tiledb::Array> array);

        template <typename T>
        SCResult slice(std::string& dim_name, std::vector<T> ranges);

        static std::shared_ptr<SCArray> from_uri(tiledb::Context& ctx, std::string uri);

        std::shared_ptr<tiledb::Array> array() { return array_; };

    private:
        std::shared_ptr<tiledb::Array> array_;
};


/**
 * Class representing a collection of arrays comprising a single-cell data matrix.
 *
 * @details
 * Contains SCArray references and methods related to individual and joint (label-based) slicing
 * of arrays in the matrix.
 */
class TILEDBSC_EXPORT SCGroup {
    public:
        SCGroup(std::string& group_uri, std::optional<SCConfig> config = std::nullopt);

        std::string group_uri() { return uri_; };

        std::shared_ptr<SCArray> X() { return X_; };

    private:
         /* ctx_ *must* precede arrays due to destructor order: SCArray holds a reference */
        tiledb::Context ctx_;

        /* group URI */
        std::string uri_;
         // TODO this needs to support potential multiplicity of arrays
         //      for layers in separate arrays.

    public:
        // TODO these should become private
        std::shared_ptr<SCArray> X_;
        std::shared_ptr<SCArray> var_;
        std::shared_ptr<SCArray> obs_;
};

}; // end namespace tiledbsc

#endif // TILEDBSC_H