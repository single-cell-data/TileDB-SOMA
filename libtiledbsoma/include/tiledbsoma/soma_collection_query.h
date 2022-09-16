/**
 * @file   soma_collection_query.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 *   This declares the soma collect API
 */

#ifndef SOMA_COLLECTION_QUERY_H
#define SOMA_COLLECTION_QUERY_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <future>
#include <mutex>

#include <tiledb/tiledb>

#include "tiledbsoma/column_buffer.h"
#include "tiledbsoma/soma.h"

namespace tiledbsoma {
using namespace tiledb;

class SOMA;            // forward declaration
class SOMAQuery;       // forward declaration
class SOMACollection;  // forward declaration

// 1. Add all SOMAs in the SOCO to a map: SOMA name -> SOMAQuery
// 2. Apply attrs, ids, and query conditions to each SOMAQuery
// 3. next_results(): Run all SOMAQuery.next_results() in parallel, return
// results.

class SOMACollectionQuery {
   public:
    /**
     * @brief Construct a new SOMACollectionQuery object with all SOMAs in the
     * provided SOMACollection.
     *
     * @param soco SOMACollection
     * @param index_alloc Number of cells allocated for obs and var
     * @param x_alloc Number of cells allocated for X/data
     */
    SOMACollectionQuery(SOMACollection* soco);

    /**
     * @brief Select obs attributes to materialize.
     *
     * @param attr_names Vector of attribute names.
     */
    void select_obs_attrs(std::vector<std::string>& attr_names) {
        for (auto& sq : soma_queries_) {
            sq->select_obs_attrs(attr_names);
        }
    }

    /**
     * @brief Select var attributes to materialize.
     *
     * @param attr_names Vector of attribute names.
     */
    void select_var_attrs(std::vector<std::string>& attr_names) {
        for (auto& sq : soma_queries_) {
            sq->select_var_attrs(attr_names);
        }
    }

    /**
     * @brief Select obs_ids to include in the query.
     *
     * @param ids Vector of obs_id values.
     */
    void select_obs_ids(std::vector<std::string>& ids) {
        for (auto& sq : soma_queries_) {
            sq->select_obs_ids(ids);
        }
    }

    /**
     * @brief Select var_ids to include in the query.
     *
     * @param ids Vector of var_id values.
     */
    void select_var_ids(std::vector<std::string>& ids) {
        for (auto& sq : soma_queries_) {
            sq->select_var_ids(ids);
        }
    }

    /**
     * @brief Set a query condition for the obs array query.
     *
     * @param qc TIleDB QueryCondition.
     */
    void set_obs_condition(QueryCondition& qc) {
        for (auto& sq : soma_queries_) {
            sq->set_obs_condition(qc);
        }
    }

    /**
     * @brief Set a query condition for the var array query.
     *
     * @param qc TIleDB QueryCondition.
     */
    void set_var_condition(QueryCondition& qc) {
        for (auto& sq : soma_queries_) {
            sq->set_var_condition(qc);
        }
    }

    /**
     * @brief Submit the query and return the first batch of results. To handle
     * incomplete queries, continue to call `next_results` until std::nullopt is
     * returned.
     *
     * @return std::optional<MultiArrayBuffers> Results or std::nullopt
     */
    std::optional<MultiArrayBuffers> next_results();

   private:
    // Map of SOMA name -> SOMAQuery
    std::vector<std::unique_ptr<SOMAQuery>> soma_queries_;

    // If true, the query has been submitted
    bool submitted_ = false;

    // Number of threads in the thread pool
    size_t threads_ = 16;
};

}  // namespace tiledbsoma

#endif
