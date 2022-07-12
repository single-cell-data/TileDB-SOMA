#ifndef SOMA_COLLECTION_QUERY_H
#define SOMA_COLLECTION_QUERY_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <future>
#include <mutex>

#include <tiledb/tiledb>

#include "tiledbsc/managed_query.h"
#include "tiledbsc/soma.h"

namespace tiledbsc {
using namespace tiledb;

class SOMA;            // forward declaration
class SOMAQuery;       // forward declaration
class SOMACollection;  // forward declaration

// 1. Add all SOMAs in the SOCO to a map: SOMA name -> SOMAQuery
// 2. Apply attrs, ids, and query conditions to each SOMAQuery
// 3. next_results(): Run all SOMAQuery.next_results() in parallel, return
// results.

// TODO: Pass alloc sizes to SOMAQuery

class SOMACollectionQuery {
   private:
    static const size_t DEFAULT_INDEX_ALLOC = 1 << 20;  // 1M cells
    static const size_t DEFAULT_X_ALLOC = 1 << 26;      // 64M cells

   public:
    /**
     * @brief Construct a new SOMACollectionQuery object with all SOMAs in the
     * provided SOMACollection.
     *
     * @param soco SOMACollection
     * @param index_alloc Number of cells allocated for obs and var
     * @param x_alloc Number of cells allocated for X/data
     */
    SOMACollectionQuery(
        SOMACollection* soco,
        size_t index_alloc = DEFAULT_INDEX_ALLOC,
        size_t x_alloc = DEFAULT_X_ALLOC);

    /**
     * @brief Select obs attributes to materialize.
     *
     * @param attr_names Vector of attribute names.
     */
    void select_obs_attrs(std::vector<std::string>& attr_names) {
        for (auto& [name, sq] : soma_queries_) {
            sq->select_obs_attrs(attr_names);
        }
    }

    /**
     * @brief Select var attributes to materialize.
     *
     * @param attr_names Vector of attribute names.
     */
    void select_var_attrs(std::vector<std::string>& attr_names) {
        for (auto& [name, sq] : soma_queries_) {
            sq->select_var_attrs(attr_names);
        }
    }

    /**
     * @brief Select obs_ids to include in the query.
     *
     * @param ids Vector of obs_id values.
     */
    void select_obs_ids(std::vector<std::string>& ids) {
        for (auto& [name, sq] : soma_queries_) {
            sq->select_obs_ids(ids);
        }
    }

    /**
     * @brief Select var_ids to include in the query.
     *
     * @param ids Vector of var_id values.
     */
    void select_var_ids(std::vector<std::string>& ids) {
        for (auto& [name, sq] : soma_queries_) {
            sq->select_var_ids(ids);
        }
    }

    /**
     * @brief Set a query condition for the obs array query.
     *
     * @param qc TIleDB QueryCondition.
     */
    void set_obs_condition(QueryCondition& qc) {
        for (auto& [name, sq] : soma_queries_) {
            sq->set_obs_condition(qc);
        }
    }

    /**
     * @brief Set a query condition for the var array query.
     *
     * @param qc TIleDB QueryCondition.
     */
    void set_var_condition(QueryCondition& qc) {
        for (auto& [name, sq] : soma_queries_) {
            sq->set_var_condition(qc);
        }
    }

    /**
     * @brief Submit the query and return the first batch of results. To handle
     * incomplete queries, continue to call `next_results` until std::nullopt is
     * returned.
     *
     * @return std::optional<
     * std::unordered_map<std::string, std::shared_ptr<ColumnBuffer>>> Results
     * or std::nullopt if the query is complete.
     */
    std::optional<
        std::unordered_map<std::string, std::shared_ptr<ColumnBuffer>>>
    next_results();

   private:
    // Map of SOMA name -> SOMAQuery
    std::unordered_map<std::string, std::unique_ptr<SOMAQuery>> soma_queries_;

    // If true, the query has been submitted
    bool submitted_ = false;

    // Number of threads in the thread pool
    size_t threads_ = 16;
};

}  // namespace tiledbsc

#endif
