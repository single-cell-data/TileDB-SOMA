#ifndef SOMA_QUERY_H
#define SOMA_QUERY_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <future>
#include <mutex>

#include <tiledb/tiledb>

#include "tiledbsc/managed_query.h"
#include "tiledbsc/soma.h"

namespace tiledbsc {
using namespace tiledb;

class SOMA;  // forward declaration

class SOMAQuery {
   public:
    /**
     * @brief Construct a new SOMAQuery object
     *
     * @param soma SOMA
     * @param name SOMA name
     */
    SOMAQuery(SOMA* soma, std::string name);

    /**
     * @brief Select obs attributes to materialize.
     *
     * @param attr_names Vector of attribute names.
     */
    void select_obs_attrs(std::vector<std::string>& attr_names) {
        mq_obs_->select_columns(attr_names);
    }

    /**
     * @brief Select var attributes to materialize.
     *
     * @param attr_names Vector of attribute names.
     */
    void select_var_attrs(std::vector<std::string>& attr_names) {
        mq_var_->select_columns(attr_names);
    }

    /**
     * @brief Select obs_ids to include in the query.
     *
     * @param ids Vector of obs_id values.
     */
    void select_obs_ids(std::vector<std::string>& ids) {
        mq_obs_->select_points<std::string>("obs_id", ids);
    }

    /**
     * @brief Select var_ids to include in the query.
     *
     * @param ids Vector of var_id values.
     */
    void select_var_ids(std::vector<std::string>& ids) {
        mq_var_->select_points<std::string>("var_id", ids);
    }

    /**
     * @brief Set a query condition for the obs array query.
     *
     * @param qc TIleDB QueryCondition.
     */
    void set_obs_condition(QueryCondition& qc) {
        mq_obs_->set_condition(qc);
    }

    /**
     * @brief Set a query condition for the var array query.
     *
     * @param qc TIleDB QueryCondition.
     */
    void set_var_condition(QueryCondition& qc) {
        mq_var_->set_condition(qc);
    }

    template <class T>
    void set_obs_condition(const std::string& attr, T value, int op) {
        auto qc = QueryCondition::create(
            *ctx_, attr, value, (tiledb_query_condition_op_t)op);
        set_obs_condition(qc);
    }

    template <class T>
    void set_var_condition(const std::string& attr, T value, int op) {
        auto qc = QueryCondition::create(
            *ctx_, attr, value, (tiledb_query_condition_op_t)op);
        set_var_condition(qc);
    }

    /**
     * @brief Submit the query and return the first batch of results. To handle
     * incomplete queries, continue to call `next_results` until std::nullopt is
     * returned.
     *
     * @return std::optional<MultiArrayBuffers> Results or std::nullopt
     */
    std::optional<MultiArrayBuffers> next_results();

    std::optional<MultiArrayBuffers> results() {
        if (results_.empty()) {
            return std::nullopt;
        }
        return results_;
    }

    /**
     * @brief Return true if the SOMA query is complete.
     */
    bool is_complete() {
        return complete_;
    }

    /**
     * @brief Return the SOMA name.
     *
     * @return std::string SOMA name.
     */
    std::string name() {
        return name_;
    }

   private:
    // SOMA name
    std::string name_;

    // TileDB context
    std::shared_ptr<Context> ctx_;

    // Managed query for the obs array
    std::unique_ptr<ManagedQuery> mq_obs_;

    // Managed query for the var array
    std::unique_ptr<ManagedQuery> mq_var_;

    // Managed query for the X array
    std::unique_ptr<ManagedQuery> mq_x_;

    // Results of last call to next_results.
    MultiArrayBuffers results_;

    // Mutex to control access to mq_x_
    std::mutex mtx_;

    // If true, the query is empty because the obs or var query was empty
    bool empty_ = false;

    // If true, the query is complete
    bool complete_ = false;

    /**
     * @brief Submit a query (obs or var) and use the results to slice
     * the X query.
     *
     * @param mq Managed query for obs or var
     * @param dim_name "obs_id" or "var_id"
     */
    void query_and_select(
        std::unique_ptr<ManagedQuery>& mq, const std::string& dim_name);
};

}  // namespace tiledbsc

#endif
