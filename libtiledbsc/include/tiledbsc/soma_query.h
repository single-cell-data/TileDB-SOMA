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

constexpr size_t DEFAULT_INDEX_ALLOC = 1 << 20;  // 1 MiB
constexpr size_t DEFAULT_X_ALLOC = 1 << 26;      // 64 MiB

class SOMA;  // forward declaration

// clang-format off
/*
    def query(
        self,
        obs_attr_names: List[str] = [],
        obs_query_string: str = None,
        obs_ids: List[str] = None,
        var_attr_names: List[str] = [],
        var_query_string: str = None,
        var_ids: List[str] = None,
    ) -> Optional[SOMASlice]:
        """
        Subselects the obs, var, and X/data using the specified queries on obs and var,
        concatenating across SOMAs in the collection.  Queries use the TileDB-Py `QueryCondition`
        API. If `obs_query_string` is `None`, the `obs` dimension is not filtered and all of `obs`
        is used; similiarly for `var`. Return value of `None` indicates an empty slice.
        If `obs_ids` or `var_ids` are not `None`, they are effectively ANDed into the query.
        For example, you can pass in a known list of `obs_ids`, then use `obs_query_string`
        to further restrict the query.
        """

        soma_slices = []
        for soma in self:
            # E.g. querying for 'cell_type == "blood"' but this SOMA doesn't have a cell_type column in
            # its obs at all.
            if obs_query_string is not None and not soma.obs.has_attr_names(
                obs_attr_names
            ):
                continue
            # E.g. querying for 'feature_name == "MT-CO3"' but this SOMA doesn't have a feature_name
            # column in its var at all.
            if var_query_string is not None and not soma.var.has_attr_names(
                var_attr_names
            ):
                continue

            soma_slice = soma.query(
                obs_query_string=obs_query_string,
                var_query_string=var_query_string,
                obs_ids=obs_ids,
                var_ids=var_ids,
            )
            if soma_slice != None:
                # print("Slice SOMA from", soma.name, soma.X.data.shape(), "to", soma_slice.ann.X.shape)
                soma_slices.append(soma_slice)

        return SOMASlice.concat(soma_slices)


    def query(
        self,
        obs_query_string: Optional[str] = None,
        var_query_string: Optional[str] = None,
        obs_ids: Optional[List[str]] = None,
        var_ids: Optional[List[str]] = None,
    ) -> SOMASlice:
        """
        Subselects the SOMA's obs, var, and X/data using the specified queries on obs and var.
        Queries use the TileDB-Py `QueryCondition` API. If `obs_query_string` is `None`,
        the `obs` dimension is not filtered and all of `obs` is used; similiarly for `var`.
        """

        slice_obs_df = self.obs.query(query_string=obs_query_string, ids=obs_ids)
        # E.g. querying for 'cell_type == "blood"' and this SOMA does have a cell_type column in its
        # obs, but no rows with cell_type == "blood".
        if slice_obs_df is None:
            return None
        obs_ids = list(slice_obs_df.index)

        slice_var_df = self.var.query(query_string=var_query_string, ids=var_ids)
        # E.g. querying for 'feature_name == "MT-CO3"' and this SOMA does have a feature_name column
        # in its var, but no rows with feature_name == "MT-CO3".
        if slice_var_df is None:
            return None
        var_ids = list(slice_var_df.index)

        # TODO:
        # do this here:
        # * raw_var
        # do these in _assemble_soma_slice:
        # * raw_X
        # * obsm
        # * varm
        # * obsp
        # * varp

        return self._assemble_soma_slice(obs_ids, var_ids, slice_obs_df, slice_var_df)

    # ----------------------------------------------------------------
    def _assemble_soma_slice(
        self,
        obs_ids,
        var_ids,
        slice_obs_df,
        slice_var_df,
    ) -> SOMASlice:
        """
        An internal method for constructing a `SOMASlice` object given query results.
        """

        X = {key: self.X[key].dim_select(obs_ids, var_ids) for key in self.X.keys()}

        return SOMASlice(
            X=X,
            obs=slice_obs_df,
            var=slice_var_df,
        )

*/
// clang-format on

class SOMAQuery {
   public:
    /**
     * @brief Construct a new SOMAQuery object
     *
     * @param soma
     */
    SOMAQuery(
        SOMA* soma,
        size_t index_alloc = DEFAULT_INDEX_ALLOC,
        size_t x_alloc = DEFAULT_X_ALLOC);

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
    // Managed query for the obs array
    std::unique_ptr<ManagedQuery> mq_obs_;

    // Managed query for the var array
    std::unique_ptr<ManagedQuery> mq_var_;

    // Managed query for the X array
    std::unique_ptr<ManagedQuery> mq_x_;

    // Mutex to control access to mq_x_
    std::mutex mtx_;

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
