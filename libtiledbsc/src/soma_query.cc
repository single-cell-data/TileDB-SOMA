#include "tiledbsc/soma_query.h"
#include "tiledbsc/logger_public.h"

namespace tiledbsc {
using namespace tiledb;

#define X_ALLOC (1 << 22)

SOMAQuery::SOMAQuery(std::shared_ptr<SOMA> soma)
    : soma_(soma) {
    mq_obs_ = std::make_unique<ManagedQuery>(soma->open_array("obs"));
    mq_var_ = std::make_unique<ManagedQuery>(soma->open_array("var"));
    mq_x_ = std::make_unique<ManagedQuery>(soma->open_array("X/data"));
}

bool SOMAQuery::query(
    std::vector<std::string>& obs_attr_names,
    std::string_view obs_query_string,
    std::vector<std::string>& obs_ids,
    std::vector<std::string>& var_attr_names,
    std::string_view var_query_string,
    std::vector<std::string>& var_ids) {
    // TODO: Add QueryCondition
    (void)obs_query_string;
    (void)var_query_string;

    // Select obs columns to query
    mq_obs_->select_columns(obs_attr_names);

    // Select var columns to query
    mq_var_->select_columns(var_attr_names);

    // Select obs_id dimension values to query
    mq_obs_->select_points<std::string>("obs_id", obs_ids);

    // Select var_id dimension values to query
    mq_var_->select_points<std::string>("var_id", var_ids);

    int total_cells = 0;
    while (auto num_cells = mq_obs_->execute()) {
        LOG_DEBUG(fmt::format("*** obs cells = {}", num_cells));
        mq_x_->select_points("obs_id", mq_obs_->strings("obs_id"));
        total_cells += num_cells;
    }
    LOG_DEBUG(fmt::format("*** total obs cells = {}", total_cells));

    total_cells = 0;
    while (auto num_cells = mq_var_->execute()) {
        LOG_DEBUG(fmt::format("*** var cells = {}", num_cells));
        mq_x_->select_points("var_id", mq_var_->strings("var_id"));
        total_cells += num_cells;
    }
    LOG_DEBUG(fmt::format("*** total var cells = {}", total_cells));

    total_cells = 0;
    while (auto num_cells = mq_x_->execute()) {
        LOG_DEBUG(fmt::format("*** X cells = {}", num_cells));
        total_cells += num_cells;
    }
    LOG_DEBUG(fmt::format("*** total X cells = {}", total_cells));

    return true;
}

bool SOMAQuery::next_results() {
    return true;
}

}  // namespace tiledbsc
