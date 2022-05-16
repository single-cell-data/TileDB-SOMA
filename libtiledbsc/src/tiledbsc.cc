#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <tiledb/tiledb>

#include <tiledbsc/tiledbsc.h>
#include <tiledbsc/util.h>

using namespace tiledbsc::util;

using namespace std;

namespace tiledbsc {

SCGroup::SCGroup(std::string& group_uri, std::optional<SCConfig> config)
    : ctx_(tiledb::Context())
    , uri_(group_uri)  // TODO use config here
{
    std::string X_path = group_uri + "/X";
    std::string var_path = group_uri + "/var";
    std::string obs_path = group_uri + "/obs";

    check_paths_exist(vector<string>({X_path, var_path, obs_path}), config);

    // TODO
    // - only load schema here
    // - think to lazy-load the array on first access
    X_ = SCArray::from_uri(ctx_, X_path);
    var_ = SCArray::from_uri(ctx_, var_path);
    obs_ = SCArray::from_uri(ctx_, obs_path);
};

SCArray::SCArray(shared_ptr<tiledb::Array> array)
    : array_(std::move(array)) {
}

shared_ptr<SCArray> SCArray::from_uri(tiledb::Context& ctx, std::string uri) {
    auto array = std::make_shared<tiledb::Array>(ctx, uri, TILEDB_READ);

    return std::make_unique<SCArray>(std::move(array));
}

};  // end namespace tiledbsc