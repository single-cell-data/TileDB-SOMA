#include <tiledbsc/common.h>
#include <tiledbsc/query_result.h>
#include <tiledbsc/sc_arrowio.h>

namespace tiledbsc {

ResultBuffers& QueryResult::buffers() {
    return buffers_;
}

bool QueryResult::contains(const std::string& name) {
    return buffers_.find(name) != buffers_.end();
}

BufferSet& QueryResult::get(const std::string& name) {
    auto bufset = buffers_.find(name);
    if (bufset == buffers_.end()) {
        throw TileDBSCError("Buffer not found for " + name);
    }
    return bufset->second;
}

size_t QueryResult::nbuffers() {
    return buffers_.size();
}

std::vector<std::string> QueryResult::names() {
    std::vector<std::string> names;
    for (auto& [name, buf] : buffers_) {
        (void)buf;
        names.push_back(name);
    }
    return names;
}

tiledbsc::arrow::ArrowPair QueryResult::to_arrow(
    std::optional<std::string> name) {
    tiledbsc::arrow::ArrowAdapter adapter(*this);

    tiledbsc::arrow::ArrowPair res;

    if (name) {
        adapter.export_array(name.value().c_str(), res.array, res.schema);
    } else {
        adapter.export_table(&res.array, &res.schema);
    }

    res.disown();
    return res;
}

}  // end namespace tiledbsc