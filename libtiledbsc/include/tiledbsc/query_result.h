#ifndef TILEDBSC_QUERY_RESULT_H
#define TILEDBSC_QUERY_RESULT_H

#include <memory>
#include <vector>

#include <stdexcept> // required for TileDB, bug

#include <tiledb/tiledb>

#include <tiledbsc/buffer_set.h>
#include <tiledbsc/sc_arrowio.h>

#include "tiledbsc_export.h"


namespace tiledbsc {

//namespace arrow {
//    struct ArrowPair;
//}

using namespace std;

using ResultBuffers = map<string,BufferSet>;

class TILEDBSC_EXPORT QueryResult {
public:
    QueryResult() = delete;
    QueryResult(ResultBuffers&& buffers) : buffers_(std::move(buffers)) {};
    ResultBuffers& buffers();

    tiledbsc::arrow::ArrowPair to_arrow(std::optional<std::string>);

    bool contains(const std::string& name);

    BufferSet& get(const std::string& name);

    size_t nbuffers();

    std::vector<std::string> names();

private:
    ResultBuffers buffers_;
};

}; // namespace tiledbsc

#endif // TILEDBSC_QUERY_RESULT_H