#include <stdexcept>  // required for TileDB, bug

#include <tiledb/tiledb>

#include <tiledbsc/buffer_set.h>
#include <tiledbsc/common.h>
#include <tiledbsc/managed_query.h>

#include <tiledbsc/ij_query.h>
#include <string>

using namespace std;
using namespace tiledbsc;

// we need to support two basic classes of types:
// - strings
// - scalars

namespace tiledbsc {

IJQuery::IJQuery(
    const shared_ptr<tiledb::Array> reference,
    pair<string, string> ref_pair,
    const shared_ptr<tiledb::Array> target,
    string dest_name)
    : array_ref_(reference)
    , ref_dim_(ref_pair.first)
    , ref_attr_(ref_pair.second)
    , array_tgt_(target)
    , tgt_dim_(dest_name) {
    check_compatible();
}

IJQuery::~IJQuery() = default;

void IJQuery::check_compatible() {
    auto schema1 = array_ref_->schema();

    if (!schema1.domain().has_dimension(ref_dim_)) {
        throw TileDBSCError(
            "Reference array does not have dimension " + ref_dim_);
    }
    if (!schema1.has_attribute(ref_attr_)) {
        throw TileDBSCError(
            "Reference array does not have attribute " + ref_attr_);
    }

    auto schema2 = array_tgt_->schema();
    if (!schema2.domain().has_dimension(tgt_dim_)) {
        throw TileDBSCError("Target array does not have attribute " + tgt_dim_);
    }
}

unique_ptr<QueryResult> IJQuery::select_from_points(generic_typed_vector v) {
    (void)v;
    return make_unique<QueryResult>(ResultBuffers());
};

};  // namespace tiledbsc