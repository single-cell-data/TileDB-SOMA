#ifndef TILEDBSC_IJQUERY_H
#define TILEDBSC_IJQUERY_H

#include "tiledbsc_export.h"

#include <memory>
#include <vector>

#include <tiledb/tiledb>

#include <tiledbsc/managed_query.h>
#include <tiledbsc/buffer_set.h>
#include <tiledbsc/generic_typed_vector.h>


using namespace tiledbsc;
using namespace tiledb;
using namespace std;

namespace tiledbsc {

/* ********************************* */
/*          IJQuery             */
/* ********************************* */

/**
 *
 * IJQuery returns the result of:
 * 1) query source array with selection
 * 2) use result of (1) to query dest array
 * 3) return result of (2)
 *
 */

class TILEDBSC_EXPORT IJQuery {
public:
    IJQuery(
        const shared_ptr<tiledb::Array> reference,
        pair<string,string> ref_pair,
        const shared_ptr<tiledb::Array> target,
        string dest_name
    );

    ~IJQuery();

    unique_ptr<QueryResult> select_from_points(generic_typed_vector v);

private:
    void check_compatible();

private:
    const shared_ptr<tiledb::Array> array_ref_;
    const string ref_dim_;
    const string ref_attr_;
    const shared_ptr<tiledb::Array> array_tgt_;
    std::string tgt_dim_;
};

}; // namespace tiledbsc

#endif // TILEDBSC_IJQUERY_H