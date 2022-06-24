#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_predicate.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <random>

#include <tiledbsc/util.h>
#include <tiledb/tiledb>
#include <tiledbsc/tiledbsc>

using namespace tiledb;
using namespace tiledbsc;
using namespace Catch::Matchers;

#ifndef TILEDBSC_SOURCE_ROOT
#define TILEDBSC_SOURCE_ROOT "not_defined"
#endif

const std::string src_path = TILEDBSC_SOURCE_ROOT;

namespace {
static auto create_array(const std::string& uri, Context& ctx) {
    // Delete array if it exists
    auto vfs = VFS(ctx);
    if (vfs.is_dir(uri)) {
        vfs.remove_dir(uri);
    }

    // Create schema
    ArraySchema schema(ctx, TILEDB_SPARSE);
    auto dim = Dimension::create(
        ctx, "d0", TILEDB_STRING_ASCII, nullptr, nullptr);
    dim.set_cell_val_num(TILEDB_VAR_NUM);

    Domain domain(ctx);
    domain.add_dimension(dim);
    schema.set_domain(domain);

    auto attr = Attribute::create<int>(ctx, "a0");
    schema.add_attribute(attr);
    schema.check();

    // Create array and open for writing
    Array::create(uri, schema);
    Array array(ctx, uri, TILEDB_WRITE);

    std::vector<std::string> d0 = {
        "a", "bb", "ccc", "dddd", "eeeee", "fffffff"};
    auto [d0_data, d0_offsets] = util::to_varlen_buffers(d0, false);

    std::vector<int> a0(d0.size());
    std::generate(a0.begin(), a0.end(), std::default_random_engine());

    // Write data to array and close the array
    Query query(ctx, array);
    query.set_layout(TILEDB_UNORDERED)
        .set_buffer("d0", d0_data)
        .set_offsets_buffer("d0", d0_offsets)
        .set_buffer("a0", a0);
    query.submit();
    array.close();

    // Open the array for reading and return a shared pointer
    return std::tuple(std::make_shared<Array>(ctx, uri, TILEDB_READ), d0, a0);
}

};  // namespace

TEST_CASE("ManagedQuery: Basic execution test") {
    std::string uri = "mem://unit-test-array";
    auto ctx = Context();
    auto [array, d0, a0] = create_array(uri, ctx);

    auto mq = ManagedQuery(array);
    mq.select_columns({"d0", "a0"});

    auto num_cells = mq.execute();

    REQUIRE(num_cells == d0.size());
    REQUIRE_THAT(d0, Equals(mq.strings("d0")));
    REQUIRE_THAT(a0, Equals(util::to_vector(mq.data<int>("a0"))));
}
