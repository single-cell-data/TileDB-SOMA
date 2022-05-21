#include <any>
#include <memory>
#include <stdexcept>
#include <vector>

#include <tiledb/tiledb>

#include <tiledbsc/buffer_set.h>
#include <tiledbsc/query_result.h>
#include <tiledbsc/sc_arrowio.h>
#include <tiledbsc/util.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace tiledb;
using namespace tiledbsc;
using namespace tiledbsc::arrow;

TEST_CASE("SCArrow int export", "[arrow][export]") {
    auto b1 = tiledbsc::BufferSet::alloc(
        "foo", TILEDB_INT32, 1023, 4, false, false);

    ArrowPair ap;
    ArrowAdapter::export_buffer(b1, ap.array, ap.schema);

    REQUIRE(ap.array->length == 1023);
    REQUIRE(ap.array->null_count == 0);
    REQUIRE(ap.array->offset == 0);
    REQUIRE(ap.array->n_buffers == 2);
    REQUIRE(ap.array->n_children == 0);
    REQUIRE(ap.array->children == nullptr);
    REQUIRE(ap.array->dictionary == nullptr);

    REQUIRE(std::string_view(ap.schema->format) == "i");
    REQUIRE(std::string_view(ap.schema->name) == "foo");
    REQUIRE(ap.schema->metadata == nullptr);
    REQUIRE(ap.schema->flags == 0);
    REQUIRE(ap.schema->n_children == 0);
};

TEST_CASE("SCArrow strings export", "[arrow][export]") {
    std::vector<string> data_orig{
        "", "abcd", "ef", "ghijk", "lmno", "p", "", "q", "rstu", "vwxyz", ""};

    auto&& [data_buf, offsets_buf] = util::to_varlen_buffers(data_orig);

    auto b = BufferSet::from_data(
        "strings", TILEDB_STRING_ASCII, 1, data_buf, offsets_buf);

    ArrowPair ap;
    ArrowAdapter::export_buffer(b, ap.array, ap.schema);

    REQUIRE(ap.array->length == 11);
    REQUIRE(ap.array->null_count == 0);
    REQUIRE(ap.array->offset == 0);
    REQUIRE(ap.array->n_buffers == 3);
    REQUIRE(ap.array->n_children == 0);
    REQUIRE(ap.array->children == nullptr);
    REQUIRE(ap.array->dictionary == nullptr);

    REQUIRE(std::string_view(ap.schema->format) == "U");
    REQUIRE(std::string_view(ap.schema->name) == "strings");
    REQUIRE(ap.schema->metadata == nullptr);
    REQUIRE(ap.schema->flags == 0);
    REQUIRE(ap.schema->n_children == 0);
};

TEST_CASE("SCArrow nullable strings export", "[arrow][export]") {
    std::vector<string> data_orig{
        "", "abcd", "ef", "ghijk", "lmno", "p", "", "q", "rstu", "vwxyz", ""};

    auto&& [data_buf, offsets_buf] = util::to_varlen_buffers(data_orig);

    std::vector<uint8_t> validity_vec_tmp{0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1};
    std::vector<byte> validity_buf(validity_vec_tmp.size());
    std::transform(
        validity_vec_tmp.begin(),
        validity_vec_tmp.end(),
        validity_buf.begin(),
        [](uint8_t v) -> byte { return (byte)v; });

    auto b = BufferSet::from_data(
        "strings", TILEDB_STRING_ASCII, 1, data_buf, offsets_buf, validity_buf);

    ArrowPair ap;
    ArrowAdapter::export_buffer(b, ap.array, ap.schema);

    REQUIRE(ap.array->length == 11);
    REQUIRE(ap.array->null_count == 7);
    REQUIRE(ap.array->offset == 0);
    REQUIRE(ap.array->n_buffers == 3);
    REQUIRE(ap.array->n_children == 0);
    REQUIRE(ap.array->children == nullptr);
    REQUIRE(ap.array->dictionary == nullptr);

    REQUIRE(std::string_view(ap.schema->format) == "U");
    REQUIRE(std::string_view(ap.schema->name) == "strings");
    REQUIRE(ap.schema->metadata == nullptr);
    // REQUIRE(ap.schema->flags == 0);
    REQUIRE(ap.schema->n_children == 0);
};
