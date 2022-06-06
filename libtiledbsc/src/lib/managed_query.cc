#include <stdexcept>  // TODO remove
#include <tiledb/tiledb>

#include <tiledbsc/buffer_set.h>
#include <tiledbsc/common.h>
#include <tiledbsc/managed_query.h>

using namespace std;
using namespace tiledbsc;

namespace tiledbsc {

/* ********************************* */
/*     Helper, private impl idiom    */
/* ********************************* */

class MQAux {
    // move alloc, execute, etc. here

    std::reference_wrapper<ManagedQuery> q;

   public:
    MQAux(ManagedQuery& q)
        : q(q) {
    }
    ~MQAux() = default;

    /* return true if this attribute or dimension should be skipped
       based on result restriction */
    bool skip_attrs_dims(const std::string& name) {
        return q.get()
            .use_attrs_dims_.value_or(std::set<std::string>())
            .count(name);
    }
};

/* ********************************* */
/*           ManagedQuery            */
/* ********************************* */

ManagedQuery::ManagedQuery(
    const std::shared_ptr<tiledb::Array> array, size_t initial_alloc)
    : array_(check_array(array)) {
    initial_alloc_ = initial_alloc;
    query_ = std::make_unique<tiledb::Query>(array->schema().context(), *array);

    if (array->schema().array_type() == TILEDB_SPARSE) {
        query_->set_layout(TILEDB_UNORDERED);
    } else {
        query_->set_layout(TILEDB_ROW_MAJOR);
    }

    impl_ = std::make_unique<MQAux>(*this);
}

ManagedQuery::~ManagedQuery() = default;

shared_ptr<tiledb::Array> ManagedQuery::check_array(
    std::shared_ptr<tiledb::Array> array) {
    if (array->query_type() != TILEDB_READ) {
        throw tiledb::TileDBError(
            "[ManagedQuery] Array not opened in read mode!");
    }

    return array;
}

void ManagedQuery::allocate_buffers() {
    auto schema = array_->schema();

    auto domain = schema.domain();
    for (const auto& dim : domain.dimensions()) {
        auto name = dim.name();

        if (impl_->skip_attrs_dims(name))
            continue;

        buffers_.emplace(name, BufferSet::from_dimension(dim, initial_alloc_));
    }

    for (const auto& [name, attr] : schema.attributes()) {
        if (impl_->skip_attrs_dims(name))
            continue;

        buffers_.emplace(name, BufferSet::from_attribute(attr, initial_alloc_));
    }
}

size_t ManagedQuery::initial_ncells() {
    /* TODO something smarter here */
    return std::min((size_t)10, initial_alloc_ / 5);
}

void ManagedQuery::set_buffers() {
    for (auto& [name, bg] : buffers_) {
        query_->set_data_buffer(
            name, (void*)bg.data_.data(), bg.data_.size() / bg.elem_nbytes());

        if (bg.isvar()) {
            auto buf = bg.offsets_.value();
            query_->set_offsets_buffer(name, (uint64_t*)buf.data(), buf.size());
        }

        if (bg.isnullable()) {
            auto validity = bg.validity_.value();
            query_->set_validity_buffer(
                name, (uint8_t*)validity.data(), validity.size());
        }
    }
}

void ManagedQuery::validate_query() {
}

void ManagedQuery::resize_result_buffers() {
    if (query_->query_status() != tiledb::Query::Status::COMPLETE) {
        throw tiledb::TileDBError(
            "internal error: attempted to resize result buffers but query not "
            "complete");
    }

    for (auto& [name, sizes] : query_->result_buffer_elements_nullable()) {
        auto [offsets_nelem, data_nelem, validity_nelem] = sizes;
        auto& bg = buffers_.at(name);

        size_t data_nbytes = data_nelem * bg.elem_nbytes();
        bg.data_.resize(data_nbytes);

        if (bg.isvar()) {
            bg.offsets_.value().resize(offsets_nelem);
        }

        if (bg.isnullable()) {
            bg.offsets_.value().resize(validity_nelem);
        }
    }
}

void ManagedQuery::complete_query() {
    query_->submit();

    // TODO run to completion with realloc

    resize_result_buffers();
}

std::unique_ptr<QueryResult> ManagedQuery::execute() {
    validate_query();

    allocate_buffers();

    set_buffers();

    complete_query();

    return make_unique<QueryResult>(std::move(buffers_));
}

};  // namespace tiledbsc