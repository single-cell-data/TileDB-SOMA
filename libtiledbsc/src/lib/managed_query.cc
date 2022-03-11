#include <stdexcept> // TODO remove
#include <tiledb/tiledb>

#include <tiledbsc/managed_query.h>
#include <tiledbsc/buffer_set.h>
#include <tiledbsc/common.h>

using namespace std;

namespace tiledbsc {

/* ********************************* */
/*           QueryResult             */
/* ********************************* */

//QueryResult::to_arrow

std::shared_ptr<BufferGroup> QueryResult::buffers() {
    return buffers_;
}

std::optional<std::shared_ptr<BufferSet>> QueryResult::get(std::string name) {
    auto bufset = buffers_->buffers.find(name);
    if (bufset == buffers_->buffers.end()) {
        return std::nullopt;
    }
    return std::make_optional(bufset->second);
}

size_t QueryResult::nbuffers() {
    return buffers_->buffers.size();
}

std::vector<std::string> QueryResult::names() {
    std::vector<std::string> names;
    for (auto& [name, buf] : buffers_->buffers) {
        (void)buf;
        names.push_back(name);
    }
    return names;
}

/* ********************************* */
/*     Helper, private impl idiom    */
/* ********************************* */

class MQAux {
    // move alloc, execute, etc. here

    std::reference_wrapper<ManagedQuery> q;

    public:
    MQAux(ManagedQuery& q) : q(q) {}
    ~MQAux() = default;

    /* return true if this attribute or dimension should be skipped
       based on result restriction */
    bool skip_attrs_dims(const std::string& name) {
        return q.get().use_attrs_dims_.value_or(std::set<std::string>()).count(name);
    }

};

/* ********************************* */
/*           ManagedQuery            */
/* ********************************* */

ManagedQuery::ManagedQuery(
    const std::shared_ptr<tiledb::Array> array,
    size_t initial_alloc
) : array_(check_array(array))
{
    initial_alloc_ = initial_alloc;
    query_ = std::make_unique<tiledb::Query>(array->schema().context(), *array);

    query_->set_layout(TILEDB_ROW_MAJOR);

    impl_ = std::make_unique<MQAux>(*this);
}

ManagedQuery::~ManagedQuery() = default;

shared_ptr<tiledb::Array> ManagedQuery::check_array(std::shared_ptr<tiledb::Array> array) {
    if (array->query_type() != TILEDB_READ) {
        throw tiledb::TileDBError("[ManagedQuery] Array not opened in read mode!");
    }

    return array;
}

void ManagedQuery::allocate_buffers() {
    auto bg = make_shared<BufferGroup>();

    auto schema = array_->schema();

    auto domain = schema.domain();
    for (const auto& dim : domain.dimensions()) {
        auto name = dim.name();

        if (impl_->skip_attrs_dims(name))
            continue;

        bg->buffers[name] = BufferSet::from_dimension(dim, initial_alloc_);
    }

    for (const auto& [name, attr] : schema.attributes()) {
        if (impl_->skip_attrs_dims(name))
            continue;

        bg->buffers[name] = BufferSet::from_attribute(attr, initial_alloc_);
    }

    buffers_ = std::move(bg);
}

void ManagedQuery::set_buffers() {
    for (auto&& [name, bg] : buffers_->buffers) {
        query_->set_data_buffer(
            name,
            (void*)bg->data.data(),
            bg->data.size() / bg->elem_nbytes()
        );

        if (bg->isvar()) {
            auto buf = bg->offsets.value();
            query_->set_offsets_buffer(
                name,
                (uint64_t*)buf.data(),
                buf.size()
            );
        }

        if (bg->isnullable()) {
            auto validity = bg->validity.value();
            query_->set_validity_buffer(
                name,
                (uint8_t*)validity.data(),
                validity.size()
            );
        }
    }
}

void ManagedQuery::validate_query() {

}

void ManagedQuery::resize_result_buffers() {
    if (query_->query_status() != tiledb::Query::Status::COMPLETE) {
        throw tiledb::TileDBError("internal error: attempted to resize result buffers but query not complete");
    }

    for (auto&& [name, sizes] : query_->result_buffer_elements_nullable()) {
        auto [offsets_nelem, data_nelem, validity_nelem] = sizes;
        auto& buf = buffers_->buffers[name];

        size_t data_nbytes = data_nelem * buf->elem_nbytes();
        buf->data.resize(data_nbytes);

        if (buf->isvar()) {
            buf->offsets.value().resize(offsets_nelem);
        }

        if (buf->isnullable()) {
            buf->offsets.value().resize(validity_nelem);
        }
    }
}

void
ManagedQuery::complete_query() {
    query_->submit();

    // TODO run to completion with realloc

    resize_result_buffers();
}

std::unique_ptr<QueryResult>
ManagedQuery::execute() {
    validate_query();

    allocate_buffers();

    set_buffers();

    complete_query();

    return make_unique<QueryResult>(std::move(buffers_));
}

}; // namespace tiledb