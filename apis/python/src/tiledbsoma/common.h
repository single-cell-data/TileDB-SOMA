#include <exception>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledb/tiledb>  // C++
#include <tiledbsoma/tiledbsoma>

using namespace std;
using namespace tiledb;
namespace py = pybind11;

#define TPY_ERROR_LOC(m) throw TileDBSOMAError(m);

namespace tiledbsoma {

py::dtype tdb_to_np_dtype(tiledb_datatype_t type, uint32_t cell_val_num);

tiledb_datatype_t np_to_tdb_dtype(py::dtype type);

bool is_tdb_str(tiledb_datatype_t type);

std::optional<py::object> to_table(
    std::optional<std::shared_ptr<ArrayBuffers>> buffers);

py::dict meta(std::map<std::string, MetadataValue> metadata_mapping);
void set_metadata(
    SOMAObject& soma_object, const std::string& key, py::array value);

class PyQueryCondition {
   private:
    Context ctx_;
    shared_ptr<QueryCondition> qc_;

   public:
    PyQueryCondition() {
        try {
            // create one global context for all query conditions
            static Context context = Context();
            ctx_ = context;
            qc_ = shared_ptr<QueryCondition>(new QueryCondition(ctx_));
        } catch (TileDBError& e) {
            TPY_ERROR_LOC(e.what());
        }
    }

    PyQueryCondition(py::object ctx) {
        (void)ctx;
        try {
            // create one global context for all query conditions
            static Context context = Context();
            ctx_ = context;
            qc_ = shared_ptr<QueryCondition>(new QueryCondition(ctx_));
        } catch (TileDBError& e) {
            TPY_ERROR_LOC(e.what());
        }
    }

    void init(
        const string& attribute_name,
        const string& condition_value,
        tiledb_query_condition_op_t op) {
        try {
            qc_->init(attribute_name, condition_value, op);
        } catch (TileDBError& e) {
            TPY_ERROR_LOC(e.what());
        }
    }

    template <typename T>
    void init(
        const string& attribute_name,
        T condition_value,
        tiledb_query_condition_op_t op) {
        try {
            qc_->init(
                attribute_name, &condition_value, sizeof(condition_value), op);
        } catch (TileDBError& e) {
            TPY_ERROR_LOC(e.what());
        }
    }

    shared_ptr<QueryCondition> ptr() {
        return qc_;
    }

    py::capsule __capsule__() {
        return py::capsule(&qc_, "qc");
    }

    template <typename T>
    static PyQueryCondition create(
        const std::string& field_name,
        const std::vector<T>& values,
        tiledb_query_condition_op_t op) {
        auto pyqc = PyQueryCondition();

        const Context ctx = std::as_const(pyqc.ctx_);

        auto set_membership_qc = QueryConditionExperimental::create(
            ctx, field_name, values, op);

        pyqc.qc_ = std::make_shared<QueryCondition>(
            std::move(set_membership_qc));

        return pyqc;
    }

    PyQueryCondition combine(
        PyQueryCondition qc,
        tiledb_query_condition_combination_op_t combination_op) const {
        auto pyqc = PyQueryCondition(nullptr, ctx_.ptr().get());

        tiledb_query_condition_t* combined_qc = nullptr;
        ctx_.handle_error(
            tiledb_query_condition_alloc(ctx_.ptr().get(), &combined_qc));

        ctx_.handle_error(tiledb_query_condition_combine(
            ctx_.ptr().get(),
            qc_->ptr().get(),
            qc.qc_->ptr().get(),
            combination_op,
            &combined_qc));

        pyqc.qc_ = std::shared_ptr<QueryCondition>(
            new QueryCondition(pyqc.ctx_, combined_qc));

        return pyqc;
    }

   private:
    PyQueryCondition(shared_ptr<QueryCondition> qc, tiledb_ctx_t* c_ctx)
        : qc_(qc) {
        ctx_ = Context(c_ctx, false);
    }

    void set_ctx(py::object ctx) {
        tiledb_ctx_t* c_ctx;
        if ((c_ctx = (py::capsule)ctx.attr("__capsule__")()) == nullptr)
            TPY_ERROR_LOC("Invalid context pointer!")

        ctx_ = Context(c_ctx, false);
    }
};
}  // namespace tiledbsoma
