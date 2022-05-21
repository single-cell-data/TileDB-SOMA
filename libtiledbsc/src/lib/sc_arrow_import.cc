/* ****************************** */
/*         Arrow Importer         */
/* ****************************** */

class ArrowImporter {
   public:
    ArrowImporter(Query* const query);
    ~ArrowImporter();

    void import_(std::string name, ArrowArray* array, ArrowSchema* schema);

   private:
    Query* const query_;
    std::vector<void*> offset_buffers_;

};  // class ArrowExporter

ArrowImporter::ArrowImporter(Query* const query)
    : query_(query) {
}

ArrowImporter::~ArrowImporter() {
    for (auto p : offset_buffers_) {
        std::free(p);
    }
}

void ArrowImporter::import_(
    std::string name, ArrowArray* arw_array, ArrowSchema* arw_schema) {
    auto typeinfo = arrow_type_to_tiledb(arw_schema);

    // buffer conversion

    if (typeinfo.cell_val_num == TILEDB_VAR_NUM) {
        assert(arw_array->n_buffers == 3);

        void* p_offsets = const_cast<void*>(arw_array->buffers[1]);
        void* p_data = const_cast<void*>(arw_array->buffers[2]);
        const uint64_t num_offsets = arw_array->length;
        uint64_t data_nbytes = 0;
        if (typeinfo.arrow_large) {
            data_nbytes = static_cast<uint64_t*>(p_offsets)[num_offsets] *
                          typeinfo.elem_size;
        } else {
            data_nbytes = static_cast<uint32_t*>(p_offsets)[num_offsets] *
                          typeinfo.elem_size;
        }

        // Set the TileDB buffer, adding `1` to `num_offsets` to account for
        // the expected, extra offset.
        query_->set_data_buffer(name, p_data, data_nbytes);
        query_->set_offsets_buffer(
            name, static_cast<uint64_t*>(p_offsets), num_offsets + 1);
    } else {
        // fixed-size attribute (not TILEDB_VAR_NUM)
        assert(arw_array->n_buffers == 2);

        void* p_data = const_cast<void*>(arw_array->buffers[1]);
        uint64_t data_num = arw_array->length;

        query_->set_data_buffer(name, static_cast<void*>(p_data), data_num);
    }
}
