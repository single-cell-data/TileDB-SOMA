

// enum for TileDB XPtr Object type using int32_t payload (for R)
enum tiledb_xptr_object : int32_t {};
const tiledb_xptr_object tiledb_xptr_default{0};
const tiledb_xptr_object tiledb_xptr_object_array{10};
const tiledb_xptr_object tiledb_xptr_object_arrayschema{20};
const tiledb_xptr_object tiledb_xptr_object_arrayschemaevolution{30};
const tiledb_xptr_object tiledb_xptr_object_attribute{40};
const tiledb_xptr_object tiledb_xptr_object_config{50};
const tiledb_xptr_object tiledb_xptr_object_context{60};
const tiledb_xptr_object tiledb_xptr_object_dimension{70};
const tiledb_xptr_object tiledb_xptr_object_domain{80};
const tiledb_xptr_object tiledb_xptr_object_filter{90};
const tiledb_xptr_object tiledb_xptr_object_filterlist{100};
const tiledb_xptr_object tiledb_xptr_object_fragmentinfo{110};
const tiledb_xptr_object tiledb_xptr_object_group{120};
const tiledb_xptr_object tiledb_xptr_object_query{130};
const tiledb_xptr_object tiledb_xptr_object_querycondition{140};
const tiledb_xptr_object tiledb_xptr_object_vfs{150};
const tiledb_xptr_object tiledb_xptr_vfs_fh_t{160};
const tiledb_xptr_object tiledb_xptr_vlc_buf_t{170};
const tiledb_xptr_object tiledb_xptr_vlv_buf_t{180};
const tiledb_xptr_object tiledb_xptr_query_buf_t{190};

// the definitions above are internal to tiledb-r but we need a new value here
// if we want tag the external pointer
const tiledb_xptr_object tiledb_arrow_array_t{300};
const tiledb_xptr_object tiledb_arrow_schema_t{310};

const tiledb_xptr_object tiledb_soma_reader_t{500};

// templated checkers for external pointer tags
template <typename T>
const int32_t XPtrTagType = tiledb_xptr_default;  // clang++ wants a value
template <>
inline const int32_t XPtrTagType<tiledb::Array> = tiledb_xptr_object_array;
template <>
inline const int32_t
    XPtrTagType<tiledb::ArraySchema> = tiledb_xptr_object_arrayschema;
template <>
inline const int32_t XPtrTagType<tiledb::ArraySchemaEvolution> =
    tiledb_xptr_object_arrayschemaevolution;
template <>
inline const int32_t
    XPtrTagType<tiledb::Attribute> = tiledb_xptr_object_attribute;
template <>
inline const int32_t XPtrTagType<tiledb::Config> = tiledb_xptr_object_config;
template <>
inline const int32_t XPtrTagType<tiledb::Context> = tiledb_xptr_object_context;
template <>
inline const int32_t
    XPtrTagType<tiledb::Dimension> = tiledb_xptr_object_dimension;
template <>
inline const int32_t XPtrTagType<tiledb::Domain> = tiledb_xptr_object_domain;
template <>
inline const int32_t XPtrTagType<tiledb::Filter> = tiledb_xptr_object_filter;
template <>
inline const int32_t
    XPtrTagType<tiledb::FilterList> = tiledb_xptr_object_filterlist;
template <>
inline const int32_t
    XPtrTagType<tiledb::FragmentInfo> = tiledb_xptr_object_fragmentinfo;
template <>
inline const int32_t XPtrTagType<tiledb::Group> = tiledb_xptr_object_group;
template <>
inline const int32_t XPtrTagType<tiledb::Query> = tiledb_xptr_object_query;
template <>
inline const int32_t
    XPtrTagType<tiledb::QueryCondition> = tiledb_xptr_object_query;
template <>
inline const int32_t XPtrTagType<tiledb::VFS> = tiledb_xptr_object_vfs;
// this need the C API for which we do not include a header
// template <> inline const int32_t XPtrTagType<vfs_fh_t>                     =
// tiledb_xptr_vfs_fh_t; template <> inline const int32_t XPtrTagType<vlc_buf_t>
// = tiledb_xptr_vlc_buf_t; template <> inline const int32_t
// XPtrTagType<vlv_buf_t>                    = tiledb_xptr_vlv_buf_t; template
// <> inline const int32_t XPtrTagType<query_buf_t>                  =
// tiledb_xptr_query_buf_t;

template <>
inline const int32_t XPtrTagType<ArrowArray> = tiledb_arrow_array_t;
template <>
inline const int32_t XPtrTagType<ArrowSchema> = tiledb_arrow_schema_t;

template <>
inline const int32_t XPtrTagType<tdbs::SOMAArray> = tiledb_soma_reader_t;

template <typename T>
Rcpp::XPtr<T> make_xptr(T* p, bool finalize = true) {
    return Rcpp::XPtr<T>(p, finalize, Rcpp::wrap(XPtrTagType<T>), R_NilValue);
}

template <typename T>
Rcpp::XPtr<T> make_xptr(SEXP p) {
    return Rcpp::XPtr<T>(
        p);  // the default XPtr ctor with deleter on and tag and prot nil
}

template <typename T>
void check_xptr_tag(Rcpp::XPtr<T> ptr) {
    if (R_ExternalPtrTag(ptr) == R_NilValue) {
        Rcpp::stop(
            "External pointer without tag, expected tag %d\n", XPtrTagType<T>);
    }
    if (R_ExternalPtrTag(ptr) != R_NilValue) {
        int32_t tag = Rcpp::as<int32_t>(R_ExternalPtrTag(ptr));
        if (XPtrTagType<T> != tag) {
            Rcpp::stop(
                "Wrong tag type: expected %d but received %d\n",
                XPtrTagType<T>,
                tag);
        }
    }
}

// in rinterface.cpp
Rcpp::XPtr<ArrowSchema> schema_owning_xptr(void);
Rcpp::XPtr<ArrowArray> array_owning_xptr(void);
