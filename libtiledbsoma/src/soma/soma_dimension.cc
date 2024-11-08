#include "soma_dimension.h"

namespace tiledbsoma {
void SOMADimension::_set_dim_ranges(
    ManagedQuery& query, const std::any& ranges) const {
}

std::any SOMADimension::_core_domain_slot() const {
    switch (dimension.type()) {
        case TILEDB_UINT8:
            return std::make_any<std::pair<uint8_t, uint8_t>>(
                dimension.domain<uint8_t>());
        case TILEDB_UINT16:
            return std::make_any<std::pair<uint16_t, uint16_t>>(
                dimension.domain<uint16_t>());
        case TILEDB_UINT32:
            return std::make_any<std::pair<uint32_t, uint32_t>>(
                dimension.domain<uint32_t>());
        case TILEDB_UINT64:
            return std::make_any<std::pair<uint64_t, uint64_t>>(
                dimension.domain<uint64_t>());
        case TILEDB_INT8:
            return std::make_any<std::pair<int8_t, int8_t>>(
                dimension.domain<int8_t>());
        case TILEDB_INT16:
            return std::make_any<std::pair<int16_t, int16_t>>(
                dimension.domain<int16_t>());
        case TILEDB_INT32:
            return std::make_any<std::pair<int32_t, int32_t>>(
                dimension.domain<int32_t>());
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS:
        case TILEDB_INT64:
            return std::make_any<std::pair<int64_t, int64_t>>(
                dimension.domain<int64_t>());
        case TILEDB_FLOAT32:
            return std::make_any<std::pair<float_t, float_t>>(
                dimension.domain<float_t>());
        case TILEDB_FLOAT64:
            return std::make_any<std::pair<double_t, double_t>>(
                dimension.domain<double_t>());
        default:
            throw TileDBSOMAError(fmt::format(
                "[SOMADimension] Unknown dimension type {}",
                impl::type_to_str(dimension.type())));
    }
}

std::any SOMADimension::_non_empty_domain_slot() const {
    switch (dimension.type()) {
        case TILEDB_UINT8:
            return std::make_any<std::pair<uint8_t, uint8_t>>(
                array.non_empty_domain<uint8_t>(dimension.name()));
        case TILEDB_UINT16:
            return std::make_any<std::pair<uint16_t, uint16_t>>(
                array.non_empty_domain<uint16_t>(dimension.name()));
        case TILEDB_UINT32:
            return std::make_any<std::pair<uint32_t, uint32_t>>(
                array.non_empty_domain<uint32_t>(dimension.name()));
        case TILEDB_UINT64:
            return std::make_any<std::pair<uint64_t, uint64_t>>(
                array.non_empty_domain<uint64_t>(dimension.name()));
        case TILEDB_INT8:
            return std::make_any<std::pair<int8_t, int8_t>>(
                array.non_empty_domain<int8_t>(dimension.name()));
        case TILEDB_INT16:
            return std::make_any<std::pair<int16_t, int16_t>>(
                array.non_empty_domain<int16_t>(dimension.name()));
        case TILEDB_INT32:
            return std::make_any<std::pair<int32_t, int32_t>>(
                array.non_empty_domain<int32_t>(dimension.name()));
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS:
        case TILEDB_INT64:
            return std::make_any<std::pair<int64_t, int64_t>>(
                array.non_empty_domain<int64_t>(dimension.name()));
        case TILEDB_FLOAT32:
            return std::make_any<std::pair<float_t, float_t>>(
                array.non_empty_domain<float_t>(dimension.name()));
        case TILEDB_FLOAT64:
            return std::make_any<std::pair<double_t, double_t>>(
                array.non_empty_domain<double_t>(dimension.name()));
        default:
            throw TileDBSOMAError(fmt::format(
                "[SOMADimension] Unknown dimension type {}",
                impl::type_to_str(dimension.type())));
    }
}

std::any SOMADimension::_core_current_domain_slot() const {
    CurrentDomain
        current_domain = tiledb::ArraySchemaExperimental::current_domain(
            ctx, array.schema());
    NDRectangle ndrect = current_domain.ndrectangle();

    switch (dimension.type()) {
        case TILEDB_UINT8:
            std::array<uint8_t, 2> domain = ndrect.range<uint8_t>(
                dimension.name());
            return std::make_any<std::pair<uint8_t, uint8_t>>(
                std::make_pair(domain[0], domain[1]));
        case TILEDB_UINT16:
            std::array<uint16_t, 2> domain = ndrect.range<uint16_t>(
                dimension.name());
            return std::make_any<std::pair<uint16_t, uint16_t>>(
                std::make_pair(domain[0], domain[1]));
        case TILEDB_UINT32:
            std::array<uint32_t, 2> domain = ndrect.range<uint32_t>(
                dimension.name());
            return std::make_any<std::pair<uint32_t, uint32_t>>(
                std::make_pair(domain[0], domain[1]));
        case TILEDB_UINT64:
            std::array<uint64_t, 2> domain = ndrect.range<uint64_t>(
                dimension.name());
            return std::make_any<std::pair<uint64_t, uint64_t>>(
                std::make_pair(domain[0], domain[1]));
        case TILEDB_INT8:
            std::array<int8_t, 2> domain = ndrect.range<int8_t>(
                dimension.name());
            return std::make_any<std::pair<int8_t, int8_t>>(
                std::make_pair(domain[0], domain[1]));
        case TILEDB_INT16:
            std::array<int16_t, 2> domain = ndrect.range<int16_t>(
                dimension.name());
            return std::make_any<std::pair<int16_t, int16_t>>(
                std::make_pair(domain[0], domain[1]));
        case TILEDB_INT32:
            std::array<int32_t, 2> domain = ndrect.range<int32_t>(
                dimension.name());
            return std::make_any<std::pair<int32_t, int32_t>>(
                std::make_pair(domain[0], domain[1]));
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS:
        case TILEDB_INT64:
            std::array<int64_t, 2> domain = ndrect.range<int64_t>(
                dimension.name());
            return std::make_any<std::pair<int64_t, int64_t>>(
                std::make_pair(domain[0], domain[1]));
        case TILEDB_FLOAT32:
            std::array<float_t, 2> domain = ndrect.range<float_t>(
                dimension.name());
            return std::make_any<std::pair<float_t, float_t>>(
                std::make_pair(domain[0], domain[1]));
        case TILEDB_FLOAT64:
            std::array<double_t, 2> domain = ndrect.range<double_t>(
                dimension.name());
            return std::make_any<std::pair<double_t, double_t>>(
                std::make_pair(domain[0], domain[1]));
        case TILEDB_STRING_UTF8:
        case TILEDB_CHAR:
            std::array<int32_t, 2> domain = ndrect.range<int32_t>(
                dimension.name());
            return std::make_any<
                std::pair<std::vector<int32_t>, std::vector<int32_t>>>(
                std::make_pair(domain[0], domain[1]));
        default:
            throw TileDBSOMAError(fmt::format(
                "[SOMADimension] Unknown dimension type {}",
                impl::type_to_str(dimension.type())));
    }
}
}  // namespace tiledbsoma
