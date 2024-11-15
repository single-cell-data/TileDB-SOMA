#ifndef SOMA_COLUMN
#define SOMA_COLUMN

#include <algorithm>
#include <any>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include <spdlog/fmt/fmt.h>
#include "enums.h"
#include "managed_query.h"
#include "utils/arrow_adapter.h"
#include "utils/common.h"

namespace tiledbsoma {

using namespace tiledb;

class SOMAColumn {
   public:
    SOMAColumn(const Context& ctx)
        : ctx(ctx){};
    SOMAColumn(const SOMAColumn&) = default;
    SOMAColumn(SOMAColumn&&) = default;
    SOMAColumn& operator=(const SOMAColumn&) = default;
    SOMAColumn& operator=(SOMAColumn&&) = default;

    /**
     * @brief SOMADImension display name
     */
    virtual std::string name() const = 0;

    // /**
    //  * Returns the domain of the dimension.
    //  */
    // ArrowTable domain() const;

    virtual std::string domain_to_str() const = 0;

    virtual bool isIndex() const = 0;

    virtual void select_columns(
        const std::unique_ptr<ManagedQuery>& query,
        bool if_not_empty = false) const = 0;

    virtual std::optional<std::vector<Dimension>> tiledb_dimensions() = 0;

    virtual std::optional<std::vector<Attribute>> tiledb_attributes() = 0;

    virtual std::optional<std::vector<Enumeration>> tiledb_enumerations() = 0;

    virtual bool has_current_domain() = 0;

    virtual soma_column_datatype_t type() const = 0;

    virtual std::optional<tiledb_datatype_t> domain_type() const = 0;

    virtual std::optional<tiledb_datatype_t> data_type() const = 0;

    virtual ArrowArray* arrow_domain_slot(
        Array& array, enum Domainish which_kind) const = 0;

    template <typename T>
    std::pair<T, T> domain_slot(Array& array, enum Domainish which_kind) const {
        switch (which_kind) {
            case Domainish::kind_core_domain:
                return core_domain_slot<T>();
            case Domainish::kind_core_current_domain:
                return core_current_domain_slot<T>(array);
            case Domainish::kind_non_empty_domain:
                return non_empty_domain_slot<T>(array);
            default:
                throw std::runtime_error(
                    "internal coding error in SOMAArray::_core_domainish_slot: "
                    "unknown kind");
        }
    }

    void set_current_domain_slot(
        NDRectangle& rectangle, const std::vector<const void*>& domain) const {
        if (!isIndex()) {
            throw TileDBSOMAError(fmt::format(
                "[SOMAColumn] Column with name {} is not an index column",
                name()));
        }

        _set_current_domain_slot(rectangle, domain);
    }

    template <typename T>
    void set_dim_point(
        const std::unique_ptr<ManagedQuery>& query, const T& point) const {
        if (!isIndex()) {
            throw TileDBSOMAError(fmt::format(
                "[SOMAColumn] Column with name {} is not an index column",
                name()));
        }

        this->_set_dim_points(query, std::vector({point}));
    }

    template <typename T>
    void set_dim_points(
        const std::unique_ptr<ManagedQuery>& query,
        const std::vector<T>& points) const {
        if (!isIndex()) {
            throw TileDBSOMAError(fmt::format(
                "[SOMAColumn] Column with name {} is not an index column",
                name()));
        }

        this->_set_dim_points(query, points);
    }

    template <typename T>
    void set_dim_ranges(
        const std::unique_ptr<ManagedQuery>& query,
        const std::vector<std::pair<T, T>>& ranges) const {
        if (!isIndex()) {
            throw TileDBSOMAError(fmt::format(
                "[SOMAColumn] Column with name {} is not an index column",
                name()));
        }

        this->_set_dim_ranges(query, ranges);
    }

    template <typename T>
    std::pair<T, T> core_domain_slot() const {
        try {
            return std::any_cast<std::pair<T, T>>(_core_domain_slot());
        } catch (const std::exception& e) {
            throw TileDBSOMAError(e.what());
        }
    }

    template <typename T>
    std::pair<T, T> non_empty_domain_slot(Array& array) const {
        try {
            return std::any_cast<std::pair<T, T>>(
                _non_empty_domain_slot(array));
        } catch (const std::exception& e) {
            throw TileDBSOMAError(e.what());
        };
    }

    template <typename T>
    std::pair<T, T> core_current_domain_slot(Array& array) const {
        try {
            return std::any_cast<std::pair<T, T>>(
                _core_current_domain_slot(array));
        } catch (const std::exception& e) {
            throw TileDBSOMAError(e.what());
        }
    }

   protected:
    virtual void _set_dim_points(
        const std::unique_ptr<ManagedQuery>& query,
        const std::any& points) const = 0;

    virtual void _set_dim_ranges(
        const std::unique_ptr<ManagedQuery>& query,
        const std::any& ranges) const = 0;

    virtual void _set_current_domain_slot(
        NDRectangle& rectangle,
        const std::vector<const void*>& domain) const = 0;

    virtual std::any _core_domain_slot() const = 0;

    virtual std::any _non_empty_domain_slot(Array& array) const = 0;

    virtual std::any _core_current_domain_slot(Array& array) const = 0;

    const Context& ctx;
};

template <>
std::pair<std::string, std::string> SOMAColumn::core_domain_slot<std::string>()
    const;

template <>
std::pair<std::string, std::string>
SOMAColumn::core_current_domain_slot<std::string>(Array& array) const;

}  // namespace tiledbsoma
#endif