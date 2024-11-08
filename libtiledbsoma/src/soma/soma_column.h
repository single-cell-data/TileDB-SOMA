#ifndef SOMA_COLUMN
#define SOMA_COLUMN

#include <algorithm>
#include <any>
#include <string>
#include <tuple>
#include <vector>

#include <spdlog/fmt/fmt.h>
#include "managed_query.h"
#include "utils/common.h"

namespace tiledbsoma {

class SOMAColumn {
   public:
    SOMAColumn(const Context& ctx, Array& array)
        : ctx(ctx)
        , array(array){};
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

    template <typename T>
    void set_dim_ranges(
        ManagedQuery& query, const std::vector<std::pair<T, T>>& ranges) const {
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
    std::pair<T, T> non_empty_domain_slot() const {
        try {
            return std::any_cast<std::pair<T, T>>(_non_empty_domain_slot());
        } catch (const std::exception& e) {
            throw TileDBSOMAError(e.what());
        };
    }

    template <typename T>
    std::pair<T, T> core_current_domain_slot() const {
        try {
            return std::any_cast<std::pair<T, T>>(_core_current_domain_slot());
        } catch (const std::exception& e) {
            throw TileDBSOMAError(e.what());
        }
    }

   protected:
    virtual void _set_dim_ranges(
        ManagedQuery& query, const std::any& ranges) const = 0;

    virtual std::any _core_domain_slot() const = 0;

    virtual std::any _non_empty_domain_slot() const = 0;

    virtual std::any _core_current_domain_slot() const = 0;

    const Context& ctx;
    Array& array;
};

}  // namespace tiledbsoma
#endif