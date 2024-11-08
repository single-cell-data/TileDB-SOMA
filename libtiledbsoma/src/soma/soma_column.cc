#include "soma_column.h"

namespace tiledbsoma {
template <>
std::pair<std::string, std::string> SOMAColumn::core_domain_slot() const {
    return std::pair<std::string, std::string>("", "");
}

template <>
std::pair<std::string, std::string> SOMAColumn::core_current_domain_slot()
    const {
    try {
        std::pair<std::string, std::string>
            current_domain = std::any_cast<std::pair<std::string, std::string>>(
                _core_current_domain_slot());

        if (current_domain.first == "" && current_domain.second == "\xff") {
            return std::pair<std::string, std::string>("", "");
        }

        return current_domain;
    } catch (const std::exception& e) {
        throw TileDBSOMAError(e.what());
    }
}
}  // namespace tiledbsoma
