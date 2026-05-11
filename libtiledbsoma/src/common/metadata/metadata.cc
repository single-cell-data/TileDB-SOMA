#include "metadata.h"
#include "../datatype/utils.h"
#include "utils.h"

#include <tiledb/tiledb>

namespace tiledbsoma::common {

namespace impl {
template <typename T, typename... Ts>
struct is_vector : std::false_type {};

template <typename T, typename... Ts>
struct is_vector<std::vector<T, Ts...>> : std::true_type {};
}  // namespace impl

MetadataCache::MetadataCache(tiledb::Array& array) {
    for (uint64_t idx = 0; idx < array.metadata_num(); ++idx) {
        std::string key;
        tiledb_datatype_t value_type;
        uint32_t value_num;
        const void* value;
        array.get_metadata_from_index(idx, &key, &value_type, &value_num, &value);

        metadata_[key] = decode_metadata(type::as<type::DataTypeFormat::SOMA>(value_type), value_num, value);
    }
}

MetadataCache::MetadataCache(tiledb::Group& group) {
    for (uint64_t idx = 0; idx < group.metadata_num(); ++idx) {
        std::string key;
        tiledb_datatype_t value_type;
        uint32_t value_num;
        const void* value;
        group.get_metadata_from_index(idx, &key, &value_type, &value_num, &value);

        metadata_[key] = decode_metadata(type::as<type::DataTypeFormat::SOMA>(value_type), value_num, value);
    }
}

bool MetadataCache::contains(const std::string& key) const {
    return metadata_.contains(key);
}

std::map<std::string, MetadataValue> MetadataCache::get() const {
    return metadata_;
}

std::optional<MetadataValue> MetadataCache::get(const std::string& key) const {
    if (!metadata_.contains(key)) {
        return std::nullopt;
    }

    return metadata_.at(key);
}

void MetadataCache::set(const std::string& key, MetadataValue value) {
    // Sanitize key
    for (size_t i = 0; i < key.size(); ++i) {
        if (key[i] == 0) {
            throw std::runtime_error("string contains NULL bytes");
        }
    }

    // Sanitize value
    std::visit(
        [](auto&& arg) {
            using T = std::decay_t<decltype(arg)>;

            if constexpr (std::is_same_v<T, std::string>) {
                for (size_t i = 0; i < arg.size(); ++i) {
                    if (arg[i] == 0) {
                        throw std::runtime_error("string contains NULL bytes");
                    }
                }
            }
        },
        value);

    auto current_state = current_state_(key);
    metadata_[key] = value;
    mods_[key] = MetadataCache::next_state_(current_state, "set");
}

void MetadataCache::del(const std::string& key) {
    auto current_state = current_state_(key);
    metadata_.erase(key);
    mods_[key] = MetadataCache::next_state_(current_state, "del");
}

void MetadataCache::write(tiledb::Array& array) {
    if (mods_.empty()) {
        return;
    }

    for (const auto& [key, state] : mods_) {
        if (state == DictMod::added || state == DictMod::updated) {
            std::visit(
                [&](auto&& arg) {
                    using T = std::decay_t<decltype(arg)>;

                    if constexpr (impl::is_vector<T>::value) {
                        auto datatype = tiledb::impl::type_to_tiledb<typename T::value_type>::tiledb_type;

                        if constexpr (std::is_same_v<T, std::vector<bool>>) {
                            std::vector<uint8_t> bytemap(arg.begin(), arg.end());
                            array.put_metadata(key, datatype, bytemap.size(), bytemap.data());
                        } else {
                            array.put_metadata(key, datatype, arg.size(), arg.data());
                        }
                    } else if constexpr (std::is_same_v<T, std::string>) {
                        array.put_metadata(key, TILEDB_STRING_UTF8, arg.size(), arg.size() ? arg.data() : nullptr);
                    } else {
                        auto datatype = tiledb::impl::type_to_tiledb<T>::tiledb_type;

                        if constexpr (std::is_same_v<T, bool>) {
                            std::uint8_t byte = static_cast<uint8_t>(arg);
                            array.put_metadata(key, datatype, 1, &byte);
                        } else {
                            array.put_metadata(key, datatype, 1, &arg);
                        }
                    }
                },
                metadata_[key]);
        } else if (state == DictMod::deleted) {
            array.delete_metadata(key);
        }
    }

    mods_.clear();
}

void MetadataCache::write(tiledb::Group& group) {
    if (mods_.empty()) {
        return;
    }

    for (const auto& [key, state] : mods_) {
        if (state == DictMod::added || state == DictMod::updated) {
            std::visit(
                [&](auto&& arg) {
                    using T = std::decay_t<decltype(arg)>;

                    if constexpr (impl::is_vector<T>::value) {
                        auto datatype = tiledb::impl::type_to_tiledb<typename T::value_type>::tiledb_type;

                        if constexpr (std::is_same_v<T, std::vector<bool>>) {
                            std::vector<uint8_t> bytemap(arg.begin(), arg.end());
                            group.put_metadata(key, datatype, bytemap.size(), bytemap.data());
                        } else {
                            group.put_metadata(key, datatype, arg.size(), arg.data());
                        }
                    } else if constexpr (std::is_same_v<T, std::string>) {
                        group.put_metadata(key, TILEDB_STRING_UTF8, arg.size(), arg.size() ? arg.data() : nullptr);
                    } else {
                        auto datatype = tiledb::impl::type_to_tiledb<T>::tiledb_type;

                        if constexpr (std::is_same_v<T, bool>) {
                            std::uint8_t byte = static_cast<uint8_t>(arg);
                            group.put_metadata(key, datatype, 1, &byte);
                        } else {
                            group.put_metadata(key, datatype, 1, &arg);
                        }
                    }
                },
                metadata_[key]);
        } else if (state == DictMod::deleted) {
            group.delete_metadata(key);
        }
    }

    mods_.clear();
}

std::size_t MetadataCache::size() const {
    return metadata_.size();
}

MetadataCache::DictMod MetadataCache::current_state_(const std::string& key) const {
    if (mods_.contains(key)) {
        return mods_.at(key);
    }

    return metadata_.contains(key) ? MetadataCache::DictMod::present : MetadataCache::DictMod::absent;
}

MetadataCache::DictMod MetadataCache::next_state_(DictMod current_state, const std::string& action) {
    const std::map<std::pair<MetadataCache::DictMod, std::string>, MetadataCache::DictMod> action_map = {
        {{MetadataCache::DictMod::absent, "set"}, MetadataCache::DictMod::added},
        {{MetadataCache::DictMod::added, "set"}, MetadataCache::DictMod::added},
        {{MetadataCache::DictMod::added, "del"}, MetadataCache::DictMod::deleted},
        {{MetadataCache::DictMod::present, "set"}, MetadataCache::DictMod::updated},
        {{MetadataCache::DictMod::present, "del"}, MetadataCache::DictMod::deleted},
        {{MetadataCache::DictMod::updated, "set"}, MetadataCache::DictMod::updated},
        {{MetadataCache::DictMod::updated, "del"}, MetadataCache::DictMod::deleted},
        {{MetadataCache::DictMod::deleted, "set"}, MetadataCache::DictMod::updated}};

    return action_map.at({current_state, action});
}
}  // namespace tiledbsoma::common