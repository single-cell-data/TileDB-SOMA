#include "tiledbsc_export.h"

#include <string>
#include <vector>
#include <optional>
#include <span>

#include <tiledbsc/common.h>

using namespace std;

namespace tiledbsc::util {

TILEDBSC_EXPORT void check_paths_exist(vector<string>, optional<SCConfig> config);

using VarlenBufferPair = pair< vector<DELEM_T>, vector<uint64_t> >;

template <typename T>
TILEDBSC_EXPORT VarlenBufferPair to_varlen_buffers(
    vector<T> data
);

/**
 * Convert bytemap vector to bitmap in-place. Initial size() elements will
 * be overwritten in  `size` bits.
 *
 * @param bytemap The bytemap to be converted.
 * @return Number of invalid (null) slots in the final output.
 *
 */
TILEDBSC_EXPORT size_t bytemap_to_bitmap_inplace(std::span<uint8_t> bytemap);
TILEDBSC_EXPORT size_t bytemap_to_bitmap_inplace(std::span<byte> bytemap);

} // end namespace tiledb::util