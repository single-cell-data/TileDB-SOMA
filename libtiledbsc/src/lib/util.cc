#include <bitset>
#include <cstring>
#include <iostream>

#include <tiledbsc/util.h>

using namespace tiledbsc;
using namespace tiledbsc::util;

namespace tiledbsc::util {

void check_paths_exist(vector<string>, optional<SCConfig> config) {
    // TODO
    (void)config;
};

template <typename T>
VarlenBufferPair to_varlen_buffers(vector<T> data) {
    size_t nbytes = 0;
    for (auto& elem : data) {
        nbytes += elem.size();
    }

    std::vector<DELEM_T> result(nbytes);
    std::vector<uint64_t> offsets(data.size() + 1);
    size_t offset = 0;
    size_t idx = 0;

    for (auto& elem : data) {
        std::memcpy(result.data() + offset, elem.data(), elem.size());
        offsets[idx++] = offset;

        offset += elem.size();
    }
    offsets[idx] = offset;

    return {result, offsets};
}

template VarlenBufferPair to_varlen_buffers(vector<string>);

//
// Convert a bytemap (uint8_t) to a bitmap *in-place*
// The bitmap is written to the beginning of the vector
// Vector size and padding cells are not changed.
//
size_t bytemap_to_bitmap_inplace(std::span<uint8_t> bytemap) {
    // Total number of null values
    size_t null_count = 0;

    // Number of 8 bit chunks process
    size_t nchunks = bytemap.size() / 8;

    // Number of remaining bits to process
    size_t nchunks_rem = bytemap.size() % 8;

    // Loop over chunks of 8 and write each set as bits in vector slot
    size_t idx{0};

    uint8_t bm = 0;
    for (size_t i = 0; i < nchunks; i++) {
        bm = 0;
        for (uint8_t j = 0; j < 8; j++) {
            auto val = bytemap[idx++];
            null_count += val == 0 ? 1 : 0;
            bm |= (val << j);
        }
        bytemap[i] |= bm;
    }

    bm = 0;

    if (nchunks_rem > 0) {
        // Process the last chunk
        for (size_t j = 0; j < nchunks_rem; j++) {
            auto val = bytemap[idx++];
            null_count += val == 0 ? 1 : 0;
            bm |= (val << j);
        }
        bytemap[nchunks] = bm;
    }

    return null_count;
}

size_t bytemap_to_bitmap_inplace(std::span<byte> bytemap) {
    return bytemap_to_bitmap_inplace(
        span<uint8_t>((uint8_t*)bytemap.data(), bytemap.size()));
};

};  // namespace tiledbsc::util