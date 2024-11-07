/**
 * @file   fastercsx.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2024 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * Fast construction of CSX from COO.
 */

#include <cmath>
#include <cstdint>
#include <forward_list>

#include "parallel_functions.h"
#include "span/span.hpp"

namespace tiledbsoma::fastercsx {

using namespace std::chrono_literals;

const auto MILLION = (1000 * 1000);
const auto MEBI = (1024 * 1024);

enum class Format { CSR, CSC };

typedef std::pair<uint64_t, uint64_t> Shape;

/**
 * Return [start,stop) range.
 *
 * @param n_elements Total size of range being split.
 * @param n_sections Number of sections to split the range into.
 * @param section Section to return
 * @returns A std::pair containing the half-open range for the section.
 */
template <typename T>
std::pair<T, T> get_split(T n_elements, T n_sections, T section) noexcept {
    auto base_size = n_elements / n_sections;
    auto extras = n_elements % n_sections;

    auto start = section * base_size + std::min(section, extras);
    auto stop = (section + 1) * base_size + std::min(section + 1, extras);

    return std::pair<T, T>(start, stop);
}

template <typename T>
size_t sum_over_size(const std::vector<tcb::span<const T>>& v) noexcept {
    return std::transform_reduce(
        v.cbegin(), v.cend(), 0ul, std::plus<>{}, [](tcb::span<const T> a) {
            return a.size();
        });
}

template <typename COO_INDEX>
struct Partition {
    size_t size;  // sum of views[n].size for all n
    std::vector<tcb::span<COO_INDEX>> views;
};

template <typename COO_INDEX>
void bin_view(
    std::vector<Partition<COO_INDEX>>& partitions,
    const tcb::span<COO_INDEX>& view) {
    // find minimum size partition and add the span to that partition.
    size_t min_idx = 0;
    for (size_t pidx = 1; pidx < partitions.size(); ++pidx) {
        if (partitions[pidx].size < partitions[min_idx].size)
            min_idx = pidx;
    }
    partitions[min_idx].size += view.size();
    partitions[min_idx].views.push_back(view);
}

template <typename COO_INDEX>
std::vector<Partition<COO_INDEX>> partition_views(
    std::vector<tcb::span<COO_INDEX>> const& Ai,
    const size_t max_partitions,
    const size_t partition_size) {
    // Trivial greedy k-way linear partition. If we sort (descending) first, it
    // would be LPT which would be better in edge cases where there are a lot of
    // variable sized input vectors. TODO.
    assert(max_partitions > 0);
    std::vector<Partition<COO_INDEX>> partitions(max_partitions);
    for (auto& view : Ai) {
        uint64_t n_partitions = std::min(
            (view.size() + partition_size - 1) / partition_size,
            max_partitions);
        for (uint64_t i = 0; i < n_partitions; ++i) {
            const auto [start, stop] = get_split(view.size(), n_partitions, i);
            bin_view(partitions, view.subspan(start, stop - start));
        }
    }

    // Erase any partitions that are zero length
    for (auto it = partitions.begin(); it != partitions.end();) {
        if (it->size == 0)
            it = partitions.erase(it);
        else
            ++it;
    }

    return partitions;
}

template <typename COO_INDEX, typename CSX_MAJOR_INDEX>
void count_rows(
    ThreadPool* const tp,
    uint64_t n_row,
    uint64_t nnz,
    std::vector<tcb::span<COO_INDEX const>> const& Ai,
    tcb::span<CSX_MAJOR_INDEX>& Bp) {
    assert(nnz == sum_over_size(Ai));
    assert(Bp.size() == n_row + 1);
    assert(tp->concurrency_level() >= 1);

    std::fill(Bp.begin(), Bp.end(), 0);

    auto partitions = partition_views(
        Ai, tp->concurrency_level(), 32 * MEBI /* heuristic (empirical) */);
    auto n_partitions = partitions.size();
    if (n_partitions > 1) {
        // for multiple partitions, allocate accumulators and perform in
        // parallel
        std::vector<std::vector<CSX_MAJOR_INDEX>> partition_counts(
            n_partitions, std::vector<CSX_MAJOR_INDEX>(n_row + 1, 0));

        auto status = parallel_for(
            tp,
            0ul,
            n_partitions,
            [&partition_counts, &partitions, &nnz, &n_row, &n_partitions](
                const uint64_t partition) {
                auto& counts = partition_counts[partition];
                for (auto& Ai_view : partitions[partition].views) {
                    for (size_t n = 0; n < Ai_view.size(); ++n) {
                        uint64_t row = Ai_view[n];
                        if ((row > n_row - 1) || (row < 0)) [[unlikely]]
                            throw std::out_of_range("Coordinate out of range.");
                        counts[row]++;
                    }
                }
                return Status::Ok();
            });
        assert(status.ok());

        // sum partitioned counts
        for (uint64_t partition = 0; partition < n_partitions; ++partition) {
            for (uint64_t row = 0; row < n_row + 1; ++row) {
                Bp[row] += partition_counts[partition][row];
            }
        }
    } else if (n_partitions > 0) {
        // for a single partition, just accumulate directly into the output
        // array
        assert(Ai.size() == 1);
        assert(
            Ai[0].data() == partitions[0].views[0].data() &&
            Ai[0].size() == partitions[0].views[0].size());
        auto& Ai_view = Ai[0];
        for (uint64_t n = 0; n < nnz; ++n) {
            uint64_t row = Ai_view[n];
            if ((row > n_row - 1) || (row < 0)) [[unlikely]]
                throw std::out_of_range("Coordinate out of range.");
            Bp[row]++;
        }
    }
    // else, zero length array
}

/**
 * @brief Sum counts in place to create index pointers. Same as a prefix sum,
 * but with the result stored offset by one position (ie., first position is
 * always)
 */
template <typename CSX_MAJOR_INDEX>
void sum_rows_to_pointers(tcb::span<CSX_MAJOR_INDEX>& Bp) {
    CSX_MAJOR_INDEX cumsum = 0;
    for (uint64_t i = 0; i < Bp.size(); i++) {
        auto temp = Bp[i];
        Bp[i] = cumsum;
        cumsum += temp;
    }
}

/**
 * Given COO matrix in <I,J,D> format (row, col, value), compress into CSX
 * <P,J,D>. To generate a CSC, just swap the major/minor axes.
 *
 * Caller must allocate output buffers, but need not initialize them.
 *
 * @param tp ThreadPool
 * @param n_row Size of matrix on major axis
 * @param n_col Size of matrix on minor axis
 * @param nnz Number of values
 * @param Ai Input, major axis coordinates
 * @param Aj Input, minor axis coordinates
 * @param Ad Input values
 * @param Bp Output, compressed major axis.
 * @param Bj Output, minor axis coordinates.
 * @param Bd Output values
 */
template <
    class VALUE,
    class COO_INDEX,
    class CSX_MINOR_INDEX,
    class CSX_MAJOR_INDEX>
void compress_coo(
    ThreadPool* const tp,
    const Shape& shape,
    uint64_t nnz,
    const std::vector<tcb::span<const COO_INDEX>>& Ai,
    const std::vector<tcb::span<const COO_INDEX>>& Aj,
    const std::vector<tcb::span<const VALUE>>& Ad,
    tcb::span<CSX_MAJOR_INDEX> Bp,
    tcb::span<CSX_MINOR_INDEX> Bj,
    tcb::span<VALUE> Bd) {
    auto [n_row, n_col] = shape;
    assert(Ai.size() == Aj.size() && Aj.size() == Ad.size());
    assert(sum_over_size(Ai) == nnz);
    assert(sum_over_size(Aj) == nnz);
    assert(sum_over_size(Ad) == nnz);
    assert(Bp.size() == n_row + 1);
    assert(Bj.size() == nnz);
    assert(Bd.size() == nnz);

    // get major axis counts. Important side effect: range checks the Ai values.
    count_rows(tp, n_row, nnz, Ai, Bp);
    sum_rows_to_pointers(Bp);
    assert(Bp[n_row] >= 0 && static_cast<uint64_t>(Bp[n_row]) == nnz);

#if 0  // OLD version

    // Copy all minor index and data values
    auto partition_bits = std::max(13L, std::lround(std::ceil(std::log2(n_row / tp->concurrency_level()))));
    auto n_partitions = (n_row + (1 << partition_bits) - 1) >> partition_bits;
    assert(partition_bits >= 0);
    assert((n_partitions > 0) || (n_row == 0));
    auto status =
        parallel_for(tp, 0ul, n_partitions, [&partition_bits, &Ai, &Bp, &Aj, &Bj, &Ad, &Bd](const uint64_t partition)
                     {
        for (uint64_t chnk = 0; chnk < Ai.size(); ++chnk)
        {
          auto Ai_ = Ai[chnk];
          auto Aj_ = Aj[chnk];
          auto Ad_ = Ad[chnk];
          for (uint64_t n = 0; n < Ai_.size(); ++n)
          {
            const uint64_t row = Ai_[n];
            if ((row >> partition_bits) != partition)
              continue;

            uint64_t dest = Bp[row];
            Bj[dest] = Aj_[n];
            Bd[dest] = Ad_[n];
            Bp[row]++;
          }
        }
        return Status::Ok(); });
    assert(status.ok());

    // Shift major axis pointers back by one slot
    for (uint64_t i = 0, last = 0; i <= n_row; i++)
    {
      uint64_t temp = Bp[i];
      Bp[i] = last;
      last = temp;
    }

#else  // unrolled to reduce read bandwidth by half

    std::vector<CSX_MAJOR_INDEX> Bp_left(Bp.begin(), Bp.end() - 1);
    std::vector<CSX_MAJOR_INDEX> Bp_right(Bp.begin() + 1, Bp.end());

    // Copy all minor index and data values
    const auto partition_bits = std::max(
                                    13L,
                                    std::lround(std::ceil(std::log2(
                                        n_row / tp->concurrency_level())))) +
                                1;
    const auto n_partitions = (n_row + (1 << partition_bits) - 1) >>
                              partition_bits;
    assert((n_partitions > 0) || (n_row == 0));
    auto status = parallel_for(
        tp,
        0ul,
        2 * n_partitions,
        [&n_partitions,
         &partition_bits,
         &Ai,
         &Bp_left,
         &Bp_right,
         &Aj,
         &Bj,
         &Ad,
         &Bd](const uint64_t partition) {
            for (uint64_t chnk = 0; chnk < Ai.size(); ++chnk) {
                auto Ai_ = Ai[chnk];
                auto Aj_ = Aj[chnk];
                auto Ad_ = Ad[chnk];

                const auto row_partition = partition / 2;
                const auto left_half = ((partition & 0x1) == 0x0);
                if (left_half) {
                    for (auto n = 0UL; n < (Ai_.size() / 2); ++n) {
                        const uint64_t row = Ai_[n];
                        if ((row >> partition_bits) != row_partition) {
                            continue;
                        }

                        uint64_t dest = Bp_left[row];
                        Bj[dest] = Aj_[n];
                        Bd[dest] = Ad_[n];
                        Bp_left[row]++;
                    }
                } else {
                    for (auto n = (Ai_.size() / 2); n < Ai_.size(); ++n) {
                        const uint64_t row = Ai_[n];
                        if ((row >> partition_bits) != row_partition) {
                            continue;
                        }

                        Bp_right[row]--;
                        uint64_t dest = Bp_right[row];
                        Bj[dest] = Aj_[n];
                        Bd[dest] = Ad_[n];
                    }
                }
            }

            return Status::Ok();
        });
    assert(status.ok());

#endif

    // now Bp,Bj,Bx form a CSX representation (with possible duplicates)
};

template <typename CSX_MINOR_INDEX, typename VALUE>
bool index_lt(
    const std::pair<CSX_MINOR_INDEX, VALUE>& a,
    const std::pair<CSX_MINOR_INDEX, VALUE>& b) {
    return a.first < b.first;
}

/**
 * Inplace sort of minor axis.
 */
template <class VALUE, class CSX_MINOR_INDEX, class CSX_MAJOR_INDEX>
void sort_indices(
    ThreadPool* const tp,
    uint64_t n_row,
    uint64_t nnz,
    const tcb::span<CSX_MAJOR_INDEX> Bp,
    tcb::span<CSX_MINOR_INDEX> Bj,
    tcb::span<VALUE> Bd) {
    assert(Bp.size() == n_row + 1);
    assert(Bj.size() == nnz);
    assert(Bd.size() == nnz);

    auto status = parallel_for(
        tp, 0ul, n_row, [&Bp, &Bj, &Bd, &nnz](uint64_t row) {
            uint64_t idx_start = Bp[row];
            uint64_t idx_end = Bp[row + 1];

            if (idx_end < idx_start || idx_end > nnz)
                throw std::overflow_error("Row pointer exceeds nnz");

            std::vector<std::pair<CSX_MINOR_INDEX, VALUE>> temp(
                idx_end - idx_start);
            for (uint64_t n = 0, idx = idx_start; idx < idx_end; ++n, ++idx) {
                temp[n] = std::make_pair(Bj[idx], Bd[idx]);
            }

            std::sort(
                temp.begin(), temp.end(), index_lt<CSX_MINOR_INDEX, VALUE>);

            for (uint64_t n = 0, idx = idx_start; idx < idx_end; ++n, ++idx) {
                Bj[idx] = temp[n].first;
                Bd[idx] = temp[n].second;
            }

            return Status::Ok();
        });

    assert(status.ok());
};

/**
 * @brief Unsafe copy from compressed to dense format, as a C-order contiguous
 * array
 *
 */
template <class VALUE, class CSX_MINOR_INDEX, class CSX_MAJOR_INDEX>
void copy_to_dense(
    ThreadPool* const tp,
    uint64_t major_start,
    uint64_t major_end,
    const Shape& shape,
    Format cm_format,
    const tcb::span<const CSX_MAJOR_INDEX>& Bp,
    const tcb::span<const CSX_MINOR_INDEX>& Bj,
    const tcb::span<const VALUE>& Bd,
    tcb::span<VALUE> out) {
    auto [n_row, n_col] = shape;
    assert(Bp.size() > ((cm_format == Format::CSR) ? n_row : n_col));
    assert(major_start <= ((cm_format == Format::CSR) ? n_row : n_col));
    assert(major_end <= ((cm_format == Format::CSR) ? n_row : n_col));
    if (major_start >= major_end)
        return;

    if (cm_format == Format::CSR) {
        uint64_t out_n_col = n_col;
        auto status = parallel_for(
            tp,
            major_start,
            major_end,
            [&Bp, &Bj, &Bd, &major_start, &out_n_col, &out](uint64_t i) {
                uint64_t p_start = Bp[i];
                uint64_t p_stop = Bp[i + 1];
                uint64_t out_row = i - major_start;
                for (uint64_t j = p_start; j < p_stop; ++j) {
                    uint64_t out_col = Bj[j];
                    uint64_t out_idx = out_row * out_n_col +
                                       out_col;  // C ordered
                    if (out_idx >= out.size()) [[unlikely]]
                        throw std::overflow_error(
                            "Out array is too small for provided coordinate "
                            "range.");
                    out[out_idx] = Bd[j];
                }
                return Status::Ok();
            });
        assert(status.ok());
    } else {
        // parallelizing this is less beneficial than the CSR case as it
        // partitions by columns, but writes in C order. This does not align
        // with CPU cache lines, and is therefore quite a bit slower than the
        // speedup achieved by the (above) CSR loop.
        //
        // We could write a Fortran-ordered output array, but then it would
        // typically need to be converted back to C order, which is even worse
        // as you move (even) more data once it is dense.
        //
        // If there is ever a requirement for use of F-ordered data, we could
        // add the output order as an option.
        //
        assert(cm_format == Format::CSC);
        uint64_t out_n_col = major_end - major_start;
        auto status = parallel_for(
            tp,
            major_start,
            major_end,
            [&Bp, &Bj, &Bd, &major_start, &out_n_col, &out](uint64_t i) {
                uint64_t p_start = Bp[i];
                uint64_t p_stop = Bp[i + 1];
                uint64_t out_col = i - major_start;
                for (uint64_t j = p_start; j < p_stop; ++j) {
                    uint64_t out_row = Bj[j];
                    auto out_idx = out_row * out_n_col + out_col;  // C ordered
                    if (out_idx >= out.size()) [[unlikely]]
                        throw std::overflow_error(
                            "Out array is too small for provided coordinate "
                            "range.");
                    out[out_idx] = Bd[j];
                }
                return Status::Ok();
            });
        assert(status.ok());
    }
};

}  // namespace tiledbsoma::fastercsx
