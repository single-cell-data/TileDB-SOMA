/**
 * @file   fastercsx.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * Fast construction of CSX from COO.
 */

#include <atomic>
#include <cmath>
#include <cstdint>
#include <format>
#include <numeric>
#include <span>

#include "parallel_functions.h"

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
std::pair<T, T> get_split_(T n_elements, T n_sections, T section) noexcept {
    auto base_size = n_elements / n_sections;
    auto extras = n_elements % n_sections;

    auto start = section * base_size + std::min(section, extras);
    auto stop = (section + 1) * base_size + std::min(section + 1, extras);

    return std::pair<T, T>(start, stop);
}

#ifndef NDEBUG
/**
 * @brief Sum size over all elements in vector - used only in asserts.
 */
template <typename T>
size_t sum_over_size_(const std::vector<std::span<const T>>& v) noexcept {
    return std::transform_reduce(
        v.cbegin(), v.cend(), 0ul, std::plus<>{}, [](std::span<const T> a) {
            return a.size();
        });
}

/**
 * @brief Return true if any of the tuple sub-arrays are not equal in size,
 * false otherwise.  Used only in asserts.
 */
template <typename Ti, typename Tj, typename Td>
bool no_ragged_chunks(
    const std::vector<std::span<const Ti>>& Ai,
    const std::vector<std::span<const Tj>>& Aj,
    const std::vector<std::span<const Td>>& Ad) {
    if ((Ai.size() != Aj.size()) || (Ai.size() != Ad.size()))
        return false;

    for (uint64_t chunk_idx = 0; chunk_idx < Ai.size(); chunk_idx++) {
        if ((Ai[chunk_idx].size() != Aj[chunk_idx].size()) ||
            (Ai[chunk_idx].size() != Ad[chunk_idx].size()))
            return false;
    }
    return true;
}
#endif

/**
 * @brief Simple partitioning scheme - given a vector of spans, and a minimum
 * partition size, create up to max_partitions.
 *
 * Used to split work across workers in a reasonably even manner. Currently a
 * greedy n-way linear partition, which works well enough for the common case
 * where most inputs are of similar sizes (which is typically the case).
 *
 * TODO: If we sort (descending) first, it will be classic LPT (Largest
 * Processing Time) partitioning which would be better in edge cases where there
 * are a lot of variable sized input vectors.
 */
template <typename COO_IDX>
struct Partition {
    size_t size;  // sum of views[n].size for all n
    std::vector<std::span<COO_IDX>> views;
};

template <typename COO_IDX>
void bin_view_(
    std::vector<Partition<COO_IDX>>& partitions,
    const std::span<COO_IDX>& view) {
    // find minimum size partition and add the span to that partition.
    size_t min_idx = 0;
    for (size_t pidx = 1; pidx < partitions.size(); ++pidx) {
        if (partitions[pidx].size < partitions[min_idx].size)
            min_idx = pidx;
    }
    partitions[min_idx].size += view.size();
    partitions[min_idx].views.push_back(view);
}

template <typename COO_IDX>
std::vector<Partition<COO_IDX>> partition_views_(
    std::vector<std::span<COO_IDX>> const& Ai,
    const size_t max_partitions,
    const size_t partition_size) {
    assert(max_partitions > 0);
    std::vector<Partition<COO_IDX>> partitions(max_partitions);
    for (auto& view : Ai) {
        size_t n_partitions = std::min(
            (view.size() + partition_size - 1) / partition_size,
            max_partitions);
        for (size_t i = 0; i < n_partitions; ++i) {
            const auto [start, stop] = get_split_(view.size(), n_partitions, i);
            bin_view_(partitions, view.subspan(start, stop - start));
        }
    }

    // Erase any partitions that are zero length
    partitions.erase(
        std::remove_if(
            partitions.begin(),
            partitions.end(),
            [](const Partition<COO_IDX>& p) { return p.size == 0; }),
        partitions.end());

    return partitions;
}

/**
 * @brief Count occurrences of all values present in Ai -- used as the basis
 * for constructing the index pointer array.
 */
template <typename COO_IDX, typename CSX_MAJOR_IDX>
void count_rows_(
    ThreadPool* const tp,
    uint64_t n_row,
    uint64_t nnz,
    std::vector<std::span<COO_IDX const>> const& Ai,
    std::span<CSX_MAJOR_IDX>& Bp) {
    assert(nnz == sum_over_size_(Ai));
    assert(Bp.size() == n_row + 1);
    assert(tp->concurrency_level() >= 1);

    std::fill(Bp.begin(), Bp.end(), 0);

    auto partitions = partition_views_(
        Ai, tp->concurrency_level(), 32 * MEBI /* heuristic (empirical) */);
    auto n_partitions = partitions.size();
    if (n_partitions > 1) {
        // for multiple partitions, allocate accumulators and perform in
        // parallel
        std::vector<std::vector<CSX_MAJOR_IDX>> partition_counts(
            n_partitions, std::vector<CSX_MAJOR_IDX>(n_row + 1, 0));

        auto status = parallel_for(
            tp,
            0ul,
            n_partitions,
            [&partition_counts, &partitions, &n_row](const uint64_t partition) {
                auto& counts = partition_counts[partition];
                for (auto& Ai_view : partitions[partition].views) {
                    for (size_t n = 0; n < Ai_view.size(); n++) {
                        auto row = Ai_view[n];
                        if ((row < 0) ||
                            (static_cast<std::make_unsigned_t<COO_IDX>>(row) >=
                             n_row)) [[unlikely]] {
                            throw std::out_of_range(std::format(
                                "First coordinate {} out of range {}.",
                                row,
                                0,
                                n_row));
                        }
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
        for (auto& Ai_view : partitions[0].views) {
            for (size_t n = 0; n < Ai_view.size(); n++) {
                auto row = Ai_view[n];
                if ((row < 0) ||
                    (static_cast<std::make_unsigned_t<COO_IDX>>(row) >= n_row))
                    [[unlikely]] {
                    throw std::out_of_range(std::format(
                        "First coordinate {} out of range {}.", row, 0, n_row));
                }
                Bp[row]++;
            }
        }
    }
    // else, zero length array
}

/**
 * @brief Sum counts in place to create index pointers. Same as a prefix sum,
 * but with the result stored offset by one position (ie., first position is
 * always zero).
 */
template <typename CSX_MAJOR_IDX>
void sum_rows_to_pointers_(std::span<CSX_MAJOR_IDX>& Bp) {
    CSX_MAJOR_IDX cumsum = 0;
    for (uint64_t i = 0; i < Bp.size(); i++) {
        auto temp = Bp[i];
        Bp[i] = cumsum;
        cumsum += temp;
    }
}

template <
    typename VALUE,
    typename COO_IDX,
    typename CSX_MINOR_IDX,
    typename CSX_MAJOR_IDX>
void compress_coo_inner_left_(
    const uint64_t& row_partition,
    const int& partition_bits,
    const uint64_t& n_col,
    std::span<COO_IDX const>& Ai_,
    std::span<COO_IDX const>& Aj_,
    std::span<VALUE const>& Ad_,
    std::span<CSX_MAJOR_IDX>& Bp,
    std::span<CSX_MINOR_IDX>& Bj,
    std::span<VALUE>& Bd) {
    for (auto n = 0UL; n < (Ai_.size() / 2); ++n) {
        const std::make_unsigned_t<COO_IDX> row = Ai_[n];
        if ((row >> partition_bits) != row_partition)
            continue;

        const auto dest = Bp[row];
        if ((Aj_[n] < 0) ||
            (static_cast<std::make_unsigned_t<COO_IDX>>(Aj_[n]) >= n_col))
            [[unlikely]] {
            throw std::out_of_range(std::format(
                "Second coordinate {} out of range {}.", Aj_[n], 0, n_col));
        }
        Bj[dest] = Aj_[n];
        Bd[dest] = Ad_[n];
        Bp[row]++;
    }
}

template <
    typename VALUE,
    typename COO_IDX,
    typename CSX_MINOR_IDX,
    typename CSX_MAJOR_IDX>
void compress_coo_inner_right_(
    unsigned int row_partition,
    unsigned int partition_bits,
    uint64_t n_col,
    std::span<COO_IDX const>& Ai_,
    std::span<COO_IDX const>& Aj_,
    std::span<VALUE const>& Ad_,
    std::span<CSX_MAJOR_IDX>& Bp,
    std::span<CSX_MINOR_IDX>& Bj,
    std::span<VALUE>& Bd) {
    for (auto n = (Ai_.size() / 2); n < Ai_.size(); ++n) {
        const std::make_unsigned_t<COO_IDX> row = Ai_[n];
        if ((row >> partition_bits) != row_partition) {
            continue;
        }

        Bp[row]--;
        const auto dest = Bp[row];
        if ((Aj_[n] < 0) ||
            (static_cast<std::make_unsigned_t<COO_IDX>>(Aj_[n]) >= n_col))
            [[unlikely]] {
            throw std::out_of_range(std::format(
                "Second coordinate {} out of range {}.", Aj_[n], 0, n_col));
        }

        Bj[dest] = Aj_[n];
        Bd[dest] = Ad_[n];
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
template <class VALUE, class COO_IDX, class CSX_MINOR_IDX, class CSX_MAJOR_IDX>
void compress_coo(
    ThreadPool* const tp,
    const Shape& shape,
    uint64_t nnz,
    const std::vector<std::span<const COO_IDX>>& Ai,
    const std::vector<std::span<const COO_IDX>>& Aj,
    const std::vector<std::span<const VALUE>>& Ad,
    std::span<CSX_MAJOR_IDX> Bp,
    std::span<CSX_MINOR_IDX> Bj,
    std::span<VALUE> Bd) {
    auto n_row = shape.first;
    auto n_col = shape.second;
    assert(Ai.size() == Aj.size() && Aj.size() == Ad.size());
    assert(sum_over_size_(Ai) == nnz);
    assert(sum_over_size_(Aj) == nnz);
    assert(sum_over_size_(Ad) == nnz);
    assert(no_ragged_chunks(Ai, Aj, Ad));
    assert(Bp.size() == n_row + 1);
    assert(Bj.size() == nnz);
    assert(Bd.size() == nnz);

    // Get major axis counts. Important side effect: range checks the Ai
    // values.
    count_rows_(tp, n_row, nnz, Ai, Bp);
    sum_rows_to_pointers_(Bp);
    assert(Bp[n_row] >= 0 && static_cast<uint64_t>(Bp[n_row]) == nnz);

    std::vector<CSX_MAJOR_IDX> Bp_left(Bp.begin(), Bp.end() - 1);
    std::vector<CSX_MAJOR_IDX> Bp_right(Bp.begin() + 1, Bp.end());
    std::span<CSX_MAJOR_IDX> Bp_left_span{Bp_left};
    std::span<CSX_MAJOR_IDX> Bp_right_span{Bp_right};

    // Parallel, lock-free copy of all minor index and data values. Partitioned
    // by contiguous major dimension (Ai) values, minimizing the write
    // contention in the indptr/Bp updates. Each partition is split into two
    // parts, which can run in parallel, by working either left or right (minor
    // dimension beginning and end of range) - this reduces the number of major
    // dimension reads by half.
    //
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
        [&partition_bits,
         &Ai,
         &Bp_left_span,
         &Bp_right_span,
         &Aj,
         &Bj,
         &Ad,
         &Bd,
         &n_col](const uint64_t partition) {
            for (uint64_t chnk = 0; chnk < Ai.size(); ++chnk) {
                auto Ai_ = Ai[chnk];
                auto Aj_ = Aj[chnk];
                auto Ad_ = Ad[chnk];

                const auto row_partition = partition / 2;
                const auto left_half = ((partition & 0x1) == 0x0);
                if (left_half) {
                    compress_coo_inner_left_(
                        row_partition,
                        partition_bits,
                        n_col,
                        Ai_,
                        Aj_,
                        Ad_,
                        Bp_left_span,
                        Bj,
                        Bd);
                } else {
                    compress_coo_inner_right_(
                        row_partition,
                        partition_bits,
                        n_col,
                        Ai_,
                        Aj_,
                        Ad_,
                        Bp_right_span,
                        Bj,
                        Bd);
                }
            }

            return Status::Ok();
        });
    assert(status.ok());
    // now Bp,Bj,Bx form a CSX representation (with possible duplicates)
};

template <typename CSX_MINOR_IDX, typename VALUE>
bool index_lt_(
    const std::pair<CSX_MINOR_IDX, VALUE>& a,
    const std::pair<CSX_MINOR_IDX, VALUE>& b) {
    return a.first < b.first;
}

/**
 * @brief In-place sort of minor axis, used to canonicalize the CSx ordering.
 *
 * Returns false if there are duplicate coordinates, true if all coordinates
 * are unique.
 */
template <class VALUE, class CSX_MINOR_IDX, class CSX_MAJOR_IDX>
bool sort_csx_indices(
    ThreadPool* const tp,
    uint64_t n_row,
    uint64_t nnz,
    const std::span<CSX_MAJOR_IDX> Bp,
    std::span<CSX_MINOR_IDX> Bj,
    std::span<VALUE> Bd) {
    assert(Bp.size() == n_row + 1);
    assert(Bj.size() == nnz);
    assert(Bd.size() == nnz);

    std::atomic<bool> no_duplicates(true);

    auto status = parallel_for(
        tp, 0ul, n_row, [&Bp, &Bj, &Bd, &nnz, &no_duplicates](uint64_t row) {
            uint64_t idx_start = Bp[row];
            uint64_t idx_end = Bp[row + 1];

            if (idx_end < idx_start || idx_end > nnz)
                throw std::overflow_error("Row pointer exceeds nnz");

            std::vector<std::pair<CSX_MINOR_IDX, VALUE>> temp(
                idx_end - idx_start);
            for (uint64_t n = 0, idx = idx_start; idx < idx_end; ++n, ++idx) {
                temp[n] = std::make_pair(Bj[idx], Bd[idx]);
            }

            std::sort(
                temp.begin(), temp.end(), index_lt_<CSX_MINOR_IDX, VALUE>);

            for (uint64_t n = 0, idx = idx_start; idx < idx_end; ++n, ++idx) {
                Bj[idx] = temp[n].first;
                Bd[idx] = temp[n].second;

                if (n > 0 && Bj[idx] == Bj[idx - 1])
                    no_duplicates = false;
            }

            return Status::Ok();
        });

    assert(status.ok());
    return no_duplicates;
};

/**
 * @brief Unsafe copy from compressed to dense format, as a C-order
 * contiguous array
 *
 */
template <class VALUE, class CSX_MINOR_IDX, class CSX_MAJOR_IDX>
void copy_csx_to_dense(
    ThreadPool* const tp,
    uint64_t major_start,
    uint64_t major_end,
    const Shape& shape,
    Format cm_format,
    const std::span<const CSX_MAJOR_IDX>& Bp,
    const std::span<const CSX_MINOR_IDX>& Bj,
    const std::span<const VALUE>& Bd,
    std::span<VALUE> out) {
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
                            "Out array is too small for provided "
                            "coordinate "
                            "range.");
                    out[out_idx] = Bd[j];
                }
                return Status::Ok();
            });
        assert(status.ok());
    } else {
        // Parallelizing this is less beneficial than the CSR case as it
        // partitions by columns, but writes in C order. This does not align
        // with CPU cache lines, and is therefore quite a bit slower than
        // the speedup achieved by the (above) CSR loop.
        //
        // We could write a Fortran-ordered ("column major") output array, but
        // then it would typically need to be converted back to C order ("row
        // major"), which is even worse as you move (even) more data once it is
        // dense.
        //
        // If there is ever a requirement for use of F-ordered data, we
        // could add the output order as an option.
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
                            "Out array is too small for provided "
                            "coordinate "
                            "range.");
                    out[out_idx] = Bd[j];
                }
                return Status::Ok();
            });
        assert(status.ok());
    }
};

}  // namespace tiledbsoma::fastercsx
