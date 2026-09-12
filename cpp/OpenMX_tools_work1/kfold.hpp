#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace spadexp::kfold {

using Vector3 = std::array<double, 3>;
using Shift3 = std::array<int, 3>;
using ReciprocalCell = std::array<Vector3, 3>;

struct FoldResult {
  Vector3 reduced;
  Shift3 shift;
};

// Fold a fractional reciprocal coordinate into the Wigner-Seitz first BZ.
// requested = reduced + shift, where shift gives integer G1/G2/G3 counts.
FoldResult fold_to_first_bz(const Vector3& requested,
                            const ReciprocalCell& reciprocal_cell);

// Expand blocks calculated at unique first-BZ k points back to the requested
// repeated-zone order. unique_index[requested] selects the source block.
std::vector<double> expand_k_blocks(
    const double* source,
    std::size_t unique_count,
    std::size_t block_size,
    const std::vector<int>& unique_index);

// Ownership-taking adapter used by the OpenMX postprocessor.
double* expand_k_blocks_owned(
    double* source,
    std::size_t unique_count,
    std::size_t block_size,
    const std::vector<int>& unique_index);

}  // namespace spadexp::kfold
