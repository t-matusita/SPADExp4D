#include "kfold.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace spadexp::kfold {

FoldResult fold_to_first_bz(const Vector3& requested,
                            const ReciprocalCell& reciprocal_cell) {
  double best_norm2 = std::numeric_limits<double>::infinity();
  const Shift3 center = {
      static_cast<int>(std::llround(requested[0])),
      static_cast<int>(std::llround(requested[1])),
      static_cast<int>(std::llround(requested[2])),
  };
  FoldResult best{};
  for (int n0 = center[0] - 3; n0 <= center[0] + 3; ++n0) {
    for (int n1 = center[1] - 3; n1 <= center[1] + 3; ++n1) {
      for (int n2 = center[2] - 3; n2 <= center[2] + 3; ++n2) {
        const Shift3 shift = {n0, n1, n2};
        const Vector3 reduced = {
            requested[0] - n0,
            requested[1] - n1,
            requested[2] - n2,
        };
        Vector3 cartesian = {0.0, 0.0, 0.0};
        for (int reciprocal_axis = 0; reciprocal_axis < 3; ++reciprocal_axis) {
          for (int cartesian_axis = 0; cartesian_axis < 3; ++cartesian_axis) {
            cartesian[cartesian_axis] +=
                reduced[reciprocal_axis]
                * reciprocal_cell[reciprocal_axis][cartesian_axis];
          }
        }
        const double norm2 =
            cartesian[0] * cartesian[0]
            + cartesian[1] * cartesian[1]
            + cartesian[2] * cartesian[2];
        if (norm2 < best_norm2 - 1.0e-14) {
          best_norm2 = norm2;
          best = {reduced, shift};
        }
      }
    }
  }
  return best;
}

std::vector<double> expand_k_blocks(
    const double* source,
    std::size_t unique_count,
    std::size_t block_size,
    const std::vector<int>& unique_index) {
  if (source == nullptr && unique_count * block_size != 0) {
    throw std::invalid_argument("k-fold source is null");
  }
  std::vector<double> expanded(unique_index.size() * block_size);
  for (std::size_t requested = 0; requested < unique_index.size(); ++requested) {
    const int source_index = unique_index[requested];
    if (source_index < 0 || static_cast<std::size_t>(source_index) >= unique_count) {
      throw std::out_of_range("k-fold unique index is outside the source array");
    }
    const double* source_block =
        source + static_cast<std::size_t>(source_index) * block_size;
    double* destination_block = expanded.data() + requested * block_size;
    std::copy_n(source_block, block_size, destination_block);
  }
  return expanded;
}

double* expand_k_blocks_owned(
    double* source,
    std::size_t unique_count,
    std::size_t block_size,
    const std::vector<int>& unique_index) {
  std::vector<double> expanded =
      expand_k_blocks(source, unique_count, block_size, unique_index);
  delete[] source;
  double* result = new double[expanded.size()];
  std::copy(expanded.begin(), expanded.end(), result);
  return result;
}

}  // namespace spadexp::kfold
