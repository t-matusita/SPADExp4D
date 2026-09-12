#include "kfold.hpp"

#include <H5Cpp.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct CoefficientDataset {
  std::string atom;
  std::string orbital;
  std::vector<hsize_t> dimensions;
  std::vector<double> direct;
  std::vector<double> folded;
};

std::vector<hsize_t> dimensions(const H5::DataSet& dataset) {
  H5::DataSpace space = dataset.getSpace();
  std::vector<hsize_t> result(space.getSimpleExtentNdims());
  space.getSimpleExtentDims(result.data());
  return result;
}

std::size_t element_count(const std::vector<hsize_t>& dims,
                          std::size_t start = 0) {
  return std::accumulate(dims.begin() + static_cast<std::ptrdiff_t>(start),
                         dims.end(), std::size_t{1},
                         [](std::size_t a, hsize_t b) {
                           return a * static_cast<std::size_t>(b);
                         });
}

std::vector<double> read_double_dataset(const H5::DataSet& dataset) {
  std::vector<double> result(element_count(dimensions(dataset)));
  dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE);
  return result;
}

spadexp::kfold::Vector3 cross(const spadexp::kfold::Vector3& a,
                             const spadexp::kfold::Vector3& b) {
  return {a[1] * b[2] - a[2] * b[1],
          a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

double dot(const spadexp::kfold::Vector3& a,
           const spadexp::kfold::Vector3& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

spadexp::kfold::ReciprocalCell load_reciprocal_cell(H5::H5File& input) {
  H5::DataSet dataset = input.openDataSet("/Input/UnitCells/Bands");
  const std::vector<double> values = read_double_dataset(dataset);
  if (values.size() != 9) {
    throw std::runtime_error("band unit cell must be 3x3");
  }
  spadexp::kfold::ReciprocalCell real{};
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t column = 0; column < 3; ++column) {
      real[row][column] = values[3 * row + column];
    }
  }
  constexpr double two_pi = 6.283185307179586476925286766559;
  const double determinant = dot(real[0], cross(real[1], real[2]));
  return {{
      [&] { auto v = cross(real[1], real[2]); for (double& x : v) x *= two_pi / determinant; return v; }(),
      [&] { auto v = cross(real[2], real[0]); for (double& x : v) x *= two_pi / determinant; return v; }(),
      [&] { auto v = cross(real[0], real[1]); for (double& x : v) x *= two_pi / determinant; return v; }(),
  }};
}

std::vector<std::string> child_names(const H5::Group& group,
                                     H5G_obj_t requested_type) {
  std::vector<std::string> result;
  for (hsize_t i = 0; i < group.getNumObjs(); ++i) {
    if (group.getObjTypeByIdx(i) == requested_type) {
      result.push_back(group.getObjnameByIdx(i));
    }
  }
  std::sort(result.begin(), result.end());
  return result;
}

std::vector<CoefficientDataset> load_coefficients(
    H5::H5File& input, std::size_t point_count,
    const std::vector<int>& unique_index) {
  H5::Group lcao = input.openGroup("/Output/LCAO");
  std::vector<CoefficientDataset> result;
  for (const std::string& atom : child_names(lcao, H5G_GROUP)) {
    H5::Group atom_group = lcao.openGroup(atom);
    for (const std::string& orbital : child_names(atom_group, H5G_DATASET)) {
      H5::DataSet dataset = atom_group.openDataSet(orbital);
      CoefficientDataset item;
      item.atom = atom;
      item.orbital = orbital;
      item.dimensions = dimensions(dataset);
      if (item.dimensions.size() != 4 || item.dimensions[0] != point_count) {
        throw std::runtime_error("unexpected LCAO dimensions at " + atom + "/" + orbital);
      }
      item.direct = read_double_dataset(dataset);
      const std::size_t block_size = element_count(item.dimensions, 1);
      item.folded = spadexp::kfold::expand_k_blocks(
          item.direct.data(), 1, block_size, unique_index);
      result.push_back(std::move(item));
    }
  }
  return result;
}

std::size_t offset(const CoefficientDataset& item, std::size_t point,
                   std::size_t band, std::size_t magnetic,
                   std::size_t digit) {
  const std::size_t band_count = item.dimensions[1];
  const std::size_t magnetic_count = item.dimensions[2];
  const std::size_t digit_count = item.dimensions[3];
  return (((point * band_count + band) * magnetic_count + magnetic)
           * digit_count + digit);
}

std::complex<double> coefficient(const std::vector<double>& values,
                                 const CoefficientDataset& item,
                                 std::size_t point, std::size_t band,
                                 std::size_t magnetic, std::size_t pair) {
  const std::size_t first_digit = 2 * pair;
  return {values[offset(item, point, band, magnetic, first_digit)],
          values[offset(item, point, band, magnetic, first_digit + 1)]};
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
              << " direct_openmx.hdf5 coefficient_report.txt\n";
    return 2;
  }

  try {
    H5::Exception::dontPrint();
    H5::H5File input(argv[1], H5F_ACC_RDONLY);
    H5::DataSet band_dataset = input.openDataSet("/Output/Band");
    const std::vector<hsize_t> band_dims = dimensions(band_dataset);
    if (band_dims.size() != 2 || band_dims[0] != 4) {
      throw std::runtime_error("the direct calculation must contain four k points");
    }
    const std::size_t point_count = 4;
    const std::size_t band_count = static_cast<std::size_t>(band_dims[1]);
    const std::vector<double> direct_bands = read_double_dataset(band_dataset);
    const spadexp::kfold::ReciprocalCell reciprocal_cell =
        load_reciprocal_cell(input);
    const std::vector<spadexp::kfold::Vector3> requested = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},
    };
    std::vector<spadexp::kfold::FoldResult> fold_results;
    std::vector<int> gamma_for_every_point;
    for (const auto& point : requested) {
      const auto result = spadexp::kfold::fold_to_first_bz(point, reciprocal_cell);
      if (std::sqrt(dot(result.reduced, result.reduced)) > 1.0e-12) {
        throw std::runtime_error("G vector did not fold to Gamma");
      }
      fold_results.push_back(result);
      gamma_for_every_point.push_back(0);
    }
    const std::vector<double> folded_bands = spadexp::kfold::expand_k_blocks(
        direct_bands.data(), 1, band_count, gamma_for_every_point);
    const std::vector<CoefficientDataset> coefficients =
        load_coefficients(input, point_count, gamma_for_every_point);

    const char* labels[] = {
        "k=(0,0,0)",
        "k=(G1,0,0)  [reciprocal coordinates (1,0,0)]",
        "k=(0,G2,0)  [reciprocal coordinates (0,1,0)]",
        "k=(0,0,G3)  [reciprocal coordinates (0,0,1)]",
    };

    std::ofstream output(argv[2]);
    if (!output) {
      throw std::runtime_error("cannot open output report");
    }
    output << std::scientific << std::setprecision(10);
    output << "# OpenMX direct coefficients vs coefficients expanded from first-BZ Gamma\n";
    output << "# Fold map: requested [Gamma,G1,G2,G3] -> unique first-BZ index [0,0,0,0]\n";
    output << "# Folding coordinates and Folded columns are produced by kfold.cpp.\n";
    output << "# requested -> reduced ; reciprocal shift\n";
    for (std::size_t point = 0; point < point_count; ++point) {
      output << "# (" << requested[point][0] << ',' << requested[point][1] << ','
             << requested[point][2] << ") -> ("
             << fold_results[point].reduced[0] << ','
             << fold_results[point].reduced[1] << ','
             << fold_results[point].reduced[2] << ") ; G=("
             << fold_results[point].shift[0] << ','
             << fold_results[point].shift[1] << ','
             << fold_results[point].shift[2] << ")\n";
    }
    output << '\n';

    output << "SUMMARY\n";
    for (std::size_t point = 0; point < point_count; ++point) {
      double band_max = 0.0;
      double coefficient_max = 0.0;
      double overlap_min = 1.0;
      for (std::size_t band = 0; band < band_count; ++band) {
        std::complex<double> inner{0.0, 0.0};
        double direct_norm = 0.0;
        double folded_norm = 0.0;
        for (const CoefficientDataset& item : coefficients) {
          const std::size_t magnetic_count = item.dimensions[2];
          const std::size_t pair_count = item.dimensions[3] / 2;
          for (std::size_t magnetic = 0; magnetic < magnetic_count; ++magnetic) {
            for (std::size_t pair = 0; pair < pair_count; ++pair) {
              const auto direct = coefficient(item.direct, item, point, band, magnetic, pair);
              const auto folded = coefficient(item.folded, item, point, band, magnetic, pair);
              inner += std::conj(direct) * folded;
              direct_norm += std::norm(direct);
              folded_norm += std::norm(folded);
              coefficient_max = std::max(coefficient_max, std::abs(direct - folded));
            }
          }
        }
        if (direct_norm > 0.0 && folded_norm > 0.0) {
          overlap_min = std::min(
              overlap_min, std::abs(inner) / std::sqrt(direct_norm * folded_norm));
        }
        band_max = std::max(
            band_max,
            std::abs(direct_bands[point * band_count + band]
                     - folded_bands[point * band_count + band]));
      }
      output << labels[point] << "\n"
             << "  max_abs_band_difference = " << band_max << "\n"
             << "  max_abs_coefficient_difference = " << coefficient_max << "\n"
             << "  minimum_normalized_overlap = " << overlap_min << "\n";
    }

    for (std::size_t point = 0; point < point_count; ++point) {
      output << "\n============================================================\n";
      output << labels[point] << "\n";
      output << "folded_k=(" << fold_results[point].reduced[0] << ','
             << fold_results[point].reduced[1] << ','
             << fold_results[point].reduced[2] << ")  G_shift=("
             << fold_results[point].shift[0] << ','
             << fold_results[point].shift[1] << ','
             << fold_results[point].shift[2] << ")\n";
      output << "atom orbital band m pair  OpenMX_Re OpenMX_Im  "
                "Fold_Re Fold_Im  Diff_Abs\n";
      for (const CoefficientDataset& item : coefficients) {
        const std::size_t magnetic_count = item.dimensions[2];
        const std::size_t pair_count = item.dimensions[3] / 2;
        for (std::size_t band = 0; band < band_count; ++band) {
          for (std::size_t magnetic = 0; magnetic < magnetic_count; ++magnetic) {
            for (std::size_t pair = 0; pair < pair_count; ++pair) {
              const auto direct = coefficient(item.direct, item, point, band, magnetic, pair);
              const auto folded = coefficient(item.folded, item, point, band, magnetic, pair);
              output << item.atom << ' ' << item.orbital << ' '
                     << band << ' ' << magnetic << ' ' << pair << "  "
                     << direct.real() << ' ' << direct.imag() << "  "
                     << folded.real() << ' ' << folded.imag() << "  "
                     << std::abs(direct - folded) << '\n';
            }
          }
        }
      }
    }
    std::cout << "Wrote " << argv[2] << '\n';
  } catch (const H5::Exception& error) {
    std::cerr << "HDF5 error: " << error.getDetailMsg() << '\n';
    return 1;
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
