#include "../gplspec/All"

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
void require(bool condition, const std::string &message) {
   if (!condition) {
      throw std::runtime_error(message);
   }
}

GPLSpec::detail::RotationOutputAngles original_output_angles(
    const std::vector<double> &vec_p1, const std::vector<double> &vec_p2) {
   double theta1 = vec_p1[0];
   double phi1 = vec_p1[1];
   double theta2 = vec_p2[0];
   double phi2 = vec_p2[1];
   double cosphi = std::cos(theta2) * std::cos(theta1) +
                   std::sin(theta2) * std::sin(theta1) * std::cos(phi2 - phi1);
   double sinphi = std::sqrt(1 - cosphi * cosphi);
   auto tmp1 = std::sin(theta2) * std::cos(theta1) * std::cos(phi2) -
               std::cos(theta2) * std::sin(theta1) * std::cos(phi1);
   auto tmp2 = std::cos(theta2) * std::sin(theta1) * std::sin(phi1) -
               std::sin(theta2) * std::cos(theta1) * std::sin(phi2);
   auto tmp3 = std::sin(theta2) * std::sin(theta1) * std::sin(phi2 - phi1);
   auto tmp4 = -std::cos(theta1) * cosphi + std::cos(theta2);
   auto tmp5 = std::cos(theta1) * sinphi;
   if (std::abs(tmp1) < std::numeric_limits<double>::epsilon()) tmp1 = 0.0;
   if (std::abs(tmp2) < std::numeric_limits<double>::epsilon()) tmp2 = 0.0;
   if (std::abs(tmp4) < std::numeric_limits<double>::epsilon()) tmp4 = 0.0;
   if (std::abs(tmp5) < std::numeric_limits<double>::epsilon()) tmp5 = 0.0;
   return {std::atan2(tmp1, tmp2), std::acos(tmp3 / sinphi),
           std::atan2(tmp4, tmp5)};
}

GPLSpec::detail::RotationOutputAngles original_slice_angles(
    const std::vector<double> &vec_p1, const std::vector<double> &vec_p2) {
   double theta1 = vec_p1[0], phi1 = vec_p1[1];
   double theta2 = vec_p2[0], phi2 = vec_p2[1];
   double cosphi = std::cos(theta2) * std::cos(theta1) +
                   std::sin(theta2) * std::sin(theta1) * std::cos(phi2 - phi1);
   double sinphi = std::sqrt(1 - cosphi * cosphi);
   auto tmp1 = std::sin(theta2) * std::cos(theta1) * std::cos(phi2) -
               std::cos(theta2) * std::sin(theta1) * std::cos(phi1);
   auto tmp2 = std::cos(theta2) * std::sin(theta1) * std::sin(phi1) -
               std::sin(theta2) * std::cos(theta1) * std::sin(phi2);
   auto tmp3 = std::sin(theta2) * std::sin(theta1) * std::sin(phi2 - phi1);
   auto tmp4 = std::cos(theta1) * cosphi - std::cos(theta2);
   auto tmp5 = std::cos(theta1) * sinphi;
   if (tmp1 < std::numeric_limits<double>::epsilon()) tmp1 = 0.0;
   if (tmp2 < std::numeric_limits<double>::epsilon()) tmp2 = 0.0;
   if (tmp4 < std::numeric_limits<double>::epsilon()) tmp4 = 0.0;
   if (tmp5 < std::numeric_limits<double>::epsilon()) tmp5 = 0.0;
   return {std::atan2(tmp1, tmp2), std::acos(tmp3 / sinphi),
           std::atan2(tmp4, tmp5)};
}

std::vector<Eigen::MatrixXcd> original_matrices(int lMax, double alpha,
                                                 double beta, double gamma) {
   auto vec_wig = std::vector<Eigen::MatrixXcd>(lMax + 1);
   auto wigtemp =
       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                        GSHTrans::Single, GSHTrans::ColumnMajor>(
           lMax, lMax, lMax, beta);
   for (int l = 0; l < lMax + 1; ++l) {
      Eigen::MatrixXcd mat_tmp = Eigen::MatrixXcd::Zero(2 * l + 1, 2 * l + 1);
      auto multval = GPLSpec::detail::RotationHarmonicNormalization(l);
      for (int m = -l; m < l + 1; ++m) {
         auto dl = wigtemp[m];
         for (int mp = -l; mp < l + 1; ++mp) {
            std::complex<double> i1(0.0, 1.0);
            auto tmpmult =
                multval * exp(i1 * (static_cast<double>(m) * gamma +
                                    static_cast<double>(mp) * alpha));
            mat_tmp(m + l, mp + l) = dl[l, mp] * tmpmult;
         }
      }
      vec_wig[l] = mat_tmp;
   }
   return vec_wig;
}

struct OrderedGrid {
   void ForwardTransformation(int, int, const std::vector<double> &spatial,
                              std::vector<std::complex<double>> &coefficients) const {
      calls.push_back(spatial.front() < 10.0 ? 1 : 2);
      for (std::size_t i = 0; i < coefficients.size(); ++i) {
         coefficients[i] = {spatial.front() + static_cast<double>(i), 0.0};
      }
   }
   mutable std::vector<int> calls;
};
}  // namespace

int main() {
   const std::vector<std::pair<std::vector<double>, std::vector<double>>> cases{
       {{0.7, 0.2}, {1.1, 2.0}},
       {{1.3, 2.8}, {0.4, -0.7}},
       {{0.0, 0.0}, {1.2, 0.0}},
       {{1.5707963267948966, 0.0}, {1.5707963267948966, 1.5707963267948966}}};
   for (const auto &[p1, p2] : cases) {
      const auto expected = original_output_angles(p1, p2);
      const auto actual = GPLSpec::detail::RotationOutputAnglesFor(p1, p2);
      require(actual.alpha == expected.alpha && actual.beta == expected.beta &&
                  actual.gamma == expected.gamma,
              "output angle arithmetic differs from the original");
      const auto slice = original_slice_angles(p1, p2);
      if (p1[0] == 0.7) {
         // This fixture tests whether the original slice and output conventions are distinguishable.
         // Production RotateSliceToEquator output is covered by the separate baseline output test.
         require(slice.alpha != actual.alpha || slice.gamma != actual.gamma,
                 "angle fixture does not distinguish original slice and output conventions");
      }
      for (int lMax : {0, 1, 3, 5}) {
         const auto expected_matrices =
             original_matrices(lMax, expected.alpha, expected.beta,
                               expected.gamma);
         const auto actual_matrices = GPLSpec::detail::RotationOutputMatrices(
             lMax, actual.alpha, actual.beta, actual.gamma);
         require(actual_matrices.size() == expected_matrices.size(),
                 "rotation matrix degree count changed");
         for (std::size_t l = 0; l < actual_matrices.size(); ++l) {
            require((actual_matrices[l].array() ==
                     expected_matrices[l].array()).all(),
                    "rotation matrix differs from the original calculation");
         }
      }
   }

   OrderedGrid grid;
   std::vector<double> mapping{2.0, 3.0}, density{20.0, 21.0};
   std::vector<std::complex<double>> map_nonnegative(6), density_nonnegative(6);
   std::vector<std::complex<double>> map_full(9), density_full(9);
   GPLSpec::detail::ForwardAndExpandRealScalarPair(
       grid, 2, mapping, density, map_nonnegative, density_nonnegative,
       map_full, density_full);
   require(grid.calls == std::vector<int>({1, 2}),
           "paired transformations changed mapping/density call order");
   std::vector<std::complex<double>> expected_map(9), expected_density(9);
   GPLSpec::detail::ExpandRealScalarCoefficients(2, map_nonnegative,
                                                 expected_map);
   GPLSpec::detail::ExpandRealScalarCoefficients(2, density_nonnegative,
                                                 expected_density);
   require(map_full == expected_map && density_full == expected_density,
           "paired transform expansion differs from the original operation");
}
