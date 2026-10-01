#include "../gplspec/src/Shared_Utilities.h"
#include "../gplspec/src/Spectral_Element_Tools.h"
#include "../gplspec/src/Spherical_Integrator.h"
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <complex>
#include <vector>

namespace {

void Require(bool condition) {
   if (!condition) {
      throw std::runtime_error("stage-03 utility comparison failed");
   }
}

int OriginalNonNegativeHarmonicIndex(int l, int m) {
   if (m < 0) {
      return (l * (l + 1)) / 2 - m;
   } else {
      return (l * (l + 1)) / 2 + m;
   }
}

void OriginalExpandRealScalarCoefficients(
    int lMax, const std::vector<std::complex<double>> &nonnegative,
    std::vector<std::complex<double>> &full) {
   int idx = 0;
   for (int idxl = 0; idxl < lMax + 1; ++idxl) {
      for (int idxm = -idxl; idxm < 0; ++idxm) {
         int idxlm = OriginalNonNegativeHarmonicIndex(idxl, idxm);
         full[idx] = std::pow(-1.0, idxm) * std::conj(nonnegative[idxlm]);
         ++idx;
      }
      for (int idxm = 0; idxm < idxl + 1; ++idxm) {
         int idxlm = OriginalNonNegativeHarmonicIndex(idxl, idxm);
         full[idx] = nonnegative[idxlm];
         ++idx;
      }
   }
}

template <class Quadrature>
Eigen::MatrixXd OriginalGaussDerivativeMatrix(const Quadrature &q,
                                               int polynomial_order) {
   Eigen::MatrixXd derivative;
   derivative.resize(polynomial_order + 1, polynomial_order + 1);
   auto pleg = Interpolation::LagrangePolynomial(q.Points().begin(),
                                                  q.Points().end());
   for (int idxi = 0; idxi < polynomial_order + 1; ++idxi) {
      for (int idxj = 0; idxj < polynomial_order + 1; ++idxj) {
         derivative(idxi, idxj) = pleg.Derivative(idxi, q.X(idxj));
      }
   }
   return derivative;
}

template <class Float>
Float OriginalIntervalMap(const Float &x, const Float &x1, const Float &x2) {
   return ((x2 - x1) * x + (x1 + x2)) * 0.5;
}

double OriginalScaledRadialNodeMap(double x, double inner_over_outer) {
   return ((1.0 - inner_over_outer) * x + (1.0 + inner_over_outer)) * 0.5;
}

}   // namespace

int main() {
   for (int lMax = 0; lMax <= 12; ++lMax) {
      const int nonnegative_size = (lMax + 1) * (lMax + 2) / 2;
      std::vector<std::complex<double>> nonnegative(nonnegative_size);
      for (int idx = 0; idx < nonnegative_size; ++idx) {
         nonnegative[idx] = {0.25 * (idx + 1), -0.125 * (idx + 2)};
      }
      std::vector<std::complex<double>> expected((lMax + 1) * (lMax + 1));
      std::vector<std::complex<double>> actual(expected.size());
      OriginalExpandRealScalarCoefficients(lMax, nonnegative, expected);
      GPLSpec::detail::ExpandRealScalarCoefficients(lMax, nonnegative, actual);
      Require(actual == expected);
      for (int l = 0; l <= lMax; ++l) {
         for (int m = -l; m <= l; ++m) {
            Require(GPLSpec::detail::NonNegativeHarmonicIndex(l, m) ==
                   OriginalNonNegativeHarmonicIndex(l, m));
         }
         Require(GPLSpec::detail::RotationHarmonicNormalization(l) ==
                std::sqrt((4.0 * 3.1415926535) / (2 * l + 1)));
      }
   }

   for (int node_count = 3; node_count <= 9; ++node_count) {
      auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(node_count);
      const int polynomial_order = node_count - 1;
      const auto expected = OriginalGaussDerivativeMatrix(q, polynomial_order);
      const auto actual =
          GPLSpec::detail::GaussDerivativeMatrix(q, polynomial_order);
      Require(actual.rows() == expected.rows());
      Require(actual.cols() == expected.cols());
      Require((actual.array() == expected.array()).all());
   }

   const std::vector<double> interval_edges{0.0, 0.2, 0.73, 1.0, 1.2};
   const std::vector<double> reference_nodes{-1.0, -0.75, 0.0, 0.625, 1.0};
   for (std::size_t idx = 0; idx < interval_edges.size() - 1; ++idx) {
      for (double x : reference_nodes) {
         const auto expected = OriginalIntervalMap(
             x, interval_edges[idx], interval_edges[idx + 1]);
         Require(GPLSpec::detail::StandardIntervalMap(
                    x, interval_edges[idx], interval_edges[idx + 1]) ==
                expected);
         Require(Radial_Tools::StandardIntervalMap(
                    x, interval_edges[idx], interval_edges[idx + 1]) ==
                expected);
         Require(GravityFunctions::StandardIntervalMap(
                    x, interval_edges[idx], interval_edges[idx + 1]) ==
                expected);
      }
   }

   for (double ratio : {0.0, 0.125, 0.5, 0.875, 1.0}) {
      for (double x : reference_nodes) {
         Require(GPLSpec::detail::ScaledRadialNodeMap(x, ratio) ==
                OriginalScaledRadialNodeMap(x, ratio));
      }
   }

   const float lower_float = 0.125F;
   const float upper_float = 1.375F;
   const double reference_double = -0.375;
   const auto mixed_expected =
       ((upper_float - lower_float) * reference_double +
        (lower_float + upper_float)) *
       0.5;
   Require(GPLSpec::detail::StandardIntervalMap(
              reference_double, lower_float, upper_float) == mixed_expected);
}
