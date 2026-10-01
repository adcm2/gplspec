#ifndef GPLSPEC_SHARED_UTILITIES_H
#define GPLSPEC_SHARED_UTILITIES_H

#include <Eigen/Dense>
#include <Interpolation/All>
#include <cmath>
#include <complex>
#include <vector>

namespace GPLSpec::detail {

// Index a non-negative-order GSHTrans scalar coefficient in the triangular
// (l,m >= 0) layout used by the real-field conversion loops.
inline int NonNegativeHarmonicIndex(int l, int m) {
   if (m < 0) {
      return (l * (l + 1)) / 2 - m;
   } else {
      return (l * (l + 1)) / 2 + m;
   }
}

// Expand GSHTrans non-negative-m coefficients into the full real-field layout:
// for each degree, negative m precedes m >= 0, with the legacy Condon-Shortley
// conjugation factor retained exactly.
inline void ExpandRealScalarCoefficients(
    int lMax, const std::vector<std::complex<double>> &nonnegative,
    std::vector<std::complex<double>> &full) {
   int idx = 0;
   for (int idxl = 0; idxl < lMax + 1; ++idxl) {
      for (int idxm = -idxl; idxm < 0; ++idxm) {
         int idxlm = NonNegativeHarmonicIndex(idxl, idxm);
         full[idx] = std::pow(-1.0, idxm) * std::conj(nonnegative[idxlm]);
         ++idx;
      }
      for (int idxm = 0; idxm < idxl + 1; ++idxm) {
         int idxlm = NonNegativeHarmonicIndex(idxl, idxm);
         full[idx] = nonnegative[idxlm];
         ++idx;
      }
   }
}

// Scalar factor used by the existing rotation-output routines.
inline double RotationHarmonicNormalization(int l) {
   return std::sqrt((4.0 * 3.1415926535) / (2 * l + 1));
}

// Map a reference coordinate x in [-1,1] onto a physical radial interval.
template <class X, class X1, class X2>
auto StandardIntervalMap(const X &x, const X1 &x1, const X2 &x2)
    -> decltype(((x2 - x1) * x + (x1 + x2)) * 0.5) {
   return ((x2 - x1) * x + (x1 + x2)) * 0.5;
}

// Map x to the physical radius divided by the outer radius, using the ratio
// inner_radius / outer_radius used by the existing gravity kernels.
inline double ScaledRadialNodeMap(double x, double inner_over_outer) {
   return ((1.0 - inner_over_outer) * x + (1.0 + inner_over_outer)) * 0.5;
}

// Construct the Gauss-point derivative matrix on the supplied quadrature.
// polynomial_order is passed explicitly to preserve existing model loop bounds.
template <class Quadrature>
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
GaussDerivativeMatrix(const Quadrature &q, int polynomial_order) {
   const int matrix_order = polynomial_order + 1;
   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> derivative;
   derivative.resize(matrix_order, matrix_order);
   auto pleg = Interpolation::LagrangePolynomial(q.Points().begin(),
                                                  q.Points().end());
   for (int idxi = 0; idxi < polynomial_order + 1; ++idxi) {
      for (int idxj = 0; idxj < polynomial_order + 1; ++idxj) {
         derivative(idxi, idxj) = pleg.Derivative(idxi, q.X(idxj));
      }
   }
   return derivative;
}

}   // namespace GPLSpec::detail

#endif
