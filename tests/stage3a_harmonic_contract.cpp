#include <gplspec/All>
#include <GSHTrans/All>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Real = double;
using Complex = std::complex<Real>;
using Grid = GSHTrans::GaussLegendreGrid<Real, GSHTrans::All, GSHTrans::All>;

void require(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error(message);
}

void near(Complex actual, Complex expected, Real tolerance,
          const std::string &message) {
  if (!std::isfinite(actual.real()) || !std::isfinite(actual.imag()) ||
      !std::isfinite(expected.real()) || !std::isfinite(expected.imag()) ||
      std::abs(actual - expected) > tolerance * std::max(1.0, std::abs(expected))) {
    std::ostringstream detail;
    detail << std::setprecision(17) << message << ": actual=" << actual
           << " expected=" << expected << " difference="
           << std::abs(actual - expected) << " tolerance="
           << tolerance * std::max(1.0, std::abs(expected));
    throw std::runtime_error(detail.str());
  }
}

void check_indices() {
  constexpr int lmax = 4;
  const GSHTrans::GSHIndices<GSHTrans::All> full(lmax, lmax, 0);
  const GSHTrans::GSHIndices<GSHTrans::NonNegative> half(lmax, lmax, 0);
  for (int l = 0; l <= lmax; ++l) {
    for (int m = -l; m <= l; ++m)
      require(full.Index(l, m) == l * l + l + m,
              "full-order GSHTrans index differs from GPLSpec layout");
    for (int m = 0; m <= l; ++m)
      require(half.Index(l, m) == l * (l + 1) / 2 + m,
              "nonnegative-order GSHTrans index differs from triangular layout");
  }
  for (int spin : {-1, 1}) {
    const GSHTrans::GSHIndices<GSHTrans::All> indices(lmax, lmax, spin);
    for (int l = 1; l <= lmax; ++l)
      for (int m = -l; m <= l; ++m)
        require(indices.Index(l, m) == l * l - 1 + l + m,
                "spin +/-1 GSHTrans index has an unexpected offset");
  }
}

void check_scalar_analytic_identities() {
  constexpr int lmax = 3;
  Grid grid(lmax, 1, FFTWpp::Estimate);
  const GSHTrans::GSHIndices<GSHTrans::NonNegative> indices(lmax, lmax, 0);
  std::vector<Real> samples(grid.FieldSize());
  std::vector<Complex> coefficients(indices.Size());
  constexpr Real tolerance = 3e-13;
  const Real pi = std::numbers::pi_v<Real>;

  for (std::ptrdiff_t itheta = 0; itheta < grid.NumberOfCoLatitudes(); ++itheta) {
    const Real theta = *(grid.CoLatitudes().begin() + itheta);
    for (std::ptrdiff_t iphi = 0; iphi < grid.NumberOfLongitudes(); ++iphi) {
      const auto offset = itheta * grid.NumberOfLongitudes() + iphi;
      const Real phi = *(grid.Longitudes().begin() + iphi);
      samples[offset] = 1.0 / std::sqrt(4.0 * pi);
    }
  }
  grid.ForwardTransformation(lmax, 0, samples, coefficients);
  near(coefficients[indices.Index(0, 0)], Complex{1.0, 0.0}, tolerance,
       "analytic Y00 did not transform to unit coefficient");
  for (std::size_t i = 1; i < coefficients.size(); ++i)
    near(coefficients[i], Complex{}, tolerance,
         "constant field leaked outside the Y00 coefficient");

  std::fill(coefficients.begin(), coefficients.end(), Complex{});
  for (std::ptrdiff_t itheta = 0; itheta < grid.NumberOfCoLatitudes(); ++itheta) {
    const Real theta = *(grid.CoLatitudes().begin() + itheta);
    for (std::ptrdiff_t iphi = 0; iphi < grid.NumberOfLongitudes(); ++iphi) {
      const auto offset = itheta * grid.NumberOfLongitudes() + iphi;
      samples[offset] = std::sqrt(3.0 / (4.0 * pi)) * std::cos(theta);
    }
  }
  grid.ForwardTransformation(lmax, 0, samples, coefficients);
  near(coefficients[indices.Index(1, 0)], Complex{1.0, 0.0}, tolerance,
       "analytic Y10 did not transform to unit coefficient");
  for (std::size_t i = 0; i < coefficients.size(); ++i)
    if (i != indices.Index(1, 0))
      near(coefficients[i], Complex{}, tolerance,
           "analytic Y10 leaked outside its coefficient");

  std::fill(coefficients.begin(), coefficients.end(), Complex{});
  for (std::ptrdiff_t itheta = 0; itheta < grid.NumberOfCoLatitudes(); ++itheta) {
    const Real theta = *(grid.CoLatitudes().begin() + itheta);
    for (std::ptrdiff_t iphi = 0; iphi < grid.NumberOfLongitudes(); ++iphi) {
      const auto offset = itheta * grid.NumberOfLongitudes() + iphi;
      const Real phi = *(grid.Longitudes().begin() + iphi);
      samples[offset] = std::sqrt(3.0 / (4.0 * pi)) *
                        std::sin(theta) * std::cos(phi);
    }
  }
  grid.ForwardTransformation(lmax, 0, samples, coefficients);
  near(coefficients[indices.Index(1, 1)], Complex{-1.0 / std::sqrt(2.0), 0.0},
       tolerance, "analytic real Y11 field has unexpected Condon-Shortley sign");
  for (std::size_t i = 0; i < coefficients.size(); ++i)
    if (i != indices.Index(1, 1))
      near(coefficients[i], Complex{}, tolerance,
           "analytic real Y11 field leaked outside its positive-m coefficient");
  std::vector<Complex> full((lmax + 1) * (lmax + 1));
  GPLSpec::detail::ExpandRealScalarCoefficients(lmax, coefficients, full);
  near(full[1 * 1 + 1 - 1], Complex{1.0 / std::sqrt(2.0), 0.0}, tolerance,
       "real-field expansion has an unexpected negative-m conjugacy sign");
  near(full[1 * 1 + 1 + 1], Complex{-1.0 / std::sqrt(2.0), 0.0}, tolerance,
       "real-field expansion changed the positive-m coefficient");
}

void check_spin_one_analytic_identities() {
  constexpr int lmax = 3;
  Grid grid(lmax, 1, FFTWpp::Estimate);
  const Real pi = std::numbers::pi_v<Real>;
  std::vector<Complex> samples(grid.FieldSize());
  for (std::ptrdiff_t itheta = 0; itheta < grid.NumberOfCoLatitudes(); ++itheta) {
    const Real theta = *(grid.CoLatitudes().begin() + itheta);
    for (std::ptrdiff_t iphi = 0; iphi < grid.NumberOfLongitudes(); ++iphi)
      samples[itheta * grid.NumberOfLongitudes() + iphi] =
          std::sqrt(3.0 / (8.0 * pi)) * std::sin(theta);
  }

  for (int spin : {-1, 1}) {
    const GSHTrans::GSHIndices<GSHTrans::All> indices(lmax, lmax, spin);
    std::vector<Complex> coefficients(indices.Size());
    grid.ForwardTransformation(lmax, spin, samples, coefficients);
    const Complex expected = spin == 1 ? Complex{1.0, 0.0}
                                       : Complex{-1.0, 0.0};
    near(coefficients[indices.Index(1, 0)], expected, 5e-13,
         "analytic spin-one l=1,m=0 coefficient has unexpected sign or scale");
    for (int l = 1; l <= lmax; ++l)
      for (int m = -l; m <= l; ++m)
        if (l != 1 || m != 0)
          near(coefficients[indices.Index(l, m)], Complex{}, 5e-13,
               "analytic spin-one field leaked outside l=1,m=0");

    std::vector<Complex> one_hot(indices.Size());
    one_hot[indices.Index(1, 0)] = 1.0;
    std::vector<Complex> reconstructed(grid.FieldSize());
    std::vector<Complex> roundtrip(indices.Size());
    grid.InverseTransformation(lmax, spin, one_hot, reconstructed);
    const Real spin_sign = spin == 1 ? 1.0 : -1.0;
    for (std::ptrdiff_t itheta = 0; itheta < grid.NumberOfCoLatitudes(); ++itheta) {
      const Real theta = *(grid.CoLatitudes().begin() + itheta);
      const Complex expected_field{
          spin_sign * std::sqrt(3.0 / (8.0 * pi)) * std::sin(theta), 0.0};
      for (std::ptrdiff_t iphi = 0; iphi < grid.NumberOfLongitudes(); ++iphi)
        near(reconstructed[itheta * grid.NumberOfLongitudes() + iphi],
             expected_field, 5e-13,
             "analytic spin-one inverse field differs at grid point");
    }
    grid.ForwardTransformation(lmax, spin, reconstructed, roundtrip);
    near(roundtrip[indices.Index(1, 0)], Complex{1.0, 0.0}, 5e-13,
         "spin +/-1 one-hot transform roundtrip changed normalization");
  }
}

void check_nyquist_and_rotation_contract() {
  constexpr int lmax = 3;
  Grid grid(lmax, lmax, FFTWpp::Estimate);
  const auto all = GSHTrans::GSHIndices<GSHTrans::All>(lmax, lmax, 0);
  std::vector<Complex> samples(grid.FieldSize());
  for (std::size_t i = 0; i < samples.size(); ++i)
    samples[i] = Complex{static_cast<Real>(i + 1), -0.25 * static_cast<Real>(i)};
  std::vector<Complex> coefficients(all.Size());
  grid.ForwardTransformation(lmax, 0, samples, coefficients);
  require(coefficients[all.Index(lmax, lmax)] == Complex{},
          "pinned complex transform no longer zeros the +lmax Nyquist slot");
  require(grid.NumberOfLongitudes() == 2 * lmax,
          "pinned longitude count no longer equals 2*lmax");

  const Real scale =
      std::sqrt(4.0 * 3.1415926535 / 3.0) *
      (std::numbers::inv_sqrtpi_v<Real> / 2.0) * std::sqrt(3.0);
  const Real half_turn = std::numbers::pi_v<Real> / 2.0;
  const auto identity = GPLSpec::detail::RotationOutputMatrices(1, 0.0, 0.0, 0.0);
  const auto z_rotation = GPLSpec::detail::RotationOutputMatrices(1, 0.0, 0.0, half_turn);
  const auto &id = identity[1];
  const auto &rot = z_rotation[1];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      near(id(i, j), i == j ? Complex{scale, 0.0} : Complex{}, 5e-13,
           "degree-one identity rotation matrix entry " + std::to_string(i) +
               "," + std::to_string(j));
      const Complex expected = i == j
          ? scale * std::polar(1.0, static_cast<Real>(i - 1) * half_turn)
          : Complex{};
      near(rot(i, j), expected, 5e-13,
           "degree-one z-axis rotation differs from analytic phase matrix");
    }
  const auto y_rotation = GPLSpec::detail::RotationOutputMatrices(
      1, 0.0, half_turn, 0.0);
  const Real alpha = 0.31;
  const Real beta = 0.71;
  const Real gamma = -0.23;
  const auto wrapped = GPLSpec::detail::RotationOutputMatrices(
      1, alpha, beta, gamma);
  const GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All,
                         GSHTrans::All, GSHTrans::Single,
                         GSHTrans::ColumnMajor>
      pinned(1, 1, 1, beta);
  for (int l = 0; l <= 1; ++l) {
    const auto normalization =
        GPLSpec::detail::RotationHarmonicNormalization(l);
    int row = 0;
    for (int m = -l; m <= l; ++m, ++row) {
      const auto d = pinned[m];
      int column = 0;
      for (int mp = -l; mp <= l; ++mp, ++column) {
        const Complex phase =
            std::exp(Complex{0.0, static_cast<Real>(m) * gamma +
                                      static_cast<Real>(mp) * alpha});
        near(wrapped[l](row, column),
             normalization * d[l][mp] * phase, 5e-13,
             "GPLSpec rotation wrapper differs from pinned GSHTrans Wigner");
      }
    }
  }
  const auto gram = y_rotation[1] * y_rotation[1].adjoint();
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      near(gram(i, j), i == j ? Complex{scale * scale, 0.0} : Complex{},
           5e-13, "degree-one y rotation violates analytic unitarity");
}
}  // namespace

int main() {
  check_indices();
  check_scalar_analytic_identities();
  check_spin_one_analytic_identities();
  check_nyquist_and_rotation_contract();
}
