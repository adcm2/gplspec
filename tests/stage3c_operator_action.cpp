#include <gplspec/All>
#include <gplspec/src/Pseudospectral_Matrix_Wrapper3D.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
using Complex = std::complex<double>;
using Grid = MatrixReplacement<Complex>::Grid;
using TensorField = GPLSpec::OperatorActions::TensorField3D;

std::size_t scalar_index(int degree, int order) {
  return static_cast<std::size_t>(degree * degree + degree + order);
}
std::size_t spin_index(int degree, int order) {
  return static_cast<std::size_t>(degree * degree + degree + order - 1);
}
double omega(int degree) {
  return std::sqrt(static_cast<double>(degree) * (degree + 1) / 2.0);
}

template <class G, class Spatial>
Complex direct_transform(const G &grid, int lmax, int spin,
                         const Spatial &spatial, int degree,
                         int order) {
  // GSHTrans reserves the last (Nyquist) coefficient as an exact zero.
  if (degree == lmax && order == lmax) return Complex{};
  const std::size_t harmonics = static_cast<std::size_t>(lmax + 1) * (lmax + 1);
  const std::size_t spin_count =
      GSHTrans::GSHIndices<GSHTrans::All>(lmax, lmax, 1).Size();
  std::vector<Complex> basis0(harmonics, Complex{}), basis_spin(spin_count, Complex{});
  FFTWpp::vector<Complex> basis_spatial(
      grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes(), Complex{});
  if (spin == 0)
    basis0[scalar_index(degree, order)] = 1.0;
  else
    basis_spin[spin_index(degree, order)] = 1.0;
  const auto &basis = spin == 0 ? basis0 : basis_spin;
  grid.InverseTransformation(lmax, spin, basis, basis_spatial);

  Complex coefficient{};
  std::size_t sample = 0;
  for (std::size_t theta = 0; theta < grid.NumberOfCoLatitudes(); ++theta) {
    for (std::size_t phi = 0; phi < grid.NumberOfLongitudes(); ++phi) {
      const double weight = grid.CoLatitudeWeights()[theta] *
                            grid.LongitudeWeights()[phi];
      coefficient += std::conj(basis_spatial[sample]) * spatial[sample] * weight;
      ++sample;
    }
  }
  return coefficient;
}

template <class Model, class TensorAt>
Eigen::VectorXcd direct_weak_action(const Model &model, const Grid &grid,
                                    const TensorAt &tensor_at,
                                    const Eigen::VectorXcd &input,
                                    bool mapped_geometry, bool boundary) {
  const auto &mesh = model.Node_Information();
  const auto quadrature = model.q();
  const auto derivative = model.GaussDerivative();
  const int lmax = grid.MaxDegree();
  const std::size_t nh = static_cast<std::size_t>(lmax + 1) * (lmax + 1);
  const std::size_t ns = GSHTrans::GSHIndices<GSHTrans::All>(lmax, lmax, 1).Size();
  const std::size_t radial_count = mesh.NumberOfGlobalRadialNodes();
  if (input.size() != static_cast<Eigen::Index>(radial_count * nh))
    throw std::runtime_error("direct reference input size mismatch");
  Eigen::VectorXcd expected = Eigen::VectorXcd::Zero(input.size());
  const std::size_t spatial_count =
      grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();

  for (std::size_t element = 0; element < mesh.NumberOfElements(); ++element) {
    const double width = mesh.ElementWidth(element);
    for (std::size_t node = 0; node < mesh.NumNodesWithinElement(); ++node) {
      const double radius = 0.5 * ((mesh.ElementUpperRadius(element) -
                                      mesh.ElementLowerRadius(element)) *
                                     quadrature.X(node) +
                                     mesh.ElementUpperRadius(element) +
                                     mesh.ElementLowerRadius(element));
      std::vector<Complex> radial_gradient(nh, Complex{}),
          minus_gradient(ns, Complex{}), plus_gradient(ns, Complex{});
      for (std::size_t local = 0; local < mesh.NumNodesWithinElement(); ++local) {
        const auto radial = mesh.GlobalNodeIndex(element, local);
        const double scale = derivative(local, node) * 2.0 / width *
                             (mapped_geometry ? radius : 1.0);
        for (int degree = 0; degree <= lmax; ++degree)
          for (int order = -degree; order <= degree; ++order)
            radial_gradient[scalar_index(degree, order)] +=
                scale * input[radial * nh + scalar_index(degree, order)];
      }
      const auto radial = mesh.GlobalNodeIndex(element, node);
      for (int degree = 1; degree <= lmax; ++degree) {
        // In the 1D isotropic tensor, spin flux has no radial coupling; the
        // spin-gradient 1/r cancels the weak angular measure r exactly.
        const double scale = omega(degree);
        for (int order = -degree; order <= degree; ++order) {
          const auto spin = spin_index(degree, order);
          const Complex value = scale * input[radial * nh + scalar_index(degree, order)];
          minus_gradient[spin] = value;
          plus_gradient[spin] = value;
        }
      }

      FFTWpp::vector<Complex> spatial_minus(spatial_count, Complex{}),
          spatial_zero(spatial_count, Complex{}),
          spatial_plus(spatial_count, Complex{});
      grid.InverseTransformation(lmax, -1, minus_gradient, spatial_minus);
      grid.InverseTransformation(lmax, 0, radial_gradient, spatial_zero);
      grid.InverseTransformation(lmax, 1, plus_gradient, spatial_plus);
      FFTWpp::vector<Complex> flux_minus(spatial_count, Complex{}),
          flux_zero(spatial_count, Complex{}),
          flux_plus(spatial_count, Complex{});
      for (std::size_t point = 0; point < spatial_count; ++point) {
        const Eigen::Matrix3cd &a = tensor_at(element, node, point);
        // Independently state the canonical tensor signs rather than calling
        // the production contraction helper.
        flux_minus[point] = -a(0, 0) * spatial_plus[point] +
                            a(0, 1) * spatial_zero[point] -
                            a(0, 2) * spatial_minus[point];
        flux_zero[point] = -a(1, 0) * spatial_plus[point] +
                           a(1, 1) * spatial_zero[point] -
                           a(1, 2) * spatial_minus[point];
        flux_plus[point] = -a(2, 0) * spatial_plus[point] +
                           a(2, 1) * spatial_zero[point] -
                           a(2, 2) * spatial_minus[point];
      }

      std::vector<Complex> q_zero(nh), q_minus(ns), q_plus(ns);
      for (int degree = 0; degree <= lmax; ++degree) {
        for (int order = -degree; order <= degree; ++order) {
          const auto h = scalar_index(degree, order);
          q_zero[h] = direct_transform(grid, lmax, 0, flux_zero, degree, order);
          if (degree > 0) {
            const auto spin = spin_index(degree, order);
            q_minus[spin] = direct_transform(grid, lmax, -1, flux_minus,
                                             degree, order);
            q_plus[spin] = direct_transform(grid, lmax, 1, flux_plus,
                                            degree, order);
          }
        }
      }

      const double weight = quadrature.W(node);
      const double radial_measure = mapped_geometry ? radius : radius * radius;
      for (std::size_t local = 0; local < mesh.NumNodesWithinElement(); ++local) {
        const std::size_t output_radial = mesh.GlobalNodeIndex(element, local);
        const double scale = derivative(local, node) * weight * radial_measure;
        for (std::size_t h = 0; h < nh; ++h)
          expected[output_radial * nh + h] -= scale * q_zero[h];
      }
      for (int degree = 1; degree <= lmax; ++degree) {
        const double scale = width * 0.5 * weight *
                            omega(degree);
        for (int order = -degree; order <= degree; ++order) {
          const auto h = scalar_index(degree, order);
          const auto spin = spin_index(degree, order);
          expected[radial * nh + h] -=
              scale * (q_minus[spin] + q_plus[spin]);
        }
      }
    }
  }

  if (boundary) {
    const std::size_t radial = radial_count - 1;
    for (int degree = 0; degree <= lmax; ++degree)
      for (int order = -degree; order <= degree; ++order) {
        const auto index = radial * nh + scalar_index(degree, order);
        expected[index] -= (degree + 1.0) * mesh.OuterRadius() * input[index];
      }
  }
  return expected;
}

void require_close(const Eigen::VectorXcd &actual,
                   const Eigen::VectorXcd &expected,
                   const char *description) {
  const double error = (actual - expected).norm();
  const double scale = expected.norm();
  const double limit = 1e-11 * scale + 1e-13 * std::max(1.0, scale);
  if (!actual.allFinite() || !expected.allFinite() || error > limit) {
    Eigen::Index bad_actual = -1, bad_expected = -1;
    for (Eigen::Index i = 0; i < actual.size(); ++i) {
      if (!std::isfinite(actual[i].real()) || !std::isfinite(actual[i].imag())) { bad_actual = i; break; }
    }
    for (Eigen::Index i = 0; i < expected.size(); ++i) {
      if (!std::isfinite(expected[i].real()) || !std::isfinite(expected[i].imag())) { bad_expected = i; break; }
    }
    Eigen::Index max_index = 0;
    (actual - expected).cwiseAbs().maxCoeff(&max_index);
    std::cerr << description << ": error=" << error << " limit=" << limit
              << " scale=" << scale << " bad_actual=" << bad_actual
              << " bad_expected=" << bad_expected << " max_index=" << max_index
              << " actual=" << actual[max_index] << " expected=" << expected[max_index] << '\n';
    throw std::runtime_error(description);
  }
}

struct SmoothMapping {
  auto RadialMapping(int) const {
    return [](double r, double theta, double phi) {
      return 0.08 * r * (1.0 - r) * std::sin(theta) * std::cos(phi);
    };
  }
};
}  // namespace

int main() {
  constexpr int lmax = 2;
  Grid grid(lmax, 2);
  const auto quadrature =
      GaussQuad::GaussLobattoLegendreQuadrature1D<double>(3);
  GeneralEarthModels::spherical_1D one_d("modeldata/PREM500", quadrature,
                                          0.5, 1.2);
  MatrixReplacement<Complex> spherical(one_d, grid);
  const auto no_tensor = [](std::size_t, std::size_t, std::size_t) -> const Eigen::Matrix3cd & {
    static const Eigen::Matrix3cd a = [] {
      Eigen::Matrix3cd value;
      value << 0.0, 0.0, -1.0,
               0.0, 1.0, 0.0,
              -1.0, 0.0, 0.0;
      return value;
    }();
    return a;
  };

  Eigen::VectorXcd x(spherical.cols());
  Eigen::VectorXd real_x(spherical.cols());
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    x[i] = Complex(0.125 * (i + 1), -0.0625 * (i + 2));
    real_x[i] = 0.03125 * (i + 3);
  }
  const auto expected = direct_weak_action(one_d, grid, no_tensor, x, false, true);
  const auto spherical_actual = spherical.Apply(x);
  require_close(spherical_actual, expected, "spherical full-complex weak action");
  (void)spherical.Apply(real_x);
  require_close(spherical.Apply(x), spherical_actual,
                "spherical repeated-call buffer reset");
  const auto expected_real = direct_weak_action(
      one_d, grid, no_tensor, real_x.cast<Complex>(), false, true);
  require_close(spherical * real_x, expected_real, "spherical real-vector action");
  const Complex alpha(0.375, -0.25);
  Eigen::VectorXcd destination = Eigen::VectorXcd::Constant(x.size(), Complex(2.0, -1.0));
  const Eigen::VectorXcd initial = destination;
  destination.noalias() += alpha * (spherical * x);
  require_close(destination, initial + alpha * expected,
                "spherical Eigen alpha accumulation");
  const auto spherical_without_boundary = spherical.WithBoundary(false);
  const auto expected_without_boundary =
      direct_weak_action(one_d, grid, no_tensor, x, false, false);
  require_close(spherical_without_boundary * x, expected_without_boundary,
                "spherical boundary-disabled action");
  Eigen::VectorXcd boundary_delta = expected - expected_without_boundary;
  Eigen::VectorXcd expected_boundary = Eigen::VectorXcd::Zero(x.size());
  const std::size_t nh = (lmax + 1) * (lmax + 1);
  const std::size_t last = one_d.Node_Information().NumberOfGlobalRadialNodes() - 1;
  for (int degree = 0; degree <= lmax; ++degree)
    for (int order = -degree; order <= degree; ++order) {
      const auto index = last * nh + scalar_index(degree, order);
      expected_boundary[index] = -(degree + 1.0) *
          one_d.Node_Information().OuterRadius() * x[index];
    }
  require_close(boundary_delta, expected_boundary, "spherical boundary delta");

  auto base = GeneralEarthModels::SimpleModels::spherical_model::HomogeneousSphere(
      1.0, 1.0, 6371000.0, 3600.0, 5.972e24);
  SmoothMapping mapping;
  GeneralEarthModels::Density3D mapped(
      base, PlanetaryModel::TomographyZeroModel(), mapping,
      2, lmax, 0.5, 1.2);
  MatrixReplacement3D<Complex> mapped_operator(mapped);
  Eigen::VectorXcd mapped_x(mapped_operator.cols());
  Eigen::VectorXd mapped_real_x(mapped_operator.cols());
  for (Eigen::Index i = 0; i < mapped_x.size(); ++i) {
    mapped_x[i] = Complex(0.0625 * (i + 2), -0.03125 * (i + 1));
    mapped_real_x[i] = 0.015625 * (i + 4);
  }
  const auto mapped_tensor = [&mapped](std::size_t e, std::size_t j,
                                       std::size_t point) -> const Eigen::Matrix3cd & {
    return mapped.LaplaceTensorRef()[e][j][point];
  };
  const auto mapped_expected = direct_weak_action(
      mapped, grid, mapped_tensor, mapped_x, true, true);
  const auto mapped_actual = mapped_operator.Apply(mapped_x);
  require_close(mapped_actual, mapped_expected,
                "mapped full-complex direct angular quadrature");
  (void)mapped_operator.Apply(mapped_real_x);
  require_close(mapped_operator.Apply(mapped_x), mapped_actual,
                "mapped repeated-call buffer reset");
  const auto mapped_expected_real = direct_weak_action(
      mapped, grid, mapped_tensor, mapped_real_x.cast<Complex>(), true, true);
  require_close(mapped_operator * mapped_real_x, mapped_expected_real,
                "mapped real-vector action");
  Eigen::VectorXcd mapped_destination =
      Eigen::VectorXcd::Constant(mapped_x.size(), Complex(-0.5, 0.25));
  const Eigen::VectorXcd mapped_initial = mapped_destination;
  mapped_destination.noalias() += alpha * (mapped_operator * mapped_x);
  require_close(mapped_destination, mapped_initial + alpha * mapped_expected,
                "mapped Eigen alpha accumulation");
  const auto mapped_without_boundary = mapped_operator.WithBoundary(false);
  const auto mapped_expected_without_boundary = direct_weak_action(
      mapped, grid, mapped_tensor, mapped_x, true, false);
  require_close(mapped_without_boundary * mapped_x,
                mapped_expected_without_boundary,
                "mapped boundary-disabled action");
  Eigen::VectorXcd mapped_boundary = Eigen::VectorXcd::Zero(mapped_x.size());
  const auto mapped_mesh = mapped.Node_Information();
  const std::size_t mapped_nh = (lmax + 1) * (lmax + 1);
  const std::size_t mapped_last = mapped_mesh.NumberOfGlobalRadialNodes() - 1;
  for (int degree = 0; degree <= lmax; ++degree)
    for (int order = -degree; order <= degree; ++order) {
      const auto index = mapped_last * mapped_nh + scalar_index(degree, order);
      mapped_boundary[index] = -(degree + 1.0) * mapped_mesh.OuterRadius() *
                               mapped_x[index];
    }
  require_close(mapped_expected - mapped_expected_without_boundary,
                mapped_boundary, "mapped boundary delta");

  GeneralEarthModels::MappingPerturbation perturbation(mapped, mapping);
  MatrixReplacement3D<Complex> perturbation_operator(mapped, perturbation, false);
  const auto perturbation_tensor = [&perturbation](std::size_t e, std::size_t j,
                                                   std::size_t point) -> const Eigen::Matrix3cd & {
    return perturbation.delta_laplace_tensor()[e][j][point];
  };
  const auto perturbation_expected = direct_weak_action(
      mapped, grid, perturbation_tensor, mapped_x, true, false);
  require_close(perturbation_operator * mapped_x, perturbation_expected,
                "mapped perturbation tensor action");

  std::cout << "Stage 3C.2 direct weak actions, mapped tensor, alpha, real/complex, "
               "and explicit boundary checks passed\n";
}
