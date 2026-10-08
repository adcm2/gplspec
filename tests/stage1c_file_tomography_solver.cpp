#include <gplspec/All>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
using Complex = std::complex<double>;

// Independently integrate the l=0,m=0 and l=1,2,m=0 source coefficients
// over the exact model grid, then apply the radial Gauss rule and element map.
// The analytic harmonics avoid calling the transform used by FindForce.
template <class Model>
Eigen::VectorXcd independent_axisymmetric_force(const Model &model) {
  const auto grid = model.GSH_Grid();
  const auto q = model.q();
  constexpr double pi = std::numbers::pi_v<double>;
  const int lmax = grid.MaxDegree();
  const int stride = (lmax + 1) * (lmax + 1);
  const int nnode = q.N() - 1;
  const int nelem = model.Node_Information().NumberOfElements();
  const int radial_size = nelem * nnode + 1;
  Eigen::VectorXcd result = Eigen::VectorXcd::Zero(stride * radial_size);

  auto add_if_supported = [&](int element, int node, int l, int m,
                              Complex value) {
    if (l <= lmax) {
      const int radial = node + nnode * element;
      const int harmonic = l * l + l + m;
      result(radial * stride + harmonic) += value;
    }
  };

  for (int element = 0; element < nelem; ++element) {
    const double width = model.Node_Information().ElementWidth(element);
    const double lower = model.Node_Information().ElementLowerRadius(element);
    const double upper = model.Node_Information().ElementUpperRadius(element);
    for (int node = 0; node < q.N(); ++node) {
      const auto density = model.DensityAtRadialNode(element, node);
      Complex angular00{0.0, 0.0};
      Complex angular10{0.0, 0.0};
      Complex angular20{0.0, 0.0};
      int theta_index = 0;
      for (const double theta : grid.CoLatitudes()) {
        const double cosine = std::cos(theta);
        const double y00 = 1.0 / std::sqrt(4.0 * pi);
        const double y10 = std::sqrt(3.0 / (4.0 * pi)) * cosine;
        const double y20 = std::sqrt(5.0 / (16.0 * pi)) *
                           (3.0 * cosine * cosine - 1.0);
        const double dmu = grid.CoLatitudeWeights()[theta_index];
        int phi_index = 0;
        for (const double unused_phi : grid.Longitudes()) {
          (void)unused_phi;
          const double dphi = grid.LongitudeWeights()[phi_index];
          const double sample = density[theta_index * grid.NumberOfLongitudes() +
                                        phi_index];
          const double weight = dmu * dphi;
          angular00 += sample * y00 * weight;
          angular10 += sample * y10 * weight;
          angular20 += sample * y20 * weight;
          ++phi_index;
        }
        ++theta_index;
      }
      const double radius =
          0.5 * (width * q.X(node) + upper + lower);
      const double radial_weight =
          4.0 * pi * model.GravitationalConstant() * q.W(node) *
          (0.5 * width) * radius * radius;
      add_if_supported(element, node, 0, 0, angular00 * radial_weight);
      add_if_supported(element, node, 1, 0, angular10 * radial_weight);
      add_if_supported(element, node, 2, 0, angular20 * radial_weight);
    }
  }
  return result;
}
}  // namespace

int main() {
  using Model = GeneralEarthModels::Density3D;
  using Operator = MatrixReplacement3D<Complex>;
  using Preconditioner = Eigen::SphericalGeometryPreconditioner<Complex>;
  using Solver = Eigen::ConjugateGradient<Operator, Eigen::Lower | Eigen::Upper,
                                          Preconditioner>;
  constexpr double tolerance = 1e-6;

  Model model = Model::SphericalThreeDimensionalPlanetFromFile(
      "modeldata/PREM500", "modeldata/S40RTS_dvs.nc", 2, 2, 0.5, 1.2);
  Operator op(model);
  if (op.rows() <= 0 || op.rows() != op.cols())
    throw std::runtime_error("file-backed tomography produced an invalid operator size");

  Eigen::VectorXcd action_input(op.cols());
  for (Eigen::Index i = 0; i < action_input.size(); ++i)
    action_input(i) = Complex(0.125 * (i + 1), -0.0625 * (i + 2));
  const Eigen::VectorXcd action = op * action_input;
  if (!action.allFinite() || action.norm() == 0.0)
    throw std::runtime_error("file-backed tomography complex action is not finite and nonzero");

  const Eigen::VectorXcd force = Gravity_Tools::FindForce(model);
  const Eigen::VectorXcd independent_force = independent_axisymmetric_force(model);
  if (force.size() != independent_force.size() || force.norm() == 0.0 ||
      !force.allFinite() || !independent_force.allFinite())
    throw std::runtime_error("file-backed source vector is invalid");

  double max_source_error = 0.0;
  double max_source_scale = 0.0;
  for (Eigen::Index radial = 0; radial < force.size() / 9; ++radial) {
    for (const int harmonic : {0, 2, 6}) {
      const Eigen::Index index = radial * 9 + harmonic;
      max_source_error = std::max(max_source_error,
                                  std::abs(force(index) - independent_force(index)));
      max_source_scale = std::max(max_source_scale,
                                  std::abs(independent_force(index)));
    }
  }
  if (max_source_error > 2e-11 * std::max(1.0, max_source_scale)) {
    std::cerr << "independent file-backed source quadrature mismatch: max_abs="
              << max_source_error << " scale=" << max_source_scale << '\n';
    throw std::runtime_error("file-backed source quadrature mismatch");
  }

  Eigen::SparseMatrix<Complex> preconditioner_matrix =
      -model.SpectralElementInformation().fullmatrix<Complex>(
          0, model.GSH_Grid().MaxDegree());
  preconditioner_matrix.makeCompressed();
  Solver solver;
  solver.compute(op);
  solver.preconditioner().addmatrix(preconditioner_matrix);
  solver.setTolerance(tolerance);
  const Eigen::VectorXcd guess = solver.preconditioner().solve(force);
  const Eigen::VectorXcd solution = solver.solveWithGuess(force, guess);
  const double absolute_residual = (force - op * solution).norm();
  const double relative_residual = absolute_residual / force.norm();
  if (solver.info() != Eigen::Success || !solution.allFinite() ||
      !std::isfinite(solver.error()) ||
      !std::isfinite(absolute_residual) || !std::isfinite(relative_residual) ||
      relative_residual > tolerance) {
    std::cerr << "file-backed solver info=" << solver.info()
              << " iterations=" << solver.iterations()
              << " reported_error=" << solver.error()
              << " absolute_residual=" << absolute_residual
              << " independent_relative_residual=" << relative_residual << '\n';
    throw std::runtime_error("file-backed gravity solve failed its diagnostics");
  }

  const auto public_solution = Gravity_Tools::FindGravitationalPotential(model, tolerance);
  const auto expected_public_solution = model.SingleEigenVectorToGeneralFormat(solution);
  if (public_solution.size() != expected_public_solution.size())
    throw std::runtime_error("public file-backed gravity workflow returned an unexpected element count");
  double public_difference = 0.0;
  double public_scale = 0.0;
  std::size_t finite_values = 0;
  for (std::size_t element = 0; element < public_solution.size(); ++element) {
    if (public_solution[element].size() != expected_public_solution[element].size())
      throw std::runtime_error("public file-backed gravity workflow returned an unexpected node count");
    for (std::size_t node = 0; node < public_solution[element].size(); ++node) {
      if (public_solution[element][node].size() != expected_public_solution[element][node].size())
        throw std::runtime_error("public file-backed gravity workflow returned an unexpected coefficient count");
      for (std::size_t coefficient = 0;
           coefficient < public_solution[element][node].size(); ++coefficient) {
        const Complex actual = public_solution[element][node][coefficient];
        const Complex expected = expected_public_solution[element][node][coefficient];
        if (!std::isfinite(actual.real()) || !std::isfinite(actual.imag()) ||
            !std::isfinite(expected.real()) || !std::isfinite(expected.imag()))
          throw std::runtime_error("public file-backed solution contains a nonfinite coefficient");
        ++finite_values;
        public_difference = std::max(public_difference, std::abs(actual - expected));
        public_scale = std::max(public_scale, std::abs(expected));
      }
    }
  }
  if (finite_values == 0 || public_difference > 1e-12 * std::max(1.0, public_scale))
    throw std::runtime_error("public file-backed solution differs from diagnosed iterative solution");

  struct SmoothMapping {
    auto RadialMapping(int) const {
      return [](double r, double theta, double phi) {
        return 0.1 * r * (1.0 - r) * std::sin(theta) * std::cos(phi);
      };
    }
  } mapping;
  auto base = GeneralEarthModels::SimpleModels::spherical_model::HomogeneousSphere(
      1.0, 1.0, 6371000.0, 3600.0, 5.972e24);
  Model mapped_model(base, PlanetaryModel::TomographyZeroModel(), mapping,
                     2, 2, 0.5, 1.2);
  Operator mapped_op(mapped_model);
  Eigen::SparseMatrix<Complex> mapped_preconditioner_matrix =
      -mapped_model.SpectralElementInformation().fullmatrix<Complex>(
          0, mapped_model.GSH_Grid().MaxDegree());
  mapped_preconditioner_matrix.makeCompressed();
  const Eigen::VectorXcd mapped_force = Gravity_Tools::FindForce(mapped_model);
  Solver mapped_solver;
  mapped_solver.compute(mapped_op);
  mapped_solver.preconditioner().addmatrix(mapped_preconditioner_matrix);
  mapped_solver.setTolerance(tolerance);
  const Eigen::VectorXcd mapped_guess =
      mapped_solver.preconditioner().solve(mapped_force);
  const Eigen::VectorXcd mapped_solution =
      mapped_solver.solveWithGuess(mapped_force, mapped_guess);
  const double mapped_absolute_residual =
      (mapped_force - mapped_op * mapped_solution).norm();
  const double mapped_relative_residual =
      mapped_absolute_residual / mapped_force.norm();
  if (mapped_solver.info() != Eigen::Success || !mapped_solution.allFinite() ||
      mapped_solver.iterations() < 2 || !std::isfinite(mapped_solver.error()) ||
      !std::isfinite(mapped_absolute_residual) ||
      !std::isfinite(mapped_relative_residual) || mapped_force.norm() == 0.0 ||
      mapped_relative_residual > tolerance) {
    std::cerr << "mapped solver info=" << mapped_solver.info()
              << " iterations=" << mapped_solver.iterations()
              << " reported_error=" << mapped_solver.error()
              << " absolute_residual=" << mapped_absolute_residual
              << " independent_relative_residual=" << mapped_relative_residual << '\n';
    throw std::runtime_error("mapped gravity solve failed its diagnostics");
  }
  const auto mapped_public_solution =
      Gravity_Tools::FindGravitationalPotential(mapped_model, tolerance);
  const auto expected_mapped_solution =
      mapped_model.SingleEigenVectorToGeneralFormat(mapped_solution);
  if (mapped_public_solution.size() != expected_mapped_solution.size())
    throw std::runtime_error("public mapped gravity workflow returned an unexpected element count");
  double mapped_public_difference = 0.0;
  double mapped_public_scale = 0.0;
  std::size_t mapped_finite_values = 0;
  for (std::size_t element = 0; element < mapped_public_solution.size(); ++element) {
    if (mapped_public_solution[element].size() != expected_mapped_solution[element].size())
      throw std::runtime_error("public mapped gravity workflow returned an unexpected node count");
    for (std::size_t node = 0; node < mapped_public_solution[element].size(); ++node) {
      if (mapped_public_solution[element][node].size() != expected_mapped_solution[element][node].size())
        throw std::runtime_error("public mapped gravity workflow returned an unexpected coefficient count");
      for (std::size_t coefficient = 0;
           coefficient < mapped_public_solution[element][node].size(); ++coefficient) {
        const Complex actual = mapped_public_solution[element][node][coefficient];
        const Complex expected = expected_mapped_solution[element][node][coefficient];
        if (!std::isfinite(actual.real()) || !std::isfinite(actual.imag()) ||
            !std::isfinite(expected.real()) || !std::isfinite(expected.imag()))
          throw std::runtime_error("public mapped solution contains a nonfinite coefficient");
        ++mapped_finite_values;
        mapped_public_difference = std::max(mapped_public_difference,
                                            std::abs(actual - expected));
        mapped_public_scale = std::max(mapped_public_scale, std::abs(expected));
      }
    }
  }
  if (mapped_finite_values == 0 ||
      mapped_public_difference > 1e-12 * std::max(1.0, mapped_public_scale))
    throw std::runtime_error("public mapped solution differs from diagnosed iterative solution");

  std::cout << "rows=" << op.rows() << " force_norm=" << force.norm()
            << " quadrature_max_abs_error=" << max_source_error
            << " iterations=" << solver.iterations()
            << " solver_error=" << solver.error()
            << " absolute_residual=" << absolute_residual
            << " independent_relative_residual=" << relative_residual
            << " public_difference=" << public_difference
            << " public_finite_values=" << finite_values
            << " mapped_rows=" << mapped_op.rows()
            << " mapped_iterations=" << mapped_solver.iterations()
            << " mapped_independent_relative_residual=" << mapped_relative_residual
            << " mapped_public_difference=" << mapped_public_difference
            << " mapped_public_finite_values=" << mapped_finite_values << '\n';
}
