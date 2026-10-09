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
  const auto &grid = model.GSH_Grid();
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
  using SphericalOperator = MatrixReplacement<Complex>;
  using SphericalGrid = SphericalOperator::Grid;
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

  const auto preconditioner_matrix =
      Gravity_Tools::detail::SphericalPreconditionerMatrix(
          model, model.GSH_Grid().MaxDegree());
  Gravity_Tools::detail::PreparedSolver<Solver, Operator> solver(
      op, preconditioner_matrix);
  bool invalid_tolerance_rejected = false;
  try {
    (void)solver.Solve(force, 0.0, true);
  } catch (const std::invalid_argument &) {
    invalid_tolerance_rejected = true;
  }
  if (!invalid_tolerance_rejected)
    throw std::runtime_error("solver accepted a non-positive tolerance");
  Preconditioner guess_preconditioner;
  guess_preconditioner.compute(op);
  guess_preconditioner.addmatrix(preconditioner_matrix);
  const Eigen::VectorXcd guess = guess_preconditioner.solve(force);
  const double guess_relative_residual = (force - op * guess).norm() / force.norm();
  const auto solve_result = solver.Solve(force, tolerance, true);
  const auto &solve_diagnostics = solve_result.final_solve();
  const Eigen::VectorXcd &solution = solve_result.solution;
  const double update_relative_norm = (solution - guess).norm() / guess.norm();
  const double independently_measured_residual = (force - op * solution).norm();
  const double absolute_residual = solve_diagnostics.residual_norm;
  const double relative_residual = solve_diagnostics.relative_residual;
  if (solve_diagnostics.status != Eigen::Success || !solution.allFinite() ||
      !std::isfinite(solve_diagnostics.recursive_error) ||
      !std::isfinite(absolute_residual) || !std::isfinite(relative_residual) ||
      !std::isfinite(guess_relative_residual) ||
      std::abs(independently_measured_residual - absolute_residual) >
          1e-13 * std::max(1.0, independently_measured_residual) ||
      relative_residual > tolerance) {
    std::cerr << "file-backed solver info=" << solve_diagnostics.status
              << " iterations=" << solve_diagnostics.iterations
              << " requested_tolerance=" << solve_diagnostics.requested_tolerance
              << " initial_guess_relative_residual=" << guess_relative_residual
              << " update_relative_norm=" << update_relative_norm
              << " recursive_error=" << solve_diagnostics.recursive_error
              << " absolute_residual=" << absolute_residual
              << " independent_relative_residual=" << relative_residual << '\n';
    throw std::runtime_error("file-backed gravity solve failed its diagnostics");
  }

  const double stricter_tolerance =
      std::min(tolerance, 0.1 * guess_relative_residual);
  const auto updated_result =
      solver.Solve(force, stricter_tolerance, true);
  const auto &updated_diagnostics = updated_result.final_solve();
  const double strict_update_norm =
      (updated_result.solution - guess).norm() / guess.norm();
  if (updated_diagnostics.status != Eigen::Success || !(strict_update_norm > 0.0) ||
      !(updated_diagnostics.residual_norm < absolute_residual))
    throw std::runtime_error("tight solve did not improve the preconditioned guess");

  const Eigen::VectorXcd zero_rhs = Eigen::VectorXcd::Zero(force.size());
  const auto zero_result = solver.Solve(zero_rhs, tolerance, true);
  const auto &zero_diagnostics = zero_result.final_solve();
  if (zero_diagnostics.status != Eigen::Success ||
      zero_diagnostics.residual_norm != 0.0 || zero_diagnostics.relative_residual != 0.0 ||
      zero_result.solution.norm() != 0.0)
    throw std::runtime_error("zero-RHS solve diagnostics are not well-defined");

  const auto public_result = Gravity_Tools::SolvePotential(model, tolerance);
  const auto public_solution =
      Gravity_Tools::ToHarmonicField(model, public_result);
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

  const auto public_diagnostic =
      Gravity_Tools::SolvePotential(model, tolerance);
  const auto &public_solve = public_diagnostic.final_solve();
  if (public_solve.status != Eigen::Success ||
      !std::isfinite(public_solve.recursive_error) ||
      public_solve.relative_residual > tolerance ||
      public_diagnostic.solves.size() != 1 ||
      (public_diagnostic.solution - solution).norm() >
          1e-12 * std::max(1.0, solution.norm()))
    throw std::runtime_error("public solve diagnostics differ from prepared solve");

  auto failed_public_result = public_diagnostic;
  failed_public_result.solves.front().status = Eigen::NumericalIssue;
  bool failed_result_rejected = false;
  try {
    (void)Gravity_Tools::ToHarmonicField(model, failed_public_result);
  } catch (const std::runtime_error &) {
    failed_result_rejected = true;
  }
  auto malformed_public_result = public_diagnostic;
  malformed_public_result.solution.conservativeResize(
      malformed_public_result.solution.size() - 1);
  bool malformed_result_rejected = false;
  try {
    (void)Gravity_Tools::ToHarmonicField(model, malformed_public_result);
  } catch (const std::invalid_argument &) {
    malformed_result_rejected = true;
  }
  if (!failed_result_rejected || !malformed_result_rejected)
    throw std::runtime_error("harmonic output conversion accepted an invalid solve result");

  Eigen::SphericalGeometryPreconditioner<Complex> identity_preconditioner;
  identity_preconditioner.compute(op);
  const Eigen::VectorXcd identity_probe = Eigen::VectorXcd::LinSpaced(
      op.cols(), Complex(1.0, 0.5), Complex(2.0, -0.25));
  if ((identity_preconditioner.solve(identity_probe) - identity_probe).norm() != 0.0 ||
      identity_preconditioner.info() != Eigen::Success)
    throw std::runtime_error("unconfigured preconditioner must apply identity");

  const auto radial_quadrature =
      GaussQuad::GaussLobattoLegendreQuadrature1D<double>(4);
  GeneralEarthModels::spherical_1D radial_model(
      "modeldata/PREM500", radial_quadrature, 0.5, 1.2);
  SphericalGrid spherical_grid(2, 2);
  SphericalOperator spherical_operator(radial_model, spherical_grid);
  const auto spherical_force = Gravity_Tools::FindForce(radial_model, 2);
  using SphericalSolver = Eigen::BiCGSTAB<
      SphericalOperator, Eigen::SphericalGeometryPreconditioner<Complex>>;
  const auto spherical_preconditioner =
      Gravity_Tools::detail::SphericalPreconditionerMatrix(radial_model, 2);
  Gravity_Tools::detail::PreparedSolver<SphericalSolver, SphericalOperator>
      spherical_solver(spherical_operator, spherical_preconditioner);
  Eigen::SphericalGeometryPreconditioner<Complex> spherical_guess_preconditioner;
  spherical_guess_preconditioner.compute(spherical_operator);
  spherical_guess_preconditioner.addmatrix(spherical_preconditioner);
  const Eigen::VectorXcd spherical_guess =
      spherical_guess_preconditioner.solve(spherical_force);
  const auto spherical_result =
      spherical_solver.Solve(spherical_force, tolerance, false);
  const auto &spherical_diagnostics = spherical_result.final_solve();
  const double spherical_update_norm =
      (spherical_result.solution - spherical_guess).norm();
  const double spherical_true_residual =
      (spherical_force - spherical_operator * spherical_result.solution).norm();
  const auto public_spherical = Gravity_Tools::SolvePotential(
      radial_model, spherical_grid, tolerance);
  if (spherical_diagnostics.status != Eigen::Success ||
      !std::isfinite(spherical_diagnostics.recursive_error) ||
      !std::isfinite(spherical_diagnostics.relative_residual) ||
      !(spherical_update_norm > 0.0) ||
      std::abs(spherical_true_residual - spherical_diagnostics.residual_norm) >
          1e-13 * std::max(1.0, spherical_true_residual) ||
      spherical_diagnostics.relative_residual > tolerance ||
      (public_spherical.solution - spherical_result.solution).norm() >
          1e-12 * std::max(1.0, spherical_result.solution.norm()))
    throw std::runtime_error("spherical BiCGSTAB solve contract failed");

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
  const Eigen::VectorXcd mapped_force = Gravity_Tools::FindForce(mapped_model);
  const auto mapped_preconditioner_matrix =
      Gravity_Tools::detail::SphericalPreconditionerMatrix(
          mapped_model, mapped_model.GSH_Grid().MaxDegree());
  Gravity_Tools::detail::PreparedSolver<Solver, Operator> mapped_solver(
      mapped_op, mapped_preconditioner_matrix);
  Preconditioner mapped_guess_preconditioner;
  mapped_guess_preconditioner.compute(mapped_op);
  mapped_guess_preconditioner.addmatrix(mapped_preconditioner_matrix);
  const Eigen::VectorXcd mapped_guess =
      mapped_guess_preconditioner.solve(mapped_force);
  const double mapped_guess_relative_residual =
      (mapped_force - mapped_op * mapped_guess).norm() / mapped_force.norm();
  const auto mapped_result = mapped_solver.Solve(mapped_force, tolerance, true);
  const auto &mapped_diagnostics = mapped_result.final_solve();
  const Eigen::VectorXcd &mapped_solution = mapped_result.solution;
  const double mapped_update_relative_norm =
      (mapped_solution - mapped_guess).norm() / mapped_guess.norm();
  const double mapped_independently_measured_residual =
      (mapped_force - mapped_op * mapped_solution).norm();
  const double mapped_absolute_residual = mapped_diagnostics.residual_norm;
  const double mapped_relative_residual = mapped_diagnostics.relative_residual;
  if (mapped_diagnostics.status != Eigen::Success || !mapped_solution.allFinite() ||
      mapped_diagnostics.iterations < 2 || !std::isfinite(mapped_diagnostics.recursive_error) ||
      !std::isfinite(mapped_absolute_residual) ||
      !std::isfinite(mapped_relative_residual) ||
      !std::isfinite(mapped_guess_relative_residual) ||
      !(mapped_update_relative_norm > 0.0) ||
      std::abs(mapped_independently_measured_residual - mapped_absolute_residual) >
          1e-13 * std::max(1.0, mapped_independently_measured_residual) ||
      mapped_force.norm() == 0.0 ||
      mapped_relative_residual > tolerance) {
    std::cerr << "mapped solver info=" << mapped_diagnostics.status
              << " iterations=" << mapped_diagnostics.iterations
              << " requested_tolerance=" << mapped_diagnostics.requested_tolerance
              << " initial_guess_relative_residual="
              << mapped_guess_relative_residual
              << " update_relative_norm=" << mapped_update_relative_norm
              << " recursive_error=" << mapped_diagnostics.recursive_error
              << " absolute_residual=" << mapped_absolute_residual
              << " independent_relative_residual=" << mapped_relative_residual << '\n';
    throw std::runtime_error("mapped gravity solve failed its diagnostics");
  }
  const auto mapped_public_solution =
      Gravity_Tools::ToHarmonicField(
          mapped_model,
          Gravity_Tools::SolvePotential(mapped_model, tolerance));
  const auto mapped_public_diagnostic =
      Gravity_Tools::SolvePotential(mapped_model, tolerance);
  const auto &mapped_public_solve = mapped_public_diagnostic.final_solve();
  if (mapped_public_solve.status != Eigen::Success ||
      !std::isfinite(mapped_public_solve.recursive_error) ||
      mapped_public_solve.relative_residual > tolerance ||
      (mapped_public_diagnostic.solution - mapped_solution).norm() >
          1e-12 * std::max(1.0, mapped_solution.norm()))
    throw std::runtime_error(
        "public mapped solve diagnostics differ from prepared solve");
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
            << " iterations=" << solve_diagnostics.iterations
            << " recursive_error=" << solve_diagnostics.recursive_error
            << " initial_guess_relative_residual=" << guess_relative_residual
            << " update_relative_norm=" << update_relative_norm
            << " tight_update_relative_norm=" << strict_update_norm
            << " absolute_residual=" << absolute_residual
            << " independent_relative_residual=" << relative_residual
            << " public_difference=" << public_difference
            << " public_finite_values=" << finite_values
            << " mapped_rows=" << mapped_op.rows()
            << " mapped_iterations=" << mapped_diagnostics.iterations
            << " mapped_recursive_error=" << mapped_diagnostics.recursive_error
            << " mapped_initial_guess_relative_residual="
            << mapped_guess_relative_residual
            << " mapped_update_relative_norm=" << mapped_update_relative_norm
            << " mapped_independent_relative_residual=" << mapped_relative_residual
            << " mapped_public_difference=" << mapped_public_difference
            << " mapped_public_finite_values=" << mapped_finite_values << '\n';
}
