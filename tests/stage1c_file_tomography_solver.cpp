#include <gplspec/All>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>

int main() {
  using Complex = std::complex<double>;
  using Model = GeneralEarthModels::Density3D;
  using Operator = MatrixReplacement3D<Complex>;
  using Preconditioner = Eigen::SphericalGeometryPreconditioner<Complex>;
  using Solver = Eigen::ConjugateGradient<Operator, Eigen::Lower | Eigen::Upper,
                                          Preconditioner>;
  constexpr double tolerance = 1e-6;

  auto model = Model::SphericalThreeDimensionalPlanetFromFile(
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

  struct SmoothMapping {
    auto RadialMapping(int) const {
      return [](double r, double theta, double phi) {
        return 0.1 * r * (1.0 - r) * std::sin(theta) * std::cos(phi);
      };
    }
  } mapping;
  auto base = GeneralEarthModels::SimpleModels::spherical_model::HomogeneousSphere(
      1.0, 1.0, 6371000.0, 3600.0, 5.972e24);
  Model solver_model(base, PlanetaryModel::TomographyZeroModel(), mapping,
                     2, 2, 0.5, 1.2);
  Operator solver_op(solver_model);
  Eigen::SparseMatrix<Complex> preconditioner_matrix =
      -solver_model.SpectralElementInformation().fullmatrix<Complex>(
          0, solver_model.GSH_Grid().MaxDegree());
  preconditioner_matrix.makeCompressed();
  const Eigen::VectorXcd force = Gravity_Tools::FindForce(solver_model);
  Solver solver;
  solver.compute(solver_op);
  solver.preconditioner().addmatrix(preconditioner_matrix);
  solver.setTolerance(tolerance);
  const Eigen::VectorXcd guess = solver.preconditioner().solve(force);
  const Eigen::VectorXcd solution = solver.solveWithGuess(force, guess);
  const double residual_norm = (force - solver_op * solution).norm();
  const double relative_residual = residual_norm / force.norm();
  if (solver.info() != Eigen::Success || !solution.allFinite() ||
      solver.iterations() < 2 || !std::isfinite(solver.error()) ||
      !std::isfinite(residual_norm) || !std::isfinite(relative_residual) ||
      force.norm() == 0.0 || relative_residual > tolerance) {
    std::cerr << "solver info=" << solver.info()
              << " iterations=" << solver.iterations()
              << " reported_error=" << solver.error()
              << " absolute_residual=" << residual_norm
              << " independent_relative_residual=" << relative_residual << '\n';
    throw std::runtime_error("mapped gravity solve failed its diagnostics");
  }

  const auto public_solution =
      Gravity_Tools::FindGravitationalPotential(solver_model, tolerance);
  const auto expected_public_solution =
      solver_model.SingleEigenVectorToGeneralFormat(solution);
  std::size_t finite_values = 0;
  double public_difference = 0.0;
  double public_scale = 0.0;
  if (public_solution.size() != expected_public_solution.size())
    throw std::runtime_error("public gravity workflow returned an unexpected element count");
  for (std::size_t element = 0; element < public_solution.size(); ++element) {
    if (public_solution[element].size() != expected_public_solution[element].size())
      throw std::runtime_error("public gravity workflow returned an unexpected node count");
    for (std::size_t node = 0; node < public_solution[element].size(); ++node) {
      if (public_solution[element][node].size() != expected_public_solution[element][node].size())
        throw std::runtime_error("public gravity workflow returned an unexpected coefficient count");
      for (std::size_t coefficient = 0;
           coefficient < public_solution[element][node].size(); ++coefficient) {
        const Complex value = public_solution[element][node][coefficient];
        const Complex expected = expected_public_solution[element][node][coefficient];
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()) ||
            !std::isfinite(expected.real()) || !std::isfinite(expected.imag()))
          throw std::runtime_error("public gravity solution contains a nonfinite coefficient");
        ++finite_values;
        public_difference = std::max(public_difference, std::abs(value - expected));
        public_scale = std::max(public_scale, std::abs(expected));
      }
    }
  }
  if (finite_values == 0 || public_difference > 1e-12 * std::max(1.0, public_scale))
    throw std::runtime_error("public gravity solution differs from diagnosed iterative solution");

  std::cout << "rows=" << op.rows() << " iterations=" << solver.iterations()
            << " solver_error=" << solver.error()
            << " absolute_residual=" << residual_norm
            << " independent_relative_residual=" << relative_residual
            << " public_difference=" << public_difference
            << " public_finite_values=" << finite_values << '\n';
}
