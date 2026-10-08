#include <gplspec/All>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Complex = std::complex<double>;
using Model = GeneralEarthModels::Density3D;
using Operator = MatrixReplacement3D<Complex>;
using Preconditioner = Eigen::SphericalGeometryPreconditioner<Complex>;
using Solver = Eigen::ConjugateGradient<Operator, Eigen::Lower | Eigen::Upper,
                                        Preconditioner>;
using Clock = std::chrono::steady_clock;

struct Trial {
  double total_seconds = 0.0;
  double setup_seconds = 0.0;
  double guess_seconds = 0.0;
  double iteration_seconds = 0.0;
  double solver_error = -1.0;
  double absolute_residual = -1.0;
  double relative_residual = -1.0;
  double checksum = 0.0;
  int solver_info = -1;
  int iterations = -1;
};

double elapsed_seconds(Clock::time_point start, Clock::time_point end) {
  return std::chrono::duration<double>(end - start).count();
}

double median(std::vector<double> values) {
  std::sort(values.begin(), values.end());
  const auto middle = values.size() / 2;
  return values.size() % 2 ? values[middle]
                            : 0.5 * (values[middle - 1] + values[middle]);
}

void emit_summary(const std::string &scenario, const std::string &metric,
                  const std::vector<double> &samples, double checksum,
                  double target, bool target_met) {
  std::vector<double> sorted = samples;
  std::sort(sorted.begin(), sorted.end());
  const double total = std::accumulate(samples.begin(), samples.end(), 0.0);
  const double mean = total / static_cast<double>(samples.size());
  std::cout << scenario << ',' << metric << ",batch," << samples.size() << ','
            << std::setprecision(17) << total << ',' << mean << ','
            << sorted.front() << ',' << median(samples) << ',' << sorted.back()
            << ",-1,-1," << target << ",-1,-1,-1,-1," << checksum << ','
            << target_met << ",0,0,0\n";
}

void emit_solve(const std::string &scenario, const std::string &metric,
                int index, const Trial &trial, double target) {
  std::cout << scenario << ',' << metric << ',' << index << ",1,"
            << std::setprecision(17) << trial.total_seconds << ','
            << trial.total_seconds << ',' << trial.total_seconds << ','
            << trial.total_seconds << ',' << trial.total_seconds << ','
            << trial.solver_info << ',' << trial.iterations << ',' << target
            << ',' << trial.solver_error << ',' << trial.absolute_residual << ','
            << trial.relative_residual << ','
            << (trial.relative_residual <= target ? 1 : 0) << ','
            << trial.checksum << ",-1,"
            << trial.setup_seconds << ',' << trial.guess_seconds << ','
            << trial.iteration_seconds << '\n';
}

struct HomogeneousMap {
  auto RadialMapping(int) const {
    return [](double r, double, double) { return 0.2 * r; };
  }
};
}  // namespace

int main(int argc, char **argv) {
  if (argc != 3)
    throw std::runtime_error("usage: phase3a_performance <mapped|file_backed> <batch_seconds>");
  const std::string scenario = argv[1];
  const double batch_seconds = std::stod(argv[2]);
  if ((scenario != "mapped" && scenario != "file_backed") ||
      !(batch_seconds > 0.0) || !std::isfinite(batch_seconds))
    throw std::runtime_error("invalid scenario or batch duration");

  const bool mapped = scenario == "mapped";
  const double tolerance = mapped ? 1e-12 : 1e-14;
  std::ostringstream library_diagnostics;
  std::streambuf *const original_stdout = std::cout.rdbuf(library_diagnostics.rdbuf());
  constexpr double length_norm = 6371000.0;
  constexpr double time_norm = 3600.0;
  constexpr double mass_norm = 5.972e24;
  auto build_model = [&]() -> std::unique_ptr<Model> {
    if (mapped) {
      auto base = GeneralEarthModels::SimpleModels::spherical_model::HomogeneousSphere(
          1.0, 1.0, length_norm, time_norm, mass_norm);
      return std::make_unique<Model>(base,
                                     PlanetaryModel::TomographyZeroModel(),
                                     HomogeneousMap{}, 5, 2, 0.1, 1.4);
    }
    return std::make_unique<Model>(Model::SphericalThreeDimensionalPlanetFromFile(
        "modeldata/prem.200", "modeldata/S40RTS_dvs.nc", 12, 20, 0.05, 1.2));
  };
  std::unique_ptr<Model> model;
  std::vector<double> setup_samples;
  double setup_total = 0.0;
  constexpr int max_setup_repeats = 1000;
  for (int i = 0; i < max_setup_repeats && setup_total < batch_seconds; ++i) {
    model.reset();
    const auto start = Clock::now();
    model = build_model();
    const double elapsed = elapsed_seconds(start, Clock::now());
    setup_samples.push_back(elapsed);
    setup_total += elapsed;
  }
  if (!model || setup_samples.empty())
    throw std::runtime_error("model setup did not produce a model");

  const int lmax = mapped ? 2 : 20;
  const int rows = static_cast<int>(model->Num_Elements() * model->Poly_Order() + 1) *
                   (lmax + 1) * (lmax + 1);
  std::cout.rdbuf(original_stdout);
  std::cout << "scenario,metric,index,repeats,total_seconds,mean_seconds,min_seconds,median_seconds,max_seconds,solver_info,iterations,requested_tolerance,solver_error,absolute_residual,relative_residual,independent_gate_pass,checksum,target_met,solver_setup_seconds,initial_guess_seconds,iteration_seconds\n";
  std::cout << scenario << ",fixture,0,1,0,0,0,0,0,-1,-1," << tolerance
            << ",-1,-1,-1,-1," << rows << ",1,0,0,0\n";
  emit_summary(scenario, "model_setup", setup_samples, rows, batch_seconds,
               setup_total >= batch_seconds);

  std::unique_ptr<Operator> op = std::make_unique<Operator>(*model);
  if (rows != op->rows())
    throw std::runtime_error("reported fixture row count differs from operator rows");
  std::vector<double> force_samples;
  double force_total = 0.0;
  Eigen::VectorXcd force;
  constexpr int max_force_repeats = 10000;
  for (int i = 0; i < max_force_repeats && force_total < batch_seconds; ++i) {
    const auto start = Clock::now();
    force = Gravity_Tools::FindForce(*model);
    const double elapsed = elapsed_seconds(start, Clock::now());
    force_samples.push_back(elapsed);
    force_total += elapsed;
  }
  if (force.size() != op->rows() || !force.allFinite() || force.norm() == 0.0)
    throw std::runtime_error("source assembly produced an invalid force vector");
  emit_summary(scenario, "source_assembly", force_samples, force.norm(),
               batch_seconds, force_total >= batch_seconds);

  Eigen::VectorXcd input(op->cols());
  for (Eigen::Index i = 0; i < input.size(); ++i)
    input(i) = Complex(std::sin(0.01 * (i + 1)), std::cos(0.02 * (i + 1)));
  const auto first_action_start = Clock::now();
  Eigen::VectorXcd action = *op * input;
  const double first_action_seconds = elapsed_seconds(first_action_start, Clock::now());
  if (!action.allFinite() || action.norm() == 0.0)
    throw std::runtime_error("first operator action is invalid");
  emit_summary(scenario, "operator_action_first", {first_action_seconds},
               action.norm(), batch_seconds, false);

  std::vector<double> action_samples;
  double action_total = 0.0;
  constexpr int max_action_repeats = 10000;
  for (int i = 0; i < max_action_repeats && action_total < batch_seconds; ++i) {
    const auto start = Clock::now();
    action = *op * input;
    const double elapsed = elapsed_seconds(start, Clock::now());
    action_samples.push_back(elapsed);
    action_total += elapsed;
  }
  if (!action.allFinite() || action.norm() == 0.0)
    throw std::runtime_error("repeated operator action is invalid");
  emit_summary(scenario, "operator_action", action_samples, action.norm(),
               batch_seconds, action_total >= batch_seconds);

  std::vector<double> preconditioner_samples;
  double preconditioner_total = 0.0;
  Eigen::SparseMatrix<Complex> preconditioner_matrix;
  constexpr int max_preconditioner_repeats = 10000;
  for (int i = 0; i < max_preconditioner_repeats &&
                  preconditioner_total < batch_seconds; ++i) {
    const auto start = Clock::now();
    preconditioner_matrix =
        -model->SpectralElementInformation().fullmatrix<Complex>(0, lmax);
    preconditioner_matrix.makeCompressed();
    const double elapsed = elapsed_seconds(start, Clock::now());
    preconditioner_samples.push_back(elapsed);
    preconditioner_total += elapsed;
  }
  if (preconditioner_matrix.rows() != op->rows() ||
      preconditioner_matrix.cols() != op->cols())
    throw std::runtime_error("preconditioner matrix has incorrect dimensions");
  emit_summary(scenario, "preconditioner_matrix_setup", preconditioner_samples,
               preconditioner_matrix.norm(), batch_seconds,
               preconditioner_total >= batch_seconds);

  double solve_total = 0.0;
  std::vector<double> solve_samples;
  bool residual_gate_pass = true;
  constexpr int max_solve_repeats = 1000;
  for (int i = 0; i < max_solve_repeats && solve_total < batch_seconds; ++i) {
    Trial trial;
    auto start = Clock::now();
    Solver solver;
    solver.compute(*op);
    solver.preconditioner().addmatrix(preconditioner_matrix);
    solver.setTolerance(tolerance);
    // For solver_lifecycle rows, setup excludes model, operator, RHS, and
    // preconditioner-matrix construction; it times Eigen setup/installation.
    auto after_setup = Clock::now();
    const Eigen::VectorXcd guess = solver.preconditioner().solve(force);
    auto after_guess = Clock::now();
    const Eigen::VectorXcd solution = solver.solveWithGuess(force, guess);
    auto after_solve = Clock::now();
    trial.setup_seconds = elapsed_seconds(start, after_setup);
    trial.guess_seconds = elapsed_seconds(after_setup, after_guess);
    trial.iteration_seconds = elapsed_seconds(after_guess, after_solve);
    trial.total_seconds = trial.setup_seconds + trial.guess_seconds +
                          trial.iteration_seconds;
    trial.solver_info = static_cast<int>(solver.info());
    trial.iterations = solver.iterations();
    trial.solver_error = solver.error();
    trial.absolute_residual = (force - (*op * solution)).norm();
    trial.relative_residual = trial.absolute_residual / force.norm();
    trial.checksum = solution.norm();
    if (solver.info() != Eigen::Success || !solution.allFinite() ||
        !std::isfinite(trial.solver_error) ||
        !std::isfinite(trial.relative_residual) ||
        !std::isfinite(trial.checksum)) {
      std::cerr << "solve failed: scenario=" << scenario
                << " solver_info=" << trial.solver_info
                << " iterations=" << trial.iterations
                << " reported_error=" << trial.solver_error
                << " independent_relative_residual=" << trial.relative_residual
                << " tolerance=" << tolerance << '\n';
      return 2;
    }
    emit_solve(scenario, "solver_lifecycle", i, trial, tolerance);
    if (trial.relative_residual > tolerance) {
      residual_gate_pass = false;
      std::cerr << "residual gate failed: scenario=" << scenario
                << " independent_relative_residual=" << trial.relative_residual
                << " requested_tolerance=" << tolerance << '\n';
    }
    solve_total += trial.total_seconds;
    solve_samples.push_back(trial.total_seconds);
  }
  emit_summary(scenario, "solver_lifecycle_batch", solve_samples,
               solve_samples.empty() ? 0.0 : solve_total / solve_samples.size(),
               batch_seconds, solve_total >= batch_seconds);

  double end_to_end_total = 0.0;
  std::vector<double> end_to_end_samples;
  constexpr int max_end_to_end_repeats = 100;
  for (int i = 0; i < max_end_to_end_repeats &&
                  end_to_end_total < batch_seconds; ++i) {
    Trial trial;
    const auto start = Clock::now();
    auto *const current_stdout = std::cout.rdbuf(library_diagnostics.rdbuf());
    std::unique_ptr<Model> sample_model = build_model();
    std::cout.rdbuf(current_stdout);
    Operator sample_operator(*sample_model);
    const Eigen::VectorXcd sample_force = Gravity_Tools::FindForce(*sample_model);
    Eigen::SparseMatrix<Complex> sample_preconditioner_matrix =
        -sample_model->SpectralElementInformation().fullmatrix<Complex>(0, lmax);
    sample_preconditioner_matrix.makeCompressed();
    Solver sample_solver;
    sample_solver.compute(sample_operator);
    sample_solver.preconditioner().addmatrix(sample_preconditioner_matrix);
    sample_solver.setTolerance(tolerance);
    // For end_to_end_solve rows, setup_seconds starts before model creation and
    // includes model, operator, RHS, matrix, and Eigen solver setup.
    auto after_setup = Clock::now();
    const Eigen::VectorXcd sample_guess =
        sample_solver.preconditioner().solve(sample_force);
    auto after_guess = Clock::now();
    const Eigen::VectorXcd sample_solution =
        sample_solver.solveWithGuess(sample_force, sample_guess);
    auto after_solve = Clock::now();
    trial.total_seconds = elapsed_seconds(start, after_solve);
    trial.setup_seconds = elapsed_seconds(start, after_setup);
    trial.guess_seconds = elapsed_seconds(after_setup, after_guess);
    trial.iteration_seconds = elapsed_seconds(after_guess, after_solve);
    trial.solver_info = static_cast<int>(sample_solver.info());
    trial.iterations = sample_solver.iterations();
    trial.solver_error = sample_solver.error();
    trial.absolute_residual = (sample_force - sample_operator * sample_solution).norm();
    trial.relative_residual = trial.absolute_residual / sample_force.norm();
    trial.checksum = sample_solution.norm();
    if (sample_solver.info() != Eigen::Success || !sample_solution.allFinite() ||
        !std::isfinite(trial.solver_error) ||
        !std::isfinite(trial.relative_residual) || !std::isfinite(trial.checksum))
      throw std::runtime_error("end-to-end solver returned nonfinite or failed output");
    emit_solve(scenario, "end_to_end_solve", i, trial, tolerance);
    if (trial.relative_residual > tolerance) {
      residual_gate_pass = false;
      std::cerr << "residual gate failed: end_to_end scenario=" << scenario
                << " independent_relative_residual=" << trial.relative_residual
                << " requested_tolerance=" << tolerance << '\n';
    }
    end_to_end_total += trial.total_seconds;
    end_to_end_samples.push_back(trial.total_seconds);
  }
  emit_summary(scenario, "end_to_end_solve_batch", end_to_end_samples,
               end_to_end_samples.empty()
                   ? 0.0
                   : end_to_end_total / end_to_end_samples.size(),
               batch_seconds, end_to_end_total >= batch_seconds);
  std::cerr << library_diagnostics.str();
  return residual_gate_pass ? 0 : 3;
}
