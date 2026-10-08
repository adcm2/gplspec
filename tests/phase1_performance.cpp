#include <gplspec/All>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>

namespace {
using Clock = std::chrono::steady_clock;
using Complex = std::complex<double>;
using Model = GeneralEarthModels::Density3D;
using Operator = MatrixReplacement3D<Complex>;
using Preconditioner = Eigen::SphericalGeometryPreconditioner<Complex>;
using Solver = Eigen::ConjugateGradient<Operator, Eigen::Lower | Eigen::Upper,
                                        Preconditioner>;
constexpr double tolerance = 1e-6;
constexpr int warm_repeats = 5;

double elapsed(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

void emit(const std::string &process, const std::string &scenario,
          const std::string &metric, int repeat, double value,
          const std::string &unit) {
  std::cout << process << ',' << scenario << ',' << metric << ',' << repeat
            << ',' << value << ',' << unit << '\n';
}

class CoutCapture {
 public:
  CoutCapture() : previous_(std::cout.rdbuf(sink_.rdbuf())) {}
  ~CoutCapture() { std::cout.rdbuf(previous_); }
  std::string text() const { return sink_.str(); }

 private:
  std::ostringstream sink_;
  std::streambuf *previous_;
};

int iteration_count(const std::string &diagnostics) {
  const std::string label = "Number of iterations:";
  const auto pos = diagnostics.find(label);
  if (pos == std::string::npos)
    throw std::runtime_error("public solver did not report its iteration count");
  return std::stoi(diagnostics.substr(pos + label.size()));
}

double maximum_norm(const auto &values) {
  double result = 0.0;
  for (const auto &element : values)
    for (const auto &node : element)
      for (const Complex value : node) {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
          throw std::runtime_error("benchmark solve returned nonfinite output");
        result = std::max(result, std::abs(value));
      }
  return result;
}

struct SmoothMapping {
  auto RadialMapping(int) const {
    return [](double r, double theta, double phi) {
      return 0.1 * r * (1.0 - r) * std::sin(theta) * std::cos(phi);
    };
  }
};

struct ZeroMapping {
  auto RadialMapping(int) const {
    return [](double, double, double) { return 0.0; };
  }
};

void benchmark_mapped(const std::string &process) {
  const auto base_start = Clock::now();
  auto base = GeneralEarthModels::SimpleModels::spherical_model::HomogeneousSphere(
      1.0, 1.0, 6371000.0, 3600.0, 5.972e24);
  emit(process, "mapped", "base_model_setup", 0, elapsed(base_start), "s");

  SmoothMapping mapping;
  std::optional<Model> model_storage;
  double model_seconds = 0.0;
  {
    CoutCapture capture;
    const auto model_start = Clock::now();
    model_storage.emplace(base, PlanetaryModel::TomographyZeroModel(), mapping,
                          2, 2, 0.5, 1.2);
    model_seconds = elapsed(model_start);
  }
  emit(process, "mapped", "model_discretization_setup", 0, model_seconds, "s");
  Model &model = *model_storage;
  Operator op(model);
  emit(process, "mapped", "operator_rows", 0, op.rows(), "count");

  Eigen::VectorXcd x(op.cols());
  for (Eigen::Index i = 0; i < x.size(); ++i)
    x(i) = Complex(0.125 * (i + 1), -0.0625 * (i + 2));
  Eigen::VectorXcd action;
  auto start = Clock::now();
  action = op * x;
  emit(process, "mapped", "matvec_first", 0, elapsed(start), "s");
  if (!action.allFinite() || action.norm() == 0.0)
    throw std::runtime_error("mapped complex operator action is invalid");
  emit(process, "mapped", "matvec_norm", 0, action.norm(), "2-norm");
  for (int repeat = 1; repeat <= warm_repeats; ++repeat) {
    start = Clock::now();
    action = op * x;
    emit(process, "mapped", "matvec_repeat", repeat, elapsed(start), "s");
  }

  for (int repeat = 0; repeat <= warm_repeats; ++repeat) {
    std::string diagnostics;
    start = Clock::now();
    decltype(Gravity_Tools::FindGravitationalPotential(model, tolerance)) solution;
    {
      CoutCapture capture;
      solution = Gravity_Tools::FindGravitationalPotential(model, tolerance);
      diagnostics = capture.text();
    }
    const double seconds = elapsed(start);
    emit(process, "mapped",
         repeat == 0 ? "public_full_solve_first" : "public_full_solve_repeat",
         repeat, seconds, "s");
    emit(process, "mapped", "public_solver_iterations", repeat,
         iteration_count(diagnostics), "count");
    emit(process, "mapped", "public_solution_max_norm", repeat,
         maximum_norm(solution), "absolute");
  }

  Eigen::SparseMatrix<Complex> preconditioner_matrix =
      -model.SpectralElementInformation().fullmatrix<Complex>(
          0, model.GSH_Grid().MaxDegree());
  preconditioner_matrix.makeCompressed();
  const Eigen::VectorXcd force = Gravity_Tools::FindForce(model);
  Solver solver;
  solver.compute(op);
  solver.preconditioner().addmatrix(preconditioner_matrix);
  solver.setTolerance(tolerance);
  const Eigen::VectorXcd guess = solver.preconditioner().solve(force);
  const Eigen::VectorXcd solution = solver.solveWithGuess(force, guess);
  const double absolute_residual = (force - op * solution).norm();
  const double relative_residual = absolute_residual / force.norm();
  if (solver.info() != Eigen::Success || solver.iterations() < 2 ||
      !solution.allFinite() || !std::isfinite(solver.error()) ||
      !std::isfinite(absolute_residual) ||
      !std::isfinite(relative_residual) || relative_residual > tolerance)
    throw std::runtime_error("mapped benchmark solver diagnostics failed");
  emit(process, "mapped", "diagnostic_iterations", 0,
       solver.iterations(), "count");
  emit(process, "mapped", "solver_reported_error", 0,
       solver.error(), "relative");
  emit(process, "mapped", "independent_absolute_residual", 0,
       absolute_residual, "2-norm");
  emit(process, "mapped", "independent_relative_residual", 0,
       relative_residual, "relative");
}

void benchmark_file_tomography(const std::string &process) {
  const auto earth_start = Clock::now();
  TestTools::EarthModel<> earth(
      GPLSpec::ModelInput::ReadPremDeck("modeldata/PREM500"));
  emit(process, "file_tomography", "earth_input_read", 0,
       elapsed(earth_start), "s");

  const auto tomography_start = Clock::now();
  Tomography tomography("modeldata/S40RTS_dvs.nc");
  emit(process, "file_tomography", "netcdf_input_read", 0,
       elapsed(tomography_start), "s");

  ZeroMapping mapping;
  std::optional<Model> model_storage;
  double model_seconds = 0.0;
  {
    CoutCapture capture;
    const auto model_start = Clock::now();
    model_storage.emplace(earth, tomography, mapping, 2, 2, 0.5, 1.2);
    model_seconds = elapsed(model_start);
  }
  emit(process, "file_tomography", "model_discretization_setup", 0,
       model_seconds, "s");
  Model &model = *model_storage;
  Operator op(model);
  emit(process, "file_tomography", "operator_rows", 0, op.rows(), "count");
  Eigen::VectorXcd x(op.cols());
  for (Eigen::Index i = 0; i < x.size(); ++i)
    x(i) = Complex(0.125 * (i + 1), -0.0625 * (i + 2));
  const auto action_start = Clock::now();
  const Eigen::VectorXcd action = op * x;
  emit(process, "file_tomography", "matvec_first", 0,
       elapsed(action_start), "s");
  if (!action.allFinite() || action.norm() == 0.0)
    throw std::runtime_error("file-backed tomography operator action is invalid");
  emit(process, "file_tomography", "matvec_norm", 0,
       action.norm(), "2-norm");
}
}  // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: phase1_performance PROCESS_ID\n";
    return 2;
  }
  try {
    std::cout << std::setprecision(17)
              << "process,scenario,metric,repeat,value,unit\n";
    const std::string process = argv[1];
    benchmark_mapped(process);
    benchmark_file_tomography(process);
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0)
      throw std::runtime_error("getrusage failed");
    emit(process, "process", "peak_rss", 0, usage.ru_maxrss, "KiB");
  } catch (const std::exception &error) {
    std::cerr << "phase1_performance: " << error.what() << '\n';
    return 1;
  }
}
