#include <gplspec/All>

#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>
#include <sstream>
#include <type_traits>

namespace {
template <class T>
void emit(std::ostream &out, const std::string &path, const T &value) {
  if constexpr (std::is_arithmetic_v<T>) {
    out << path << ',' << value << '\n';
  } else if constexpr (requires { value.rows(); value.cols(); value(0, 0); }) {
    for (Eigen::Index i = 0; i < value.rows(); ++i)
      for (Eigen::Index j = 0; j < value.cols(); ++j)
        emit(out, path + '/' + std::to_string(i) + '/' + std::to_string(j),
             value(i, j));
  } else if constexpr (requires { value.real(); value.imag(); }) {
    out << path << ',' << value.real() << ',' << value.imag() << '\n';
  } else if constexpr (requires { value.begin(); value.end(); }) {
    std::size_t i = 0;
    for (const auto &entry : value)
      emit(out, path + '/' + std::to_string(i++), entry);
  } else {
    static_assert(sizeof(T) == 0, "Unsupported stage05 reference value");
  }
}

struct ZeroMapping {
  auto RadialMapping(int) const {
    return [](double, double, double) { return 0.0; };
  }
};
}  // namespace

int main() {
  using Complex = std::complex<double>;
  const double length = 6371000.0;
  const double time = 3600.0;
  const double mass = 5.972e24;
  auto material = GeneralEarthModels::SimpleModels::spherical_model::
      HomogeneousSphere(1.0, 1.0, length, time, mass);
  ZeroMapping zero;
  GeneralEarthModels::Density3D model(
      material, PlanetaryModel::TomographyZeroModel(), zero, 3, 2, 0.5, 1.2);
  const std::string coefficients =
      std::string(GPLSPEC_SOURCE_DIR) +
      "/tests/reference/stage05-file-perturbation-coefficients.dat";
  GeneralEarthModels::MappingPerturbation perturbation(model, coefficients, 0,
                                                        2);

  std::cout << std::setprecision(17);
  emit(std::cout, "file_perturbation/dxi", perturbation.dxi());
  emit(std::cout, "file_perturbation/dxilm", perturbation.dxilm());
  emit(std::cout, "file_perturbation/df", perturbation.df());
  emit(std::cout, "file_perturbation/da", perturbation.da());
  emit(std::cout, "file_perturbation/source", Gravity_Tools::FindForce(model));

  MatrixReplacement3D<Complex> op(model, perturbation);
  Eigen::VectorXcd x(op.cols()), y(op.cols());
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    x(i) = Complex(0.125 * (i + 1), -0.0625 * (i + 2));
    y(i) = Complex(-0.03125 * (i + 3), 0.09375 * (i + 1));
  }
  const Eigen::VectorXcd ax = op * x;
  const Eigen::VectorXcd ay = op * y;
  emit(std::cout, "file_perturbation/operator_x", x);
  emit(std::cout, "file_perturbation/operator_y", ax);
  emit(std::cout, "file_perturbation/operator_ay", ay);
  const Complex xay = x.dot(ay);
  const Complex ax_y = ax.dot(y);
  emit(std::cout, "file_perturbation/xAy", xay);
  emit(std::cout, "file_perturbation/Ax_y", ax_y);
  std::ostringstream solver_diagnostics;
  auto *original_stdout = std::cout.rdbuf(solver_diagnostics.rdbuf());
  auto solution = Gravity_Tools::FindGravitationalPotentialPerturbation(
      model, perturbation, 1e-6);
  std::cout.rdbuf(original_stdout);
  emit(std::cout, "file_perturbation/solution", solution);
}
