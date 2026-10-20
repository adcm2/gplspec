#include <gplspec/All>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <type_traits>
#include <complex>
#include <numbers>

namespace {
template<class T>
void emit(std::ostream& out, const std::string& path, const T& value) {
  if constexpr (std::is_arithmetic_v<T>) {
    out << path << ',' << value << '\n';
  } else if constexpr (requires { value.rows(); value.cols(); value(0,0); }) {
    for (Eigen::Index i = 0; i < value.rows(); ++i)
      for (Eigen::Index j = 0; j < value.cols(); ++j)
        emit(out, path + '/' + std::to_string(i) + '/' + std::to_string(j), value(i,j));
  } else if constexpr (requires { value.real(); value.imag(); }) {
    out << path << ',' << value.real() << ',' << value.imag() << '\n';
  } else if constexpr (requires { value.begin(); value.end(); }) {
    std::size_t i = 0;
    for (const auto& item : value) emit(out, path + '/' + std::to_string(i++), item);
  } else {
    static_assert(sizeof(T) == 0, "Unsupported stage00 reference value");
  }
}

template<class Model>
void write_legacy_referential_field(
    const std::filesystem::path &path, const Model &model,
    const Eigen::VectorXcd &solution) {
  std::ofstream file(path);
  const auto last_element = static_cast<std::size_t>(model.Num_Elements() - 1);
  GPLSpec::FieldEvaluation::ForEachRadialHarmonicSample(
      model, solution, GPLSpec::FieldEvaluation::Quantity::Potential,
      [&](const GPLSpec::FieldEvaluation::Sample &sample) {
        if (sample.radial_node != 0 &&
            !(sample.element == last_element &&
              sample.radial_node == model.Poly_Order())) return;
        if (sample.degree == 0 && sample.order == 0)
          file << std::setprecision(16) << sample.reference_radius_m;
        file << ';' << std::setprecision(16) << sample.value.real()
             << ';' << sample.value.imag();
        if (sample.degree == model.GSH_Grid().MaxDegree() &&
            sample.order == sample.degree)
          file << '\n';
      });
}

template<class Model>
auto capture(std::ostream& out, const std::string& name, Model& model,
             const std::filesystem::path &output_dir) {
  emit(out, name + "/density", model.Density());
  emit(out, name + "/mapping", model.Mapping());
  emit(out, name + "/jacobian", model.Jacobian());
  emit(out, name + "/inverse_f", model.InverseFRef());
  emit(out, name + "/laplace_tensor", model.LaplaceTensorRef());
  emit(out, name + "/source", Gravity_Tools::FindForce(model));
  MatrixReplacement3D<std::complex<double>> op(model);
  Eigen::VectorXcd x(op.cols());
  for (Eigen::Index i = 0; i < x.size(); ++i)
    x(i) = {0.125 * (i + 1), -0.0625 * (i + 2)};
  emit(out, name + "/operator_x", x);
  emit(out, name + "/operator_y", op * x);
  auto result = Gravity_Tools::SolvePotential(model, 1e-6);
  auto field = Gravity_Tools::ToHarmonicField(model, result);
  emit(out, name + "/solution", field);
  auto metadata = GPLSpec::FieldBundle::WithSolve(
      result, "stage00_reference", name, "CG", 1e-6);
  metadata.parameters = {{"length_norm_m", std::to_string(model.LengthNorm())},
      {"mass_norm_kg", std::to_string(model.MassNorm())},
      {"time_norm_s", std::to_string(model.TimeNorm())}};
  GPLSpec::FieldBundle::WriteAngular(
      output_dir / "field_bundles" / name, model, result.solution,
      GPLSpec::FieldEvaluation::Quantity::Potential, std::move(metadata));
  return result.solution;
}

template<class Model, class Map>
void capture_perturbation(std::ostream& out, const std::string& name,
                          Model& model, const Map& mapping,
                          const std::filesystem::path &output_dir) {
  GeneralEarthModels::MappingPerturbation perturbation(model, mapping);
  emit(out, name + "/da", perturbation.delta_laplace_tensor());
  MatrixReplacement3D<std::complex<double>> op(model, perturbation);
  Eigen::VectorXcd x(op.cols());
  for (Eigen::Index i = 0; i < x.size(); ++i)
    x(i) = {0.125 * (i + 1), -0.0625 * (i + 2)};
  emit(out, name + "/operator_x", x);
  emit(out, name + "/operator_y", op * x);
  const auto result = Gravity_Tools::SolveMappingPerturbation(
      model, perturbation, 1e-6);
  emit(out, name + "/solution", Gravity_Tools::ToHarmonicField(model, result));
  auto metadata = GPLSpec::FieldBundle::WithSolve(
      result, "stage00_reference", name, "CG", 1e-6);
  GPLSpec::FieldBundle::WriteAngular(
      output_dir / "field_bundles" / name, model, result.solution,
      GPLSpec::FieldEvaluation::Quantity::Potential, std::move(metadata));
}

struct SmoothMapping {
  auto RadialMapping(int) const {
    return [](double r, double theta, double phi) {
      return 0.01 * r * (1.0 - r) * std::sin(theta) * std::cos(phi);
    };
  }
};
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: stage00_reference OUTPUT_DIR ACTUAL.csv\n";
    return 2;
  }
  std::filesystem::create_directories(argv[1]);
  std::ofstream out(argv[2]);
  if (!out) return 3;
  out << std::setprecision(17);
  using namespace GeneralEarthModels;
  const double length = 6371000.0, time = 3600.0, mass = 5.972e24;
  auto sphere = Density3D::SphericalHomogeneousPlanet(1.0, 1.0, length, time,
                                                       mass, 0.5, 1.2);
  auto sphere_solution = capture(out, "homogeneous", sphere, argv[1]);
  const std::string output_dir = std::string(argv[1]) + "/homogeneous";
  std::filesystem::create_directories(output_dir);
  write_legacy_referential_field(
      std::filesystem::path(output_dir) / "MatrixSolution.out", sphere,
      sphere_solution);

  auto layered = Density3D::SphericalLayeredPlanet({0.5, 1.0}, {2.0, 1.0},
                                                    length, time, mass, 0.5, 1.2);
  capture(out, "layered", layered, argv[1]);

  std::vector<double> depths{0.0, 6371.0}, longitudes{0.0, 90.0, 180.0, 270.0},
      latitudes{-90.0, 0.0, 90.0}, values;
  for (std::size_t d = 0; d < depths.size(); ++d)
    for (double lat : latitudes)
      for (double lon : longitudes)
        values.push_back(0.2 * std::sin(lat * std::numbers::pi / 180.0) *
                         std::cos(lon * std::numbers::pi / 180.0));
  Tomography lateral_density(depths, longitudes, latitudes, values);
  auto base = SimpleModels::spherical_model::HomogeneousSphere(
      1.0, 1.0, length, time, mass);
  struct ZeroMapping { auto RadialMapping(int) const { return [](double, double, double) { return 0.0; }; } } zero;
  Density3D lateral(base, lateral_density, zero, 3, 2, 0.5, 1.2);
  capture(out, "lateral_density", lateral, argv[1]);

  SmoothMapping smooth;
  Density3D mapped(base, PlanetaryModel::TomographyZeroModel(), smooth,
                   3, 2, 0.5, 1.2);
  capture(out, "smooth_mapping", mapped, argv[1]);
  capture_perturbation(out, "mapping_perturbation", mapped, smooth, argv[1]);
  return out.good() ? 0 : 4;
}
