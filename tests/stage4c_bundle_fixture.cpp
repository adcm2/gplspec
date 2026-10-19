#include "../gplspec/All"

#include <complex>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

int main(int argc, char **argv) {
   if (argc != 13) throw std::invalid_argument("expected bundle output directories and model input");
   using namespace GeneralEarthModels;
   const auto model = Density3D::SphericalHomogeneousPlanet(
       1.0, 1.0, 6.371e6, 3600.0, 5.972e24, 0.5, 1.2);
   const auto &mesh = model.Node_Information();
   const int lmax = model.GSH_Grid().MaxDegree();
   const std::size_t stride = static_cast<std::size_t>((lmax + 1) * (lmax + 1));
   std::vector<std::complex<double>> solution(
       static_cast<std::size_t>(mesh.NumberOfGlobalRadialNodes()) * stride);
   for (int node = 0; node < mesh.NumberOfGlobalRadialNodes(); ++node)
      solution[static_cast<std::size_t>(node) * stride] = std::sqrt(4.0 * std::numbers::pi);

   GPLSpec::FieldBundle::RunMetadata metadata;
   metadata.workflow = "stage4c_roundtrip";
   metadata.dataset = "homogeneous_potential";
   metadata.solver_method = "CG";
   metadata.solver_tolerance = 1.0e-12;
   metadata.inputs.push_back({"input.txt", std::string(64, 'a')});
   metadata.parameters = {{"radius_norm", "1"}, {"ball_radius_norm", "1.2"}};
   std::vector<double> scalar_profile(
       static_cast<std::size_t>(mesh.NumberOfGlobalRadialNodes()), 2.0);
   std::vector<std::vector<std::complex<double>>> boundary_harmonics(
       static_cast<std::size_t>(mesh.NumberOfElements() + 1),
       std::vector<std::complex<double>>(stride));
   for (auto &row : boundary_harmonics) row[0] = {1.0, 0.0};
   const GPLSpec::FieldBundle::Writer writer(argv[1], model, std::move(metadata));
   writer.Write([&](auto &&emit) {
      GPLSpec::FieldEvaluation::ForEachRadialHarmonicSample(
          model, solution, GPLSpec::FieldEvaluation::Quantity::Potential,
          std::forward<decltype(emit)>(emit));
      GPLSpec::FieldEvaluation::ForEachScalarRadialProfile(
          model, scalar_profile, GPLSpec::FieldEvaluation::Quantity::Potential,
          std::forward<decltype(emit)>(emit));
      GPLSpec::FieldEvaluation::ForEachElementBoundaryHarmonics(
          model, boundary_harmonics,
          GPLSpec::FieldEvaluation::Quantity::SensitivityKernel,
          std::forward<decltype(emit)>(emit));
   });

   const GPLSpec::FieldEvaluation::Rotation rotation{0.2, 0.7, -0.1};
   GPLSpec::FieldBundle::RunMetadata rotated_metadata;
   rotated_metadata.workflow = "stage4c_roundtrip";
   rotated_metadata.dataset = "rotated_full_grid";
   rotated_metadata.solver_method = "CG";
   rotated_metadata.solver_tolerance = 1.0e-12;
   rotated_metadata.inputs.push_back({"input.txt", std::string(64, 'a')});
   rotated_metadata.rotation = rotation;
   GPLSpec::FieldBundle::Writer(argv[2], model, std::move(rotated_metadata)).Write(
       [&](auto &&emit) {
          GPLSpec::FieldEvaluation::ForEachAngularGridSample(
              model, solution, GPLSpec::FieldEvaluation::Quantity::Potential,
              emit, &rotation);
          GPLSpec::FieldEvaluation::ForEachAngularGridSample(
              model, std::span<const std::complex<double>>{},
              GPLSpec::FieldEvaluation::Quantity::ReferentialDensity,
              emit, &rotation);
          GPLSpec::FieldEvaluation::ForEachAngularGridSample(
              model, std::span<const std::complex<double>>{},
              GPLSpec::FieldEvaluation::Quantity::PhysicalDensity,
              emit, &rotation);
       });

   GPLSpec::FieldBundle::RunMetadata surface_metadata;
   surface_metadata.workflow = "stage4c_roundtrip";
   surface_metadata.dataset = "surface_radius";
   const GPLSpec::FieldBundle::Writer surface_writer(
       argv[3], model, std::move(surface_metadata));
   surface_writer.Write([&](auto &&emit) {
      GPLSpec::FieldEvaluation::Sample sample;
      sample.quantity = GPLSpec::FieldEvaluation::Quantity::SurfaceRadius;
      sample.representation = GPLSpec::FieldEvaluation::Representation::AngularGrid;
      sample.radial_node = -1;
      sample.theta_rad = 0.4;
      sample.phi_rad = 1.2;
      sample.reference_radius_m = 6.371e6;
      sample.physical_radius_m = 6.4e6;
      sample.value = {sample.physical_radius_m, 0.0};
      emit(sample);
   });

   const std::vector<double> slice_azimuths{0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
   GPLSpec::FieldBundle::WriteSlice(
       argv[6], model, std::span<const std::complex<double>>{},
       GPLSpec::FieldEvaluation::Quantity::SurfaceRadius, rotation,
       slice_azimuths, [] {
          GPLSpec::FieldBundle::RunMetadata metadata;
          metadata.workflow = "stage4c_roundtrip";
          metadata.dataset = "surface_radius_slice";
          return metadata;
       }());

   const auto phobos_data =
       GPLSpec::ModelInput::ReadPhobosHarmonics(argv[11]);
   constexpr double length_norm = 10000.0;
   constexpr double mass_norm = length_norm * length_norm * length_norm * 1342.0;
   const Density3D full_model(1860.0 / 1342.0, phobos_data, 5, 4,
                              length_norm, 3600.0, mass_norm, 0.1, 1.5, 0);
   const Density3D truncated_model(1860.0 / 1342.0, phobos_data, 0, 5, 4,
                                   length_norm, 3600.0, mass_norm, 0.1, 1.5, 0);
   const auto &full_mesh = full_model.Node_Information();
   const auto full_stride = static_cast<std::size_t>(5 * 5);
   std::vector<std::complex<double>> full_solution(
       static_cast<std::size_t>(full_mesh.NumberOfGlobalRadialNodes()) *
       full_stride);
   GPLSpec::FieldEvaluation::ForEachAngularGridSample(
       full_model, full_solution,
       GPLSpec::FieldEvaluation::Quantity::Potential,
       [&](const GPLSpec::FieldEvaluation::Sample &sample) {
          const auto &grid = full_model.GSH_Grid();
          const auto theta = grid.CoLatitudes();
          const auto phi = grid.Longitudes();
          const auto nearest = [](const auto &values, double target) {
             return static_cast<std::size_t>(std::min_element(
                 values.begin(), values.end(), [target](double a, double b) {
                    return std::abs(a - target) < std::abs(b - target);
                 }) - values.begin());
          };
          const std::size_t spatial = nearest(theta, sample.theta_rad) * phi.size() +
                                     nearest(phi, sample.phi_rad);
          const double expected_radius = sample.reference_radius_m +
              full_model.Mapping_Point(static_cast<int>(sample.element),
                                       sample.radial_node,
                                       static_cast<int>(spatial)) *
                  full_model.LengthNorm();
          const double tolerance = 1.0e-12 * std::max(1.0, std::abs(expected_radius));
          if (std::abs(sample.physical_radius_m - expected_radius) > tolerance ||
              std::abs(sample.x_m - expected_radius * std::sin(sample.theta_rad) *
                                         std::cos(sample.phi_rad)) > tolerance ||
              std::abs(sample.y_m - expected_radius * std::sin(sample.theta_rad) *
                                         std::sin(sample.phi_rad)) > tolerance ||
              std::abs(sample.z_m - expected_radius * std::cos(sample.theta_rad)) > tolerance)
             throw std::runtime_error("full-model mapped coordinates differ from direct geometry");
       });
   GPLSpec::FieldBundle::RunMetadata geometry_metadata;
   geometry_metadata.workflow = "stage4c_geometry";
   geometry_metadata.dataset = "full_model_geometry";
   geometry_metadata.rotation = rotation;
   GPLSpec::FieldBundle::WriteAngular(
       argv[7], full_model, full_solution,
       GPLSpec::FieldEvaluation::Quantity::Potential, geometry_metadata,
       &rotation);
   geometry_metadata.dataset = "truncated_model_geometry";
   GPLSpec::FieldBundle::WriteAngular(
       argv[8], truncated_model, full_solution,
       GPLSpec::FieldEvaluation::Quantity::Potential, geometry_metadata,
       &rotation);
   geometry_metadata.dataset = "full_model_geometry_slice";
   GPLSpec::FieldBundle::WriteSlice(
       argv[9], full_model, full_solution,
       GPLSpec::FieldEvaluation::Quantity::Potential, rotation,
       slice_azimuths, geometry_metadata);
   geometry_metadata.dataset = "truncated_model_geometry_slice";
   GPLSpec::FieldBundle::WriteSlice(
       argv[10], truncated_model, full_solution,
       GPLSpec::FieldEvaluation::Quantity::Potential, rotation,
       slice_azimuths, geometry_metadata);
   GPLSpec::FieldBundle::WriteEquatorialSlice(
       argv[12], model, std::span<const std::complex<double>>{},
       GPLSpec::FieldEvaluation::Quantity::SurfaceRadius, slice_azimuths,
       [] {
          GPLSpec::FieldBundle::RunMetadata metadata;
          metadata.workflow = "stage4c_roundtrip";
          metadata.dataset = "surface_radius_equatorial_slice";
          return metadata;
       }());

   const GPLSpec::FieldEvaluation::Rotation conflicting_rotation{0.4, 0.5, 0.6};
   GPLSpec::FieldBundle::RunMetadata angular_metadata;
   angular_metadata.workflow = "stage4c_roundtrip";
   angular_metadata.dataset = "write_angular_rotation";
   angular_metadata.rotation = conflicting_rotation;
   GPLSpec::FieldBundle::WriteAngular(
       argv[4], model, solution,
       GPLSpec::FieldEvaluation::Quantity::Potential, angular_metadata,
       &rotation);
   angular_metadata.dataset = "write_angular_unrotated";
   GPLSpec::FieldBundle::WriteAngular(
       argv[5], model, solution,
       GPLSpec::FieldEvaluation::Quantity::Potential, angular_metadata);
}
