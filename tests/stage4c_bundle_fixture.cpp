#include "../gplspec/All"

#include <complex>
#include <filesystem>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

int main(int argc, char **argv) {
   if (argc != 3) throw std::invalid_argument("expected regular and rotated bundle directories");
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
}
