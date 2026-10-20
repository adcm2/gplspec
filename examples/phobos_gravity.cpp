#include <GaussQuad/All>
#include <gplspec/All>
#include <gplspec/Test>
#include <gplspec/Timer>
#include <PlanetaryModel/All>
#include <TomographyModels/All>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <ranges>
#include <sstream>
#include <cmath>
#include <numbers>

int
main() {
   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;
   using namespace SimpleModels;
   using FLOAT = double;
   std::string pathtofile = "modeldata/TABLEA1.DAT";
   Timer timer1;
   timer1.start();
   double physdensity = 1860.0;
   double lengthnorm = 10000.0;
   double usedensity = 1342.0;
   double massnorm = std::pow(lengthnorm, 3.0) * usedensity;
   double timenorm = 3600.0;
   double scaledensity = physdensity / massnorm * std::pow(lengthnorm, 3.0);
   int lMax = 100;
   int npoly = 5;
   const auto phobos_data =
       GPLSpec::ModelInput::ReadPhobosHarmonics(pathtofile);
   Density3D phobos(scaledensity, phobos_data, npoly, lMax, lengthnorm, timenorm,
                    massnorm, 0.1, 1.5);
   timer1.stop("Time to construct");

   // std::cout << phobos.GSH_Grid().NumberOfCoLatitudes() << " " <<
   // phobos.GSH_Grid().NumberOfLongitudes() << "\n"; for (auto idx:
   // phobos.GSH_Grid().CoLatitudes()){
   //    std::cout << idx << "\n";
   // }
   // testing model
   // std::cout << "Referential density: " << phobos.Density_Point(0, 0, 0)
   //           << "\n";
   // std::cout << "Physical density: "
   //           << phobos.Density_Point(0, 0, 0) / phobos.Jacobian_Point(0, 0,
   //           0)
   //           << "\n";
   std::cout << "Volume: " << phobos.Volume() / std::pow(10.0, 9.0) << "km^3"
             << "\n";
   std::cout << "Density: " << phobos.Mass() / (phobos.Volume() * 1000.0)
             << " g/cm^3\n";

   // get gravitational field
   timer1.start();
   auto potential_result = SolvePotential(phobos, std::pow(10.0, -12.0));
   timer1.stop("Time for gravity");
   timer1.start();

   // get gravitational field
   timer1.start();
   std::vector<int> sensitivity_l{2}, sensitivity_m{0};
   std::vector<double> sensitivity_weight{1.0};
   auto sensitivity_result = SolveSensitivityKernel(phobos, 2, 0);
   timer1.stop("Time for sensitivity kernel");

   /////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////
   // test rotation
   //  find Euler angles
   double theta1 = std::numbers::pi_v<double> / 2.0;
   double theta2 = std::numbers::pi_v<double> / 2.0;
   theta2 *= 0.0;
   double phi1 = std::numbers::pi_v<double> / 2.0;
   phi1 = 0.0;
   double phi2 = std::numbers::pi_v<double> / 2.0;
   phi2 = 0.0;

   // vectors containing theta,phi information
   std::vector<double> vec_ang1{theta1, phi1}, vec_ang2{theta2, phi2};
   std::vector<double> vec_ang12{theta1, 0.0}, vec_ang22{theta1, theta1};

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
   std::string pathtofolder2 = "./work/Phobos/Spherical";
   std::string pathtofolder3 = "./work/Phobos/Sensitivity";

   auto field_bundle = GPLSpec::FieldBundle::WithSolve(
       potential_result, "phobos_gravity", "potential_density_power", "CG",
       potential_result.final_solve().requested_tolerance);
   field_bundle.inputs.push_back(
       GPLSpec::FieldBundle::InputPath(pathtofile));
   field_bundle.parameters = {{"density_kg_m-3", GPLSpec::FieldBundle::ParameterNumber(physdensity)},
       {"length_norm_m", GPLSpec::FieldBundle::ParameterNumber(lengthnorm)},
       {"mass_norm_kg", GPLSpec::FieldBundle::ParameterNumber(massnorm)},
       {"time_norm_s", GPLSpec::FieldBundle::ParameterNumber(timenorm)},
       {"polynomial_order", std::to_string(npoly)},
       {"degree_max", std::to_string(lMax)},
       {"max_step_norm", "0.1"}, {"exterior_radius_norm", "1.5"}};
   const auto field_bundle_path =
       std::filesystem::path(pathtofolder2) / "field_bundle";
   GPLSpec::FieldBundle::Writer(field_bundle_path, phobos, field_bundle).Write(
       [&](auto &&emit) {
          GPLSpec::FieldEvaluation::ForEachAngularGridSample(
              phobos, potential_result.solution,
              GPLSpec::FieldEvaluation::Quantity::Potential, emit);
          GPLSpec::FieldEvaluation::ForEachAngularGridSample(
              phobos, std::span<const std::complex<double>>{},
              GPLSpec::FieldEvaluation::Quantity::ReferentialDensity, emit);
          GPLSpec::FieldEvaluation::ForEachAngularGridSample(
              phobos, std::span<const std::complex<double>>{},
              GPLSpec::FieldEvaluation::Quantity::PhysicalDensity, emit);
          GPLSpec::FieldEvaluation::ForEachRadialHarmonicSample(
              phobos, potential_result.solution,
              GPLSpec::FieldEvaluation::Quantity::DegreePower, emit);
       });
   auto sensitivity_bundle = GPLSpec::FieldBundle::WithSolve(
       sensitivity_result, "phobos_gravity", "sensitivity_selector_2_0", "CG",
       sensitivity_result.final_solve().requested_tolerance);
   sensitivity_bundle.inputs = field_bundle.inputs;
   sensitivity_bundle.parameters = field_bundle.parameters;
   sensitivity_bundle.selector_l = sensitivity_l;
   sensitivity_bundle.selector_m = sensitivity_m;
   sensitivity_bundle.selector_weight = sensitivity_weight;
   GPLSpec::FieldBundle::Writer(
       std::filesystem::path(pathtofolder3) / "field_bundle", phobos,
       sensitivity_bundle).Write([&](auto &&emit) {
      GPLSpec::FieldEvaluation::ForEachRadialHarmonicSample(
          phobos, sensitivity_result.solution,
          GPLSpec::FieldEvaluation::Quantity::SensitivityKernel, emit);
      GPLSpec::FieldEvaluation::ForEachAngularGridSample(
          phobos, sensitivity_result.solution,
          GPLSpec::FieldEvaluation::Quantity::SensitivityKernel, emit);
   });
   const auto rotation_angles = GPLSpec::detail::RotationOutputAnglesFor(
       vec_ang12, vec_ang22);
   const GPLSpec::FieldEvaluation::Rotation rotation{
       rotation_angles.alpha, rotation_angles.beta, rotation_angles.gamma};
   std::vector<double> slice_azimuths;
   for (double longitude : phobos.GSH_Grid().Longitudes())
      slice_azimuths.push_back(longitude);
   auto physical_rotated_grid = GPLSpec::FieldBundle::WithSolve(
       potential_result, "phobos_gravity", "physical_rotated_angular_grid", "CG",
       potential_result.final_solve().requested_tolerance);
   GPLSpec::FieldBundle::CopyInputProvenance(field_bundle,
                                              physical_rotated_grid);
   physical_rotated_grid.rotation = rotation;
   GPLSpec::FieldBundle::WriteAngular(
       std::filesystem::path(pathtofolder2) / "physical_rotated_grid_bundle",
       phobos, potential_result.solution,
       GPLSpec::FieldEvaluation::Quantity::Potential,
       std::move(physical_rotated_grid), &rotation);
   auto equatorial_metadata = GPLSpec::FieldBundle::WithSolve(
       potential_result, "phobos_gravity", "equatorial_slice", "CG",
       potential_result.final_solve().requested_tolerance);
   GPLSpec::FieldBundle::CopyInputProvenance(field_bundle,
                                              equatorial_metadata);
   GPLSpec::FieldBundle::WriteEquatorialSlice(
       std::filesystem::path(pathtofolder2) / "equatorial_slice_bundle", phobos,
       potential_result.solution, GPLSpec::FieldEvaluation::Quantity::Potential,
       slice_azimuths,
       std::move(equatorial_metadata));
   auto physical_slice_metadata = GPLSpec::FieldBundle::WithSolve(
       potential_result, "phobos_gravity", "physical_rotated_slice", "CG",
       potential_result.final_solve().requested_tolerance);
   GPLSpec::FieldBundle::CopyInputProvenance(field_bundle,
                                              physical_slice_metadata);
   GPLSpec::FieldBundle::WriteSlice(
       std::filesystem::path(pathtofolder2) / "physical_rotated_slice_bundle", phobos,
       potential_result.solution, GPLSpec::FieldEvaluation::Quantity::Potential,
       rotation, slice_azimuths, std::move(physical_slice_metadata));
   const auto referential_rotation_angles =
       GPLSpec::detail::RotationOutputAnglesFor(vec_ang1, vec_ang2);
   const GPLSpec::FieldEvaluation::Rotation referential_rotation{
       referential_rotation_angles.alpha, referential_rotation_angles.beta,
       referential_rotation_angles.gamma};
   auto referential_rotated_grid = GPLSpec::FieldBundle::WithSolve(
       potential_result, "phobos_gravity", "referential_rotated_angular_grid", "CG",
       potential_result.final_solve().requested_tolerance);
   GPLSpec::FieldBundle::CopyInputProvenance(field_bundle,
                                              referential_rotated_grid);
   referential_rotated_grid.rotation = referential_rotation;
   GPLSpec::FieldBundle::WriteAngular(
       std::filesystem::path(pathtofolder2) / "referential_rotated_grid_bundle",
       phobos, potential_result.solution,
       GPLSpec::FieldEvaluation::Quantity::Potential,
       std::move(referential_rotated_grid), &referential_rotation);
   auto referential_rotated_metadata = GPLSpec::FieldBundle::WithSolve(
       potential_result, "phobos_gravity", "referential_rotated_slice", "CG",
       potential_result.final_solve().requested_tolerance);
   GPLSpec::FieldBundle::CopyInputProvenance(field_bundle,
                                              referential_rotated_metadata);
   GPLSpec::FieldBundle::WriteSlice(
       std::filesystem::path(pathtofolder2) / "referential_rotated_slice_bundle",
       phobos, potential_result.solution,
       GPLSpec::FieldEvaluation::Quantity::Potential, referential_rotation,
       slice_azimuths, std::move(referential_rotated_metadata));
   GPLSpec::FieldBundle::WriteSlice(
       std::filesystem::path(pathtofolder3) / "sensitivity_slice_bundle", phobos,
       sensitivity_result.solution,
       GPLSpec::FieldEvaluation::Quantity::SensitivityKernel, rotation,
       slice_azimuths, sensitivity_bundle);
   auto density_rotation = GPLSpec::FieldBundle::RunMetadata{};
   density_rotation.workflow = "phobos_gravity";
   density_rotation.dataset = "rotated_density_traces";
   GPLSpec::FieldBundle::CopyInputProvenance(field_bundle,
                                              density_rotation);
   density_rotation.parameters = field_bundle.parameters;
   density_rotation.rotation = rotation;
   GPLSpec::FieldBundle::Writer(
       std::filesystem::path(pathtofolder2) / "rotated_density_slice_bundle",
       phobos, density_rotation).Write([&](auto &&emit) {
      GPLSpec::FieldEvaluation::ForEachRotatedSliceSample(
          phobos, std::span<const std::complex<double>>{},
          GPLSpec::FieldEvaluation::Quantity::ReferentialDensity, rotation,
          std::span<const double>(slice_azimuths), emit);
      GPLSpec::FieldEvaluation::ForEachRotatedSliceSample(
          phobos, std::span<const std::complex<double>>{},
          GPLSpec::FieldEvaluation::Quantity::PhysicalDensity, rotation,
          std::span<const double>(slice_azimuths), emit);
   });
   timer1.stop("Physical rotated");

   return 0;
}
