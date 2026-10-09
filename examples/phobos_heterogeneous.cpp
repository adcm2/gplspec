/**
 * @file phobos_heterogeneous.cpp
 * @brief Example for computing gravitational field and sensitivity kernel for
 * Phobos using heterogeneous density.
 *
 * This example reads a density model for Phobos, computes gravitational
 * potential solutions, sensitivity kernels, and outputs results to files in
 * various formats.
 */

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

/**
 * @brief Main routine for Phobos gravity and sensitivity kernel computation.
 *
 * - Reads density model data.
 * - Constructs Density3D object for Phobos.
 * - Computes gravitational potential and sensitivity kernel.
 * - Outputs results to files.
 */
int
main() {
   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;
   using namespace SimpleModels;
   using FLOAT = double;

   // Path to density model data file
   std::string pathtofile = "modeldata/TABLEA1.DAT";

   Timer timer1;

   // Physical and model parameters
   double physdensity = 1860.0;
   double lengthnorm = 10000.0;
   double usedensity = 1342.0;
   double massnorm = std::pow(lengthnorm, 3.0) * usedensity;
   double timenorm = 3600.0;
   double scaledensity = physdensity / massnorm * std::pow(lengthnorm, 3.0);
   int lMax = 100;
   int npoly = 5;
   bool nhomo = true;

   // Construct Density3D model for Phobos
   timer1.start();
   const auto phobos_data =
       GPLSpec::ModelInput::ReadPhobosHarmonics(pathtofile);
   Density3D phobos(scaledensity, phobos_data, npoly, lMax, lengthnorm, timenorm,
                    massnorm, 0.1, 1.5, nhomo);
   timer1.stop("Time to construct");

   // Compute gravitational potential solutions
   timer1.start();
   auto potential_result = SolvePotential(phobos, std::pow(10.0, -5.0));
   auto stdvec_potsol = ToHarmonicField(phobos, potential_result);
   timer1.stop("Time for gravity");

   timer1.start();
   auto stdvec_potsol2 = ToHarmonicField(
       phobos, SolvePotential(phobos, std::pow(10.0, -3.0)));
   timer1.stop("Time for gravity 2");

   // Compute spherical harmonic sensitivity kernel
   timer1.start();
   std::vector<int> sensitivity_l{2}, sensitivity_m{0};
   std::vector<double> sensitivity_weight{1.0};
   auto sensitivity_result = SolveSensitivityKernel(phobos, 2, 0);
   auto stdvec_senskernel = ToHarmonicField(phobos, sensitivity_result);
   timer1.stop("Time for sensitivity kernel");

   // Compute power spectrum
   auto stdpower = phobos.PowerSTD(stdvec_potsol);

   // Rotation angles for output
   double theta1 = std::numbers::pi_v<double> / 2.0;
   double theta2 = 0.0;
   double phi1 = 0.0;
   double phi2 = 0.0;

   std::vector<double> vec_ang1{theta1, phi1}, vec_ang2{theta2, phi2};
   std::vector<double> vec_ang12{theta1, 0.0}, vec_ang22{theta1, theta1};

   // Output paths
   std::string pathtofolder1 = "./work/Phobos/Cartesian";
   std::string pathtofolder2 = "./work/Phobos/Spherical";
   std::string pathtofolder3 = "./work/Phobos/Sensitivity";
   std::string pathtofile1 =
       pathtofolder2 + "/MatrixSolutionReferentialRotated.out";
   std::string pathtofile2 =
       pathtofolder2 + "/MatrixSolutionPhysicalRotated.out";
   std::string pathtofile3 = pathtofolder2 + "/MatrixSolutionRotated.out";
   std::string pathtofile4 =
       pathtofolder2 + "/PhysicalDensitySolutionRotated.out";
   std::string pathtofile5 =
       pathtofolder2 + "/ReferentialDensitySolutionRotated.out";
   std::string pathpower = "./work/Phobos/Spherical/PowerSpectrum.out";

   // Output results
   phobos.CartesianOutputAtElement(pathtofolder1, stdvec_potsol);
   phobos.PhysicalOutputAtElement(pathtofolder2, stdvec_potsol);
   phobos.ReferentialOutputAtElement(pathtofolder3, stdvec_senskernel, true);
   phobos.ReferentialOutputSlice(pathtofolder3, stdvec_senskernel, true);
   phobos.ReferentialOutputRotated(pathtofile1, vec_ang12, vec_ang22,
                                   stdvec_potsol);
   phobos.PhysicalOutputRotated(pathtofile2, vec_ang12, vec_ang22,
                                stdvec_potsol);
   phobos.ModelDensityOutputRotated(pathtofile4, vec_ang12, vec_ang22, true);
   phobos.ModelDensityOutputRotated(pathtofile5, vec_ang12, vec_ang22);

   // Output power spectrum to file
   std::ofstream file(pathpower);
   for (const auto &row : stdpower) {
      for (const auto &col : row) {
         file << std::setprecision(16) << col[0];
         for (std::size_t k = 1; k < col.size(); ++k) {
            file << ";" << col[k];
         }
         file << std::endl;
      }
   }

   auto field_bundle = GPLSpec::FieldBundle::WithSolve(
       potential_result, "phobos_heterogeneous", "potential_density_power", "CG",
       potential_result.final_solve().requested_tolerance);
   field_bundle.inputs.push_back(
       GPLSpec::FieldBundle::InputPath(pathtofile));
   field_bundle.parameters = {{"density_kg_m-3", std::to_string(physdensity)},
       {"length_norm_m", std::to_string(lengthnorm)},
       {"mass_norm_kg", std::to_string(massnorm)},
       {"time_norm_s", std::to_string(timenorm)},
       {"polynomial_order", std::to_string(npoly)},
       {"degree_max", std::to_string(lMax)},
       {"max_step_norm", "0.1"}, {"exterior_radius_norm", "1.5"},
       {"homogeneous_flag", nhomo ? "true" : "false"}};
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
       sensitivity_result, "phobos_heterogeneous", "sensitivity_selector_2_0",
       "CG", sensitivity_result.final_solve().requested_tolerance);
   sensitivity_bundle.inputs = field_bundle.inputs;
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
   auto rotated_grid = GPLSpec::FieldBundle::WithSolve(
       potential_result, "phobos_heterogeneous", "rotated_potential_angular_grid",
       "CG", potential_result.final_solve().requested_tolerance);
   rotated_grid.inputs = field_bundle.inputs;
   rotated_grid.rotation = rotation;
   GPLSpec::FieldBundle::WriteAngular(
       std::filesystem::path(pathtofolder2) / "rotated_grid_bundle", phobos,
       potential_result.solution, GPLSpec::FieldEvaluation::Quantity::Potential,
       std::move(rotated_grid), &rotation);
   std::vector<double> slice_azimuths;
   for (double longitude : phobos.GSH_Grid().Longitudes())
      slice_azimuths.push_back(longitude);
   auto rotated_potential = GPLSpec::FieldBundle::WithSolve(
       potential_result, "phobos_heterogeneous", "rotated_equatorial_slice", "CG",
       potential_result.final_solve().requested_tolerance);
   rotated_potential.inputs = field_bundle.inputs;
   GPLSpec::FieldBundle::WriteSlice(
       std::filesystem::path(pathtofolder2) / "rotated_slice_bundle", phobos,
       potential_result.solution, GPLSpec::FieldEvaluation::Quantity::Potential,
       rotation, slice_azimuths, rotated_potential);
   auto rotated_density_grid = GPLSpec::FieldBundle::RunMetadata{};
   rotated_density_grid.workflow = "phobos_heterogeneous";
   rotated_density_grid.dataset = "rotated_density_angular_grid";
   rotated_density_grid.inputs = field_bundle.inputs;
   rotated_density_grid.rotation = rotation;
   GPLSpec::FieldBundle::Writer(
       std::filesystem::path(pathtofolder2) / "rotated_density_grid_bundle",
       phobos, rotated_density_grid).Write([&](auto &&emit) {
      GPLSpec::FieldEvaluation::ForEachAngularGridSample(
          phobos, std::span<const std::complex<double>>{},
          GPLSpec::FieldEvaluation::Quantity::ReferentialDensity, emit, &rotation);
      GPLSpec::FieldEvaluation::ForEachAngularGridSample(
          phobos, std::span<const std::complex<double>>{},
          GPLSpec::FieldEvaluation::Quantity::PhysicalDensity, emit, &rotation);
   });
   GPLSpec::FieldBundle::WriteSlice(
       std::filesystem::path(pathtofolder3) / "sensitivity_slice_bundle", phobos,
       sensitivity_result.solution,
       GPLSpec::FieldEvaluation::Quantity::SensitivityKernel, rotation,
       slice_azimuths, sensitivity_bundle);
   auto density_rotation = GPLSpec::FieldBundle::RunMetadata{};
   density_rotation.workflow = "phobos_heterogeneous";
   density_rotation.dataset = "rotated_density_traces";
   density_rotation.inputs = field_bundle.inputs;
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
   file.close();

   return 0;
}
