#include <GaussQuad/All>
#include <gplspec/All>
#include <gplspec/Test>
#include <gplspec/Timer>
#include <PlanetaryModel/All>
#include <TomographyModels/All>
#include <GSHTrans/All>
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
   using INT = std::ptrdiff_t;
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
   int lsmall = 0;
   std::cout << "lmax for base calculation: ";
   std::cin >> lsmall;
   int npoly = 5;
   using MRange = GSHTrans::All;
   using NRange = GSHTrans::All;
   using Grid = GSHTrans::GaussLegendreGrid<double, MRange, NRange>;
   timer1.start();
   const auto phobos_data =
       GPLSpec::ModelInput::ReadPhobosHarmonics(pathtofile);
   Density3D phobosfull(scaledensity, phobos_data, npoly, lMax, lengthnorm,
                        timenorm, massnorm, 0.1, 1.5, 0);
   timer1.stop("Time to construct full model");
   timer1.start();
   Density3D phobos(scaledensity, phobos_data, lsmall, npoly, lMax, lengthnorm,
                    timenorm, massnorm, 0.1, 1.5, 0);
   timer1.stop("Time to construct small model");

   /////////////////////////////////////////////////////////////////////////////
   // exact

   // get gravitational field
   timer1.start();
   auto full_result = SolvePotential(phobosfull, std::pow(10.0, -6.0));
   auto stdvec_potsol = ToHarmonicField(phobosfull, full_result);
   timer1.stop("Time for gravity with full model");

   timer1.start();
   auto small_result = SolvePotential(phobos, std::pow(10.0, -6.0));
   auto stdvec_potsol_small = ToHarmonicField(phobos, small_result);
   timer1.stop("Time for gravity with small model");
   /////////////////////////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////////////////////
   // find perturbed solution
   // define map
   timer1.start();
   MappingPerturbation testpertmap(phobos, phobos_data, lsmall + 1, lMax);
   timer1.stop("Time for map to be defined");

   // get solution
   timer1.start();
   auto mapping_result = SolveMappingPerturbation(
       phobos, testpertmap, std::pow(10.0, -6.0), std::pow(10.0, -6.0));
   auto stdvec_potsol_perturb = ToHarmonicField(phobos, mapping_result);
   timer1.stop("Time for perturbation calculation");

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

   /////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////
   // output test
   std::string pathtofolder1 = "./work/Phobos/Cartesian";
   std::string pathtofolder2 = "./work/Phobos/Spherical";
   std::string pathtofolder3 = "./work/Phobos/Rotated";
   std::string pathtofolder4 = "./work/Phobos/Perturbed";
   std::string pathtofile1 = pathtofolder3 + "/MatrixSolution.out";
   phobosfull.CartesianOutputAtElement(pathtofolder1, stdvec_potsol);
   phobosfull.PhysicalOutputAtElement(pathtofolder2, stdvec_potsol);
   phobos.PhysicalOutputRotated(pathtofile1, vec_ang1, vec_ang2, stdvec_potsol);
   phobos.PhysicalOutputAtElement(pathtofolder4, stdvec_potsol_perturb);

   const std::vector<std::pair<std::string, std::string>> run_parameters = {
       {"physical_density_kg_m-3", std::to_string(physdensity)},
       {"scaled_density", std::to_string(scaledensity)},
       {"polynomial_order", std::to_string(npoly)},
       {"degree_max", std::to_string(lMax)},
       {"base_degree", std::to_string(lsmall)},
       {"max_step_norm", "0.1"}, {"exterior_radius_norm", "1.5"},
       {"solve_tolerance", "1e-6"},
       {"first_tolerance", "1e-6"}, {"correction_tolerance", "1e-6"}};

   auto full_bundle = GPLSpec::FieldBundle::WithSolve(
       full_result, "phobos_gravity_perturbation", "full_degree_potential", "CG",
       full_result.final_solve().requested_tolerance);
   full_bundle.inputs.push_back(
       GPLSpec::FieldBundle::InputPath(pathtofile));
   full_bundle.parameters = run_parameters;
   GPLSpec::FieldBundle::WriteAngular(
       std::filesystem::path(pathtofolder2) / "field_bundle", phobosfull,
       full_result.solution, GPLSpec::FieldEvaluation::Quantity::Potential,
       std::move(full_bundle));
   auto base_bundle = GPLSpec::FieldBundle::WithSolve(
       small_result, "phobos_gravity_perturbation", "base_small_model", "CG",
       small_result.final_solve().requested_tolerance);
   base_bundle.inputs.push_back(
       GPLSpec::FieldBundle::InputPath(pathtofile));
   base_bundle.parameters = run_parameters;
   GPLSpec::FieldBundle::WriteAngular(
       std::filesystem::path(pathtofolder2) / "base_field_bundle", phobos,
       small_result.solution, GPLSpec::FieldEvaluation::Quantity::Potential,
       std::move(base_bundle));
   auto perturb_bundle = GPLSpec::FieldBundle::WithSolve(
       mapping_result, "phobos_gravity_perturbation", "mapping_perturbation", "CG",
       mapping_result.final_solve().requested_tolerance);
   perturb_bundle.inputs.push_back(
       GPLSpec::FieldBundle::InputPath(pathtofile));
   perturb_bundle.parameters = run_parameters;
   GPLSpec::FieldBundle::WriteAngular(
       std::filesystem::path(pathtofolder4) / "field_bundle", phobos,
       mapping_result.solution, GPLSpec::FieldEvaluation::Quantity::Potential,
       std::move(perturb_bundle));

   const auto rotation_angles =
       GPLSpec::detail::RotationOutputAnglesFor(vec_ang1, vec_ang2);
   const GPLSpec::FieldEvaluation::Rotation rotation{
       rotation_angles.alpha, rotation_angles.beta, rotation_angles.gamma};
   std::vector<double> slice_azimuths;
   for (double longitude : phobos.GSH_Grid().Longitudes())
      slice_azimuths.push_back(longitude);
   auto rotated_metadata = GPLSpec::FieldBundle::WithSolve(
       full_result, "phobos_gravity_perturbation", "full_rotated_slice", "CG",
       full_result.final_solve().requested_tolerance);
   rotated_metadata.inputs.push_back(
       GPLSpec::FieldBundle::InputPath(pathtofile));
   rotated_metadata.parameters = run_parameters;
   rotated_metadata.rotation = rotation;
   GPLSpec::FieldBundle::WriteAngular(
       std::filesystem::path(pathtofolder3) / "full_rotated_grid_bundle",
       phobos, full_result.solution,
       GPLSpec::FieldEvaluation::Quantity::Potential,
       std::move(rotated_metadata), &rotation);
   auto rotated_slice_metadata = GPLSpec::FieldBundle::WithSolve(
       full_result, "phobos_gravity_perturbation", "full_rotated_slice", "CG",
       full_result.final_solve().requested_tolerance);
   rotated_slice_metadata.inputs.push_back(
       GPLSpec::FieldBundle::InputPath(pathtofile));
   rotated_slice_metadata.parameters = run_parameters;
   rotated_slice_metadata.rotation = rotation;
   GPLSpec::FieldBundle::WriteSlice(
       std::filesystem::path(pathtofolder3) / "full_rotated_slice_bundle",
       phobos, full_result.solution,
       GPLSpec::FieldEvaluation::Quantity::Potential, rotation,
       slice_azimuths, std::move(rotated_slice_metadata));

   //    std::cout << "Value: " << valtest << "\n";
   return 0;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// // test rotation
// //  find Euler angles
// double theta1 = std::numbers::pi_v<double> / 2.0;
// double theta2 = std::numbers::pi_v<double> / 2.0;
// theta2 *= 0.5;
// double phi1 = std::numbers::pi_v<double> / 2.0;
// phi1 = 0.0;
// double phi2 = std::numbers::pi_v<double> / 2.0;
// phi2 = 0.0;

// double cosphi = std::cos(theta2) * std::cos(theta1) +
//                 std::sin(theta2) * std::sin(theta1) * std::cos(phi2 - phi1);
// double sinphi = std::sqrt(1 - cosphi * cosphi);
// auto tmp1 = std::sin(theta2) * std::cos(theta1) * std::cos(phi2) -
//             std::cos(theta2) * std::sin(theta1) * std::cos(phi1);
// auto tmp2 = std::cos(theta2) * std::sin(theta1) * std::sin(phi1) -
//             std::sin(theta2) * std::cos(theta1) * std::sin(phi2);
// auto tmp3 = std::sin(theta2) * std::sin(theta1) * std::sin(phi2 - phi1);
// auto tmp4 = std::cos(theta1) * cosphi - std::cos(theta2);
// auto tmp5 = std::cos(theta1) * sinphi;

// if (std::abs(tmp1) < std::numeric_limits<double>::epsilon()) {
//    tmp1 = 0.0;
// }
// if (std::abs(tmp2) < std::numeric_limits<double>::epsilon()) {
//    tmp2 = 0.0;
// }
// if (std::abs(tmp4) < std::numeric_limits<double>::epsilon()) {
//    tmp4 = 0.0;
// }
// if (std::abs(tmp5) < std::numeric_limits<double>::epsilon()) {
//    tmp5 = 0.0;
// }

// double alpha = std::atan2(tmp1, tmp2);
// double beta = std::acos(tmp3 / sinphi);
// double gamma = std::atan2(tmp4, tmp5);

// // std::cout << "\ntmp1: " << tmp1 << ", tmp2: " << tmp2 << ", alpha: " <<
// // alpha
// //           << "\n";

// // timer1.stop("Time to construct");

// // vectors containing theta,phi information
// std::vector<double> vec_ang1{theta1, phi1}, vec_ang2{theta2, phi2};

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
