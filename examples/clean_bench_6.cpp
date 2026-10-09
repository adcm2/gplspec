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

int
main() {

   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;
   using namespace SimpleModels;
   using FLOAT = double;

   // declaring variables:
   double maxstep = 0.05;
   double ballrad = 1.2;
   int npoly = 12;
   //  int lmax = 1;
   int lmax;
   double h;
   double radius = 1.0;
   double density = 1.0;
   double lengthnorm = EarthModels::EarthConstants<FLOAT>().LengthNorm();
   double timenorm = EarthModels::EarthConstants<FLOAT>().TimeNorm();
   double massnorm = EarthModels::EarthConstants<FLOAT>().MassNorm();
   std::string pathtoprem = "modeldata/prem.200";
   std::string pathtotomo = "modeldata/S40RTS_dvs.nc";

   Timer timer1;

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////
   // declare model and find the potential
   h = 2;
   lmax = 20;

   // make model
   timer1.start();
   Density3D testtomo = Density3D::SphericalThreeDimensionalPlanetFromFile(
       pathtoprem, pathtotomo, npoly, lmax, maxstep, ballrad);
   timer1.stop("Time to make model");

   //////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////
   // testing preconditioner
   Eigen::SparseMatrix<std::complex<double>> testmat =
       -testtomo.SpectralElementInformation().fullmatrix<std::complex<double>>(
           testtomo.Node_Information(), 0,
           testtomo.GSH_Grid().MaxDegree());
   testmat.makeCompressed();
   //    std::cout << "\n\n" << testmat.block(0, 26, 200, 1) << "\n\n";
   //    std::cout << "\n\n" << testmat.block(0, 0, 30, 1) << "\n\n";

   //////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////

   // do PSM solution
   timer1.start();
   auto potential_result = SolvePotential(testtomo, std::pow(10.0, -14.0));
   auto stdvec_potsol = ToHarmonicField(testtomo, potential_result);
   timer1.stop("Time for 3D");

   // solution using spherical integration method:
   timer1.start();
   auto vec_integral_potential = GravitationalSphericalIntegral(testtomo);
   timer1.stop("Time for integral");

   // sensitivity kernel
   timer1.start();
   std::vector<int> vec_l{0, 1, 2, 3, 4}, vec_m{0, 0, 2, 3, 0};
   std::vector<double> multval{1.0, 3.0, 5.0, 7.0, 9.0};
   auto sensitivity_result = SolveSensitivityKernel(testtomo, vec_l, vec_m, multval);
   auto stdvec_senskernel = ToHarmonicField(testtomo, sensitivity_result);
   timer1.stop("Time for sensitivity kernel");

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
   std::string pathtofolder = "./work/Bench6";
   std::string pathtofolderlm = "./work/Bench6/lm";
   testtomo.PhysicalOutputAtElement(pathtofolder, stdvec_potsol);
   testtomo.PhysicalOutputAtElement(pathtofolder, vec_integral_potential);

   testtomo.ReferentialOutputAtElement(pathtofolderlm, stdvec_potsol);
   testtomo.OutputAtElement(pathtofolderlm, vec_integral_potential);
   testtomo.ReferentialOutputAtElement(pathtofolderlm, stdvec_senskernel, true);
   auto bundle = GPLSpec::FieldBundle::WithSolve(
       potential_result, "clean_bench_6", "prem_s40rts", "CG",
       potential_result.final_solve().requested_tolerance);
   bundle.inputs.push_back(GPLSpec::FieldBundle::InputPath(pathtoprem));
   bundle.inputs.push_back(GPLSpec::FieldBundle::InputPath(pathtotomo));
   bundle.parameters = {{"polynomial_order", std::to_string(npoly)},
       {"degree_max", std::to_string(lmax)}, {"max_step_norm", std::to_string(maxstep)},
       {"exterior_radius_norm", std::to_string(ballrad)}, {"lateral_degree", "2"},
       {"lateral_amplitude", std::to_string(h)},
       {"length_norm_m", std::to_string(lengthnorm)},
       {"mass_norm_kg", std::to_string(massnorm)},
       {"time_norm_s", std::to_string(timenorm)}};
   bundle.selector_l = vec_l;
   bundle.selector_m = vec_m;
   bundle.selector_weight = multval;
   const auto bundle_dir = std::filesystem::path(pathtofolder) / "field_bundle";
   GPLSpec::FieldBundle::Writer(bundle_dir, testtomo, bundle).Write([&](auto &&emit) {
      GPLSpec::FieldEvaluation::ForEachAngularGridSample(
          testtomo, potential_result.solution,
          GPLSpec::FieldEvaluation::Quantity::Potential, emit);
      GPLSpec::FieldEvaluation::ForEachRadialHarmonicSample(
          testtomo, std::span<const std::complex<double>>{},
          GPLSpec::FieldEvaluation::Quantity::ReferentialDensity, emit);
      GPLSpec::FieldEvaluation::ForEachRadialHarmonicSample(
          testtomo, std::span<const std::complex<double>>{},
          GPLSpec::FieldEvaluation::Quantity::PhysicalDensity, emit);
      GPLSpec::FieldEvaluation::ForEachRadialHarmonicSample(
          testtomo, potential_result.solution,
          GPLSpec::FieldEvaluation::Quantity::DegreePower, emit);
   });
   auto integral_metadata = GPLSpec::FieldBundle::RunMetadata{};
   integral_metadata.workflow = "clean_bench_6";
   integral_metadata.dataset = "prem_s40rts_spherical_integral";
   integral_metadata.inputs = bundle.inputs;
   integral_metadata.parameters = bundle.parameters;
   GPLSpec::FieldBundle::Writer(
       std::filesystem::path(pathtofolderlm) / "integral_bundle", testtomo,
       std::move(integral_metadata)).Write([&](auto &&emit) {
      GPLSpec::FieldEvaluation::ForEachElementBoundaryHarmonics(
          testtomo, vec_integral_potential,
          GPLSpec::FieldEvaluation::Quantity::Potential, emit);
   });
   auto harmonic_metadata = GPLSpec::FieldBundle::WithSolve(
       potential_result, "clean_bench_6", "referential_potential_harmonics", "CG",
       potential_result.final_solve().requested_tolerance);
   harmonic_metadata.inputs = bundle.inputs;
   harmonic_metadata.parameters = bundle.parameters;
   GPLSpec::FieldBundle::WriteRadialHarmonics(
       std::filesystem::path(pathtofolderlm) / "harmonic_bundle", testtomo,
       potential_result.solution, GPLSpec::FieldEvaluation::Quantity::Potential,
       std::move(harmonic_metadata));
   auto sensitivity_metadata = GPLSpec::FieldBundle::WithSolve(
       sensitivity_result, "clean_bench_6", "selected_sensitivity", "BiCGSTAB",
       sensitivity_result.final_solve().requested_tolerance);
   sensitivity_metadata.inputs = bundle.inputs;
   sensitivity_metadata.selector_l = vec_l;
   sensitivity_metadata.selector_m = vec_m;
   sensitivity_metadata.selector_weight = multval;
   GPLSpec::FieldBundle::WriteRadialHarmonics(
       std::filesystem::path(pathtofolderlm) / "sensitivity_bundle", testtomo,
       sensitivity_result.solution, GPLSpec::FieldEvaluation::Quantity::SensitivityKernel,
       std::move(sensitivity_metadata));

   return 0;
}
