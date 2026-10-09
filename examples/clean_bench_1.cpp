/**
 * @file clean_bench_1.cpp
 * @brief Benchmark for gravitational potential calculation of a homogeneous
 * sphere.
 *
 * This example demonstrates:
 * - Construction of a homogeneous spherical density model.
 * - Calculation of gravitational potential using two methods.
 * - Output of results for further analysis.
 */

#include <gplspec/All>

/**
 * @brief Main routine for the homogeneous sphere gravity benchmark.
 *
 * - Constructs a homogeneous sphere model.
 * - Computes gravitational potential using direct and integral methods.
 * - Outputs results to files.
 */
int
main() {
   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;

   // Physical and model parameters
   // Normalised radius, density, length, time and mass normalisations, the
   // maximum mesh step size and the ball radius
   double nrad = 1.0, nrho = 1.0, lengthnorm = 6371000.0, timenorm = 3600.0,
          massnorm = 5.972e24, maxstep = 0.1, ballrad = 1.2;

   Timer timer1;

   // Construct homogeneous sphere model and compute potential
   timer1.start();
   Density3D testsphere = Density3D::SphericalHomogeneousPlanet(
       nrad, nrho, lengthnorm, timenorm, massnorm, maxstep, ballrad);
   auto potential_result = SolvePotential(testsphere);
   auto stdvec_potsol = ToHarmonicField(testsphere, potential_result);
   timer1.stop("Time for 3D");

   // Compute potential using spherical integration method
   timer1.start();
   auto vec_integral_potential = GravitationalSphericalIntegral(testsphere);
   timer1.stop("Time for integral");

   // Compute exact solution for homogeneous sphere
   auto vec_exactsol = HomogeneousSphereIntegral(testsphere);

   // Output results to files
   std::string pathtofolder = "./work/Bench1";
   testsphere.ReferentialOutputAtElement(pathtofolder, stdvec_potsol);
   testsphere.OutputAtElement(pathtofolder, vec_integral_potential);
   testsphere.OutputAtElement(pathtofolder, vec_exactsol);
   auto bundle = GPLSpec::FieldBundle::WithSolve(
       potential_result, "clean_bench_1", "homogeneous_sphere", "CG",
       potential_result.final_solve().requested_tolerance);
   GPLSpec::FieldBundle::Writer(std::filesystem::path(pathtofolder) / "field_bundle",
                                testsphere, bundle).Write([&](auto &&emit) {
      GPLSpec::FieldEvaluation::ForEachAngularGridSample(
          testsphere, potential_result.solution,
          GPLSpec::FieldEvaluation::Quantity::Potential, emit);
      GPLSpec::FieldEvaluation::ForEachRadialHarmonicSample(
          testsphere, potential_result.solution,
          GPLSpec::FieldEvaluation::Quantity::Potential, emit);
      GPLSpec::FieldEvaluation::ForEachAngularGridSample(
          testsphere, std::span<const std::complex<double>>{},
          GPLSpec::FieldEvaluation::Quantity::ReferentialDensity, emit);
   });
   auto integral_metadata = GPLSpec::FieldBundle::RunMetadata{};
   integral_metadata.workflow = "clean_bench_1";
   integral_metadata.dataset = "spherical_integral_potential";
   integral_metadata.parameters = bundle.parameters;
   GPLSpec::FieldBundle::Writer(
       std::filesystem::path(pathtofolder) / "integral_bundle", testsphere,
       std::move(integral_metadata)).Write([&](auto &&emit) {
      GPLSpec::FieldEvaluation::ForEachElementBoundaryHarmonics(
          testsphere, vec_integral_potential,
          GPLSpec::FieldEvaluation::Quantity::Potential, emit);
   });
   auto exact_metadata = GPLSpec::FieldBundle::RunMetadata{};
   exact_metadata.workflow = "clean_bench_1";
   exact_metadata.dataset = "analytic_homogeneous_sphere";
   exact_metadata.parameters = bundle.parameters;
   GPLSpec::FieldBundle::Writer(
       std::filesystem::path(pathtofolder) / "exact_bundle", testsphere,
       std::move(exact_metadata)).Write([&](auto &&emit) {
      GPLSpec::FieldEvaluation::ForEachElementBoundaryMonopole(
          testsphere, vec_exactsol, GPLSpec::FieldEvaluation::Quantity::Potential,
          emit);
   });

   return 0;
}
