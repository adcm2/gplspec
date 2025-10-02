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
   double normrad = 1.0;            ///< Normalized radius
   double normrho = 1.0;            ///< Normalized density
   double lengthnorm = 6371000.0;   ///< Length normalization (meters)
   double timenorm = 3600.0;        ///< Time normalization (seconds)
   double massnorm = 5.972e24;      ///< Mass normalization (kg)
   double maxstep = 0.01;           ///< Maximum step size for mesh
   double ballrad = 1.2;            ///< Ball radius for computation

   Timer timer1;

   // Construct homogeneous sphere model and compute potential
   timer1.start();
   Density3D testsphere = Density3D::SphericalHomogeneousPlanet(
       normrad, normrho, lengthnorm, timenorm, massnorm, maxstep, ballrad);
   auto stdvec_potsol = FindGravitationalPotential(testsphere);
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

   return 0;
}
