/**
 * @file clean_bench_2.cpp
 * @brief Example demonstrating and benchmarking gravitational potential
 * calculations.
 *
 * This example compares two methods for calculating the gravitational potential
 * of a 1D planetary model (PREM):
 * 1. The particle relabelling method using radial spectral elements
 * 2. A direct spherical integration method
 *
 * The execution time for each method is measured and printed to the console.
 * The results from both methods are then saved to the 'work/Bench2' directory
 * for comparison.
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
 * @brief Main entry point for the second benchmark.
 * @return int Exit status.
 */
int
main() {

   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;

   /** @defgroup bench2_step1 Step 1: Setup Model Parameters
    *  @{
    */
   // Maximum step size for the model.
   double maxstep = 0.1;

   // Radius of the bounding sphere.
   double ballrad = 1.2;

   // Polynomial degree for spectral elements.
   int npoly = 5;

   // Maximum spherical harmonic degree.
   int lmax = 2;

   // Path to the PREM model data file.
   std::string pathtoprem = "modeldata/prem.200";

   // Timer to measure performance of calculations.
   Timer timer1;
   /** @} */

   /** @defgroup bench2_step2 Step 2: Calculate Potential with 3D
    * Spectral-Element Method
    *  @{
    */
   timer1.start();

   // Declare model and find the potential
   /// @brief 3D density model initialized from a 1D planetary file (PREM).
   Density3D testprem = Density3D::OneDimensionalPlanetFromFile(
       pathtoprem, npoly, lmax, maxstep, ballrad);

   /// @brief Gravitational potential calculated using the 3D spectral-element
   /// method.
   auto stdvec_potsol = FindGravitationalPotential(testprem);
   timer1.stop("Time for 3D");
   /** @} */

   /** @defgroup bench2_step3 Step 3: Calculate Potential with Spherical
    * Integration Method
    *  @{
    */
   timer1.start();
   /// @brief Gravitational potential calculated using direct spherical
   /// integration.
   auto vec_integral_potential = GravitationalSphericalIntegral(testprem);
   timer1.stop("Time for integral");
   /** @} */

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   /** @defgroup bench2_step4 Step 4: Output Results for Comparison
    *  @{
    */
   /// @brief Output directory for the results.
   std::string pathtofolder = "./work/Bench2";
   // Output from the 3D spectral-element method
   testprem.ReferentialOutputAtElement(pathtofolder, stdvec_potsol);
   // Output from the spherical integration method
   testprem.OutputAtElement(pathtofolder, vec_integral_potential);
   /** @} */

   return 0;
}
