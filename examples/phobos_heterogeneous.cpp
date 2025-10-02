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
   TestTools::PhobosRead(pathtofile);

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
   Density3D phobos(scaledensity, pathtofile, npoly, lMax, lengthnorm, timenorm,
                    massnorm, 0.1, 1.5, nhomo);
   timer1.stop("Time to construct");

   // Compute gravitational potential solutions
   timer1.start();
   auto stdvec_potsol =
       FindGravitationalPotential(phobos, std::pow(10.0, -5.0));
   timer1.stop("Time for gravity");

   timer1.start();
   auto stdvec_potsol2 =
       FindGravitationalPotential(phobos, std::pow(10.0, -3.0));
   timer1.stop("Time for gravity 2");

   // Compute spherical harmonic sensitivity kernel
   timer1.start();
   auto stdvec_senskernel = SphericalHarmonicSensitivityKernel(phobos, 2, 0);
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
   file.close();

   return 0;
}
