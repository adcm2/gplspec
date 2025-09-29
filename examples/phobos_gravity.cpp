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
   TestTools::PhobosRead(pathtofile);
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
   Density3D phobos(scaledensity, pathtofile, npoly, lMax, lengthnorm, timenorm,
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
   auto stdvec_potsol =
       FindGravitationalPotential(phobos, std::pow(10.0, -12.0));
   timer1.stop("Time for gravity");
   timer1.start();
   auto stdvec_potsol2 =
       FindGravitationalPotential(phobos, std::pow(10.0, -3.0));
   timer1.stop("Time for gravity 2");

   // get gravitational field
   timer1.start();
   auto stdvec_senskernel = SphericalHarmonicSensitivityKernel(phobos, 2, 0);
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
   std::string pathtofolder1 = "./work/Phobos/Cartesian";
   std::string pathtofolder2 = "./work/Phobos/Spherical";
   std::string pathtofolder3 = "./work/Phobos/Sensitivity";
   std::string pathtofile1 =
       pathtofolder2 + "/MatrixSolutionReferentialRotated.out";
   std::string pathtofile2 = pathtofolder2 + "/MatrixSolutionPhysical.out";
   std::string pathtofile3 = pathtofolder2 + "/MatrixSolutionRotated.out";
   std::string pathtofile4 =
       pathtofolder2 + "/PhysicalDensitySolutionRotated.out";
   std::string pathtofile5 =
       pathtofolder2 + "/ReferentialDensitySolutionRotated.out";

   timer1.start();
   phobos.CartesianOutputAtElement(pathtofolder1, stdvec_potsol);
   timer1.stop("First output");

   timer1.start();
   phobos.PhysicalOutputAtElement(pathtofolder2, stdvec_potsol);
   timer1.stop("Second output");

   timer1.start();
   phobos.ReferentialOutputAtElement(pathtofolder3, stdvec_senskernel, true);
   timer1.stop("Third output");

   timer1.start();
   phobos.ReferentialOutputSlice(pathtofolder3, stdvec_senskernel, true);
   timer1.stop("Fifth output");

   timer1.start();
   phobos.PhysicalOutputSlice(pathtofolder2, stdvec_potsol);
   timer1.stop("Sixth output");

   timer1.start();
   phobos.ReferentialOutputRotated(pathtofile1, vec_ang1, vec_ang2,
                                   stdvec_potsol);
   timer1.stop("Referential rotated");

   timer1.start();
   phobos.PhysicalOutputRotated(pathtofile2, vec_ang12, vec_ang22,
                                stdvec_potsol);
   timer1.stop("Physical rotated");

   return 0;
}
