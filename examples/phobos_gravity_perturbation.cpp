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
   int lsmall = 0;
   std::cout << "lmax for base calculation: ";
   std::cin >> lsmall;
   int npoly = 5;
   using MRange = GSHTrans::All;
   using NRange = GSHTrans::All;
   using Grid = GSHTrans::GaussLegendreGrid<double, MRange, NRange>;
   Density3D phobosfull(scaledensity, pathtofile, npoly, lMax, lengthnorm,
                        timenorm, massnorm, 0.1, 1.5, 0);
   Density3D phobos(scaledensity, pathtofile, lsmall, npoly, lMax, lengthnorm,
                    timenorm, massnorm, 0.1, 1.5, 0);

   /////////////////////////////////////////////////////////////////////////////
   // exact

   // get gravitational field
   timer1.start();
   auto stdvec_potsol =
       FindGravitationalPotential(phobosfull, std::pow(10.0, -12.0));
   timer1.stop("Time for gravity with full model");

   timer1.start();
   auto stdvec_potsol_small =
       FindGravitationalPotential(phobos, std::pow(10.0, -6.0));
   timer1.stop("Time for gravity with small model");
   /////////////////////////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////////////////////
   // find perturbed solution
   // define map
   timer1.start();
   MappingPerturbation testpertmap(phobos, pathtofile, lsmall + 1, lMax);
   timer1.stop("Time for map to be defined");

   // get solution
   timer1.start();
   auto stdvec_potsol_perturb = FindGravitationalPotentialPerturbation(
       phobos, testpertmap, std::pow(10.0, -6.0), std::pow(10.0, -6.0));
   // auto stdvec_potsol_perturb =
   //     FindGravitationalPotential(phobosfull, std::pow(10.0, -3.0));
   timer1.stop("Time for perturbation calculation");

   /////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////
   // output test
   std::string pathtofolder1 = "./work/Phobos/Cartesian";
   std::string pathtofolder2 = "./work/Phobos/Spherical";
   std::string pathtofolder3 = "./work/Phobos/Rotated";
   std::string pathtofolder4 = "./work/Phobos/Perturbed";
   //    phobos.PhysicalOutputAtElement(pathtofolder1, stdvec_potsol);
   //    phobos.CartesianOutputAtElement(pathtofolder1, stdvec_potsol);
   phobos.PhysicalOutputAtElement(pathtofolder2, stdvec_potsol);
   //    phobos.PhysicalOutputRotated(pathtofolder3, vec_ang1, vec_ang2,
   //                                 stdvec_potsol);
   phobos.PhysicalOutputAtElement(pathtofolder4, stdvec_potsol_perturb);

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
