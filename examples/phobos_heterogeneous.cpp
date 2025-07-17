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
   bool nhomo = true;
   Density3D phobos(scaledensity, pathtofile, npoly, lMax, lengthnorm, timenorm,
                    massnorm, 0.1, 1.5, nhomo);
   //    Density3D phobos2(scaledensity, pathtofile, npoly, 128, lengthnorm,
   //    timenorm,
   //                      massnorm, 0.1, 1.5);
   timer1.stop("Time to construct");

   int numlayers = phobos.Node_Information().NumberOfElements();
   for (int idx = 0; idx < numlayers; ++idx) {
      std::cout << phobos.Node_Information().LayerNumber(idx) << "\n";
   }

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
   // std::cout << "Volume: " << phobos.Volume() / std::pow(10.0, 9.0) << "
   // km^3"
   //           << "\n";
   // std::cout << "Density: " << phobos.Mass() / (phobos.Volume() * 1000.0)
   //           << " g/cm^3\n";

   // get gravitational field
   timer1.start();
   auto stdvec_potsol =
       FindGravitationalPotential(phobos, std::pow(10.0, -12.0));
   timer1.stop("Time for gravity");
   timer1.start();
   auto stdvec_potsol2 =
       FindGravitationalPotential(phobos, std::pow(10.0, -3.0));
   timer1.stop("Time for gravity 2");

   //    for (int idx = 0; idx < 30; ++idx) {
   //       std::cout << stdvec_potsol[10][0][idx] << "\n";
   //    }
   // get gravitational field
   timer1.start();
   auto stdvec_senskernel = SphericalHarmonicSensitivityKernel(phobos, 2, 0);
   timer1.stop("Time for sensitivity kernel");

   auto stdpower = phobos.PowerSTD(stdvec_potsol);

   // mean potential on surface
   // auto coefficientnumberpos =
   //     GSHTrans::GSHIndices<GSHTrans::NonNegative>(lMax, lMax, 0).size();
   // std::vector<std::complex<double>> vec_shpot(coefficientnumberpos, 0.0);
   // phobos.GSH_Grid().ForwardTransformation(lMax, 0, stdvec_potsol[10][0],
   //                                         vec_shpot);
   // std::cout << "Average on surface: "
   //           << stdvec_potsol[10][0][0] / sqrt(4.0 * std::numbers::pi) *
   //                  phobos.PotentialNorm()
   //           << "\n";

   /////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////
   // test rotation
   //  find Euler angles
   double theta1 = std::numbers::pi_v<double> / 2.0;
   //    theta1 = 0.001;
   double theta2 = std::numbers::pi_v<double> / 2.0;
   theta2 *= 0.0;
   double phi1 = std::numbers::pi_v<double> / 2.0;
   phi1 = 0.0;
   double phi2 = std::numbers::pi_v<double> / 2.0;
   phi2 = 0.0;

   // vectors containing theta,phi information
   std::vector<double> vec_ang1{theta1, phi1}, vec_ang2{theta2, phi2};
   std::vector<double> vec_ang12{theta1, 0.0}, vec_ang22{theta1, theta1};

   //    auto stdvec_potrot =
   //        phobos.RotateSliceToEquator(vec_ang1, vec_ang2, stdvec_potsol);

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
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
   //    phobos.PhysicalOutputAtElement(pathtofolder1, stdvec_potsol);
   phobos.CartesianOutputAtElement(pathtofolder1, stdvec_potsol);
   phobos.PhysicalOutputAtElement(pathtofolder2, stdvec_potsol);
   //   phobos.ReferentialOutputAtElement(pathtofolder3, stdvec_senskernel,
   //   true);
   phobos.ReferentialOutputSlice(pathtofolder3, stdvec_senskernel, true);
   //       phobos.PhysicalOutputSlice(pathtofolder2, stdvec_potsol);
   //    phobos.PhysicalOutputRotated(pathtofile3, vec_ang1, vec_ang2,
   //    stdvec_potsol);
   phobos.ReferentialOutputRotated(pathtofile1, vec_ang12, vec_ang22,
                                   stdvec_potsol);
   phobos.PhysicalOutputRotated(pathtofile2, vec_ang12, vec_ang22,
                                stdvec_potsol);
   phobos.ModelDensityOutputRotated(pathtofile4, vec_ang12, vec_ang22, true);
   phobos.ModelDensityOutputRotated(pathtofile5, vec_ang12, vec_ang22);
   // phobos.ReferentialOutputRotated(pathtofile2, vec_ang12,
   //    vec_ang22,
   //                                    stdvec_potsol);

   std::string pathpower = "./work/Phobos/Spherical/PowerSpectrum.out";
   auto file = std::ofstream(pathpower);

   for (int i = 0; i < stdpower.size(); ++i) {
      for (int j = 0; j < stdpower[i].size(); ++j) {
         file << std::setprecision(16) << stdpower[i][j][0];
         for (int k = 1; k < stdpower[i][j].size(); ++k) {
            file << ";" << stdpower[i][j][k];
         }
         file << std::endl;
      }
   }
   file.close();

   return 0;
}
