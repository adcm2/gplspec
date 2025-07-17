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
int
main() {

   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;

   // declaring variables:
   double normrad = 1.0;
   double normrho = 1.0;
   double lengthnorm = 6371000.0;
   double timenorm = 3600.0;
   double massnorm = 5.972 * std::pow(10.0, 24.0);
   double maxstep = 0.01;
   double ballrad = 1.2;

   Timer timer1;
   timer1.start();

   // declare model and find the potential
   Density3D testsphere = Density3D::SphericalHomogeneousPlanet(
       normrad, normrho, lengthnorm, timenorm, massnorm, maxstep, ballrad);
   //    Eigen::VectorXcd vec_potsol = FindGravitationalPotential(testsphere);

   //    // convert to "standard" format
   //    auto stdvec_potsol =
   //    testsphere.SingleEigenVectorToGeneralFormat(vec_potsol);
   auto stdvec_potsol = FindGravitationalPotential(testsphere);
   timer1.stop("Time for 3D");

   // solution using spherical integration method:
   timer1.start();
   auto vec_integral_potential = GravitationalSphericalIntegral(testsphere);
   timer1.stop("Time for integral");

   auto vec_exactsol = HomogeneousSphereIntegral(testsphere);

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
   std::string pathtofolder = "./work/Bench1";
   testsphere.ReferentialOutputAtElement(pathtofolder, stdvec_potsol);
   testsphere.OutputAtElement(pathtofolder, vec_integral_potential);
   testsphere.OutputAtElement(pathtofolder, vec_exactsol);

   return 0;
}
