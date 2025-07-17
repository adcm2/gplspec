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
   double maxstep = 0.01;
   double ballrad = 1.2;
   int npoly = 5;
   int lmax = 2;
   std::string pathtoprem = "modeldata/prem.200";

   Timer timer1;
   timer1.start();

   // declare model and find the potential
   Density3D testprem = Density3D::OneDimensionalPlanetFromFile(
       pathtoprem, npoly, lmax, maxstep, ballrad);

   //    Eigen::VectorXcd vec_potsol = FindGravitationalPotential(testprem);
   //    auto stdvec_potsol =
   //    testprem.SingleEigenVectorToGeneralFormat(vec_potsol);
   auto stdvec_potsol = FindGravitationalPotential(testprem);
   timer1.stop("Time for 3D");

   // solution using spherical integration method:
   timer1.start();
   auto vec_integral_potential = GravitationalSphericalIntegral(testprem);
   timer1.stop("Time for integral");

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
   std::string pathtofolder = "./work/Bench2";
   testprem.ReferentialOutputAtElement(pathtofolder, stdvec_potsol);
   testprem.OutputAtElement(pathtofolder, vec_integral_potential);

   return 0;
}
