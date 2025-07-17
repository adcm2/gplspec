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
   double ballrad = 1.4;
   int npoly = 8;
   int lmax = 2;
   int lmax2;
   double h;
   double radius = 1.0;
   double density = 1.0;
   double lengthnorm = EarthModels::EarthConstants<FLOAT>().LengthNorm();
   double timenorm = EarthModels::EarthConstants<FLOAT>().TimeNorm();
   double massnorm = EarthModels::EarthConstants<FLOAT>().MassNorm();

   Timer timer1;

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////
   // declare model and find the potential
   // mapping class
   class inp_map {
    public:
      inp_map() {};
      inp_map(const double h, const double r) : _h{h}, _prad{r} {};
      auto RadialMapping(int i) const {
         auto lambdamap = [hmult = _h, pr = _prad](double r, double theta,
                                                   double phi) {
            double rscale = r / pr;
            // return hmult * std::sin(theta) * std::sin(phi) * r * (1.0 -
            // rscale);
            return hmult * std::sin(theta) * std::cos(theta) * std::sin(phi) *
                   r * (1.0 - rscale);
         };
         return lambdamap;
      }

    private:
      double _h = 0.0;
      double _prad = 1.0;
   };

   h = 0.2;
   lmax2 = 10;
   //  std::cout << "Type in h: ";
   //  std::cin >> h;
   //  std::cout << "Type in lmax: ";
   //  std::cin >> lmax2;
   std::cout << "Type in maxstep: ";
   std::cin >> maxstep;

   timer1.start();
   // model with mapping
   inp_map map1(h, radius);
   //  spherical_model inp_sphere = spherical_model::HomogeneousSphere(
   //      radius, density, lengthnorm, timenorm, massnorm);
   Density3D testconstruct(radius, density, radius, map1, npoly, lmax2,
                           lengthnorm, timenorm, massnorm, maxstep, ballrad);
   std::cout << testconstruct.Node_Information().NumberOfLayers() << "\n";
   auto stdvec_potsol =
       FindGravitationalPotential(testconstruct, std::pow(10.0, -12.0));
   timer1.stop("Time for 3D");

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // benchmarking using an equivalent spherical homogeneous
   //  equivalent spherical model
   //  std::cout << "\n\n HELLO HELLO \n\n";
   Density3D equivalent_sphere = Density3D::SphericalHomogeneousPlanet(
       radius, density, lengthnorm, timenorm, massnorm, maxstep, ballrad);
   //  std::cout << "\n\n HELLO HELLO \n\n";
   std::vector<std::vector<double>> inp_radii2, vec_exactsol2;
   for (int idx = 0; idx < testconstruct.GSH_Grid().NumberOfCoLatitudes() *
                               testconstruct.GSH_Grid().NumberOfLongitudes();
        ++idx) {
      auto inp_radii = testconstruct.PhysicalRadius_Line(idx);
      auto vec_exactsol =
          HomogeneousSphereIntegral(equivalent_sphere, inp_radii);
      inp_radii2.push_back(inp_radii);
      vec_exactsol2.push_back(vec_exactsol);
   }
   //  for (int idx = 0; idx < vec_exactsol2[0].size(); ++idx) {
   //     std::cout << std::setprecision(16) << "idx: " << vec_exactsol2[0][idx]
   //               << "\n";
   //  }
   std::cout << std::setprecision(16)
             << -2.0 * 3.1415926535897932 *
                    equivalent_sphere.GravitationalConstant() * density *
                    radius * radius
             << "\n";
   std::cout << std::setprecision(16)
             << stdvec_potsol[0][0][0] / (std::sqrt(4.0 * 3.1415926535897932))
             << "\n";
   std::cout << std::setprecision(16)
             << -2.0 * 3.1415926535897932 *
                    equivalent_sphere.GravitationalConstant() * density *
                    radius * radius * equivalent_sphere.PotentialNorm()
             << "\n";

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
   std::string pathtofolder = "./work/Bench5";
   testconstruct.PhysicalOutputAtElement(pathtofolder, stdvec_potsol);
   equivalent_sphere.OutputAtElement(pathtofolder, inp_radii2, vec_exactsol2);

   //    testprem.OutputAtElement(pathtofolder, vec_exactsol);

   return 0;
}
