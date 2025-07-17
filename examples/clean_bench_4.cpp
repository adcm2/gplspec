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
   using namespace SimpleModels;
   using FLOAT = double;

   // declaring variables:
   double maxstep = 0.1;
   double ballrad = 1.2;
   int npoly = 5;
   int lmax = 2;
   double h = 0.02;
   double radius = 1.0;
   double density = 1.0;
   double lengthnorm = EarthModels::EarthConstants<FLOAT>().LengthNorm();
   double timenorm = EarthModels::EarthConstants<FLOAT>().TimeNorm();
   double massnorm = EarthModels::EarthConstants<FLOAT>().MassNorm();
   std::string pathtoprem = "modeldata/prem.200";

   Timer timer1;
   timer1.start();

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
            return hmult * r * (1.0 - rscale);
         };
         return lambdamap;
      }

    private:
      double _h = 0.0;
      double _prad = 1.0;
   };

   // model with mapping
   inp_map map1(h, radius);
   spherical_model inp_sphere = spherical_model::HomogeneousSphere(
       radius, density, lengthnorm, timenorm, massnorm);
   Density3D homogeneous_sphere_with_map(inp_sphere,
                                         PlanetaryModel::TomographyZeroModel(),
                                         map1, npoly, lmax, maxstep, ballrad);
   Density3D testconstruct(radius, density, radius, map1, 5, 2, lengthnorm,
                           timenorm, massnorm, maxstep, ballrad);
   auto stdvec_potsol =
       FindGravitationalPotential(testconstruct, std::pow(10.0, -12.0));
   timer1.stop("Time for 3D");

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // benchmarking using an equivalent spherical homogeneous
   //  equivalent spherical model
   Density3D equivalent_sphere = Density3D::SphericalHomogeneousPlanet(
       radius, density, lengthnorm, timenorm, massnorm, maxstep, ballrad);

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

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
   std::string pathtofolder = "./work/Bench4";
   testconstruct.PhysicalOutputAtElement(pathtofolder, stdvec_potsol);
   equivalent_sphere.OutputAtElement(pathtofolder, inp_radii2, vec_exactsol2);

   //    testprem.OutputAtElement(pathtofolder, vec_exactsol);

   return 0;
}
