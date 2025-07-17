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
   double ballrad = 1.4;
   int npoly = 5;
   int lmax = 2;
   double h = 0.2;
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
      inp_map(const double h) : _h{h} {};
      auto RadialMapping(int i) const {
         auto lambdamap = [hmult = _h](double r, double theta, double phi) {
            return hmult * r;
         };
         return lambdamap;
      }

    private:
      double _h = 0.0;
   };

   // model with mapping
   inp_map map1(h);
   spherical_model inp_sphere = spherical_model::HomogeneousSphere(
       radius, density, lengthnorm, timenorm, massnorm);
   Density3D homogeneous_sphere_with_map(inp_sphere,
                                         PlanetaryModel::TomographyZeroModel(),
                                         map1, npoly, lmax, maxstep, ballrad);
   auto stdvec_potsol = FindGravitationalPotential(homogeneous_sphere_with_map,
                                                   std::pow(10.0, -12.0));
   timer1.stop("Time for 3D");
   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // benchmarking using an equivalent spherical homogeneu
   //  equivalent spherical model
   Density3D equivalent_sphere = Density3D::SphericalHomogeneousPlanet(
       radius * (1.0 + h), density / std::pow(1.0 + h, 3.0), lengthnorm,
       timenorm, massnorm, maxstep, ballrad);

   auto inp_radii = homogeneous_sphere_with_map.PhysicalRadius_Line(0);
   auto vec_exactsol = HomogeneousSphereIntegral(equivalent_sphere, inp_radii);

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
   std::string pathtofolder = "./work/Bench3";
   homogeneous_sphere_with_map.PhysicalOutputAtElement(pathtofolder,
                                                       stdvec_potsol);
   equivalent_sphere.OutputAtElement(pathtofolder, inp_radii, vec_exactsol);

   //    testprem.OutputAtElement(pathtofolder, vec_exactsol);

   return 0;
}
