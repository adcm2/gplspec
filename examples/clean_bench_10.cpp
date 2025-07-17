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
   double maxstep = 0.1;
   double ballrad = 1.5;
   int npoly = 5;
   //  int lmax = 1;
   double normrad = 1.0;
   double normrho = 1.0;
   int lmax;
   double h;
   double radius = 1.0;
   double density = 1.0;
   double lengthnorm = EarthModels::EarthConstants<FLOAT>().LengthNorm();
   double timenorm = EarthModels::EarthConstants<FLOAT>().TimeNorm();
   double massnorm = EarthModels::EarthConstants<FLOAT>().MassNorm();
   std::string pathtoprem = "modeldata/prem.200";
   std::string pathtotomo = "modeldata/S40RTS_dvs.nc";

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
            // return hmult * std::cos(theta) * r * (1.0 - rscale);
            return hmult * std::cos(theta) * r;
         };
         return lambdamap;
      }

    private:
      double _h = 0.0;
      double _prad = 1.0;
   };

   h = 0.0;
   lmax = 30;

   std::cout << "Type in h: ";
   std::cin >> h;
   //    std::cout << "Type in lmax: ";
   //    std::cin >> lmax;
   inp_map map0(0.0, 1.0);
   inp_map map1(h, 1.0);

   // make model
   timer1.start();
   SimpleModels::spherical_model _sphere1D =
       SimpleModels::spherical_model::HomogeneousSphere(
           normrad, normrho, lengthnorm, timenorm, massnorm);
   Density3D testtomo(_sphere1D, PlanetaryModel::TomographyZeroModel(), map0,
                      npoly, lmax, maxstep, ballrad);
   Density3D testtomo2(_sphere1D, PlanetaryModel::TomographyZeroModel(), map1,
                       npoly, lmax, maxstep, ballrad);
   timer1.stop("Time to make model");

   //    // output grid information
   //    std::cout << "lmax: " << testtomo.GSH_Grid().MaxDegree()
   //              << ", theta: " << testtomo.GSH_Grid().NumberOfCoLatitudes()
   //              << ", phi: " << testtomo.GSH_Grid().NumberOfLongitudes() <<
   //              "\n\n";
   //    for (auto idx : testtomo.GSH_Grid().CoLatitudes()) {
   //       std::cout << idx << "\n";
   //    }
   //    double rand3 = 0.0;
   // length of coefficients for YLM

   //    auto coefficientnumber =
   //        GSHTrans::GSHIndices<GSHTrans::All>(lmax, lmax, 0).size();
   //    auto rho_spatial =
   //        std::vector<std::complex<double>>(2 * (lmax + 1) * lmax, 0.0);
   //    auto rho_ylm = std::vector<std::complex<double>>(coefficientnumber,
   //    0.0); timer1.start(); for (int idx = 0; idx < 100; ++idx) {
   //       testtomo.GSH_Grid().ForwardTransformation(lmax, 0, rho_spatial,
   //       rho_ylm);
   //    }
   //    timer1.stop("Time to forward transform");
   //    timer1.start();
   //    for (int idx = 0; idx < 100; ++idx) {
   //       testtomo.GSH_Grid().InverseTransformation(lmax, 0, rho_ylm,
   //       rho_spatial);
   //    }
   //    timer1.stop("Time to inverse transform");

   // do PSM solution
   timer1.start();
   auto stdvec_potsol =
       FindGravitationalPotential(testtomo, std::pow(10.0, -12.0));
   auto stdvec_potsol_exact =
       FindGravitationalPotential(testtomo2, std::pow(10.0, -12.0));
   timer1.stop("Time for 3D");

   // find solution with boundary perturbation theory
   auto stdvec_potsol_perturb = FindGravitationalPotentialClassicalPerturbation(
       testtomo, map1, std::pow(10.0, -12.0));

   // find perturbed solution
   MappingPerturbation testpertmap(testtomo, map1);
   auto stdvec_potsol_perturb2 = FindGravitationalPotentialPerturbation(
       testtomo, testpertmap, std::pow(10.0, -12.0));

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
   std::string pathtofolder1 = "./work/Bench10/exact";
   std::string pathtofolder2 = "./work/Bench10/perturbed";
   std::string pathtofolder3 = "./work/Bench10/unperturbed";
   std::string pathtofolder4 = "./work/Bench10/perturbed2";
   testtomo.PhysicalOutputAtElement(pathtofolder1, stdvec_potsol_exact);
   testtomo.PhysicalOutputAtElement(pathtofolder2, stdvec_potsol_perturb);
   testtomo.PhysicalOutputAtElement(pathtofolder3, stdvec_potsol);
   testtomo.PhysicalOutputAtElement(pathtofolder4, stdvec_potsol_perturb2);

   return 0;
}
