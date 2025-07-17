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
            // return hmult * std::sin(theta) * std::sin(phi) * r * (1.0 -
            // rscale);
            // return hmult * std::sin(theta) * std::cos(theta) * std::sin(phi)
            // *
            //        r * (1.0 - rscale);
            // return hmult * std::cos(theta) * r * (1.0 - rscale);
            return hmult * r;
         };
         return lambdamap;
      }

    private:
      double _h = 0.0;
      double _prad = 1.0;
   };

   h = 0.0;
   lmax = 20;

   std::cout << "Type in h: ";
   std::cin >> h;
   inp_map map1(h, 1.0);
   inp_map map0(0.0, 1.0);
   // std::cout << "Type in lmax: ";
   // std::cin >> lmax2;

   // make model
   timer1.start();
   // Density3D testtomo = Density3D::SphericalThreeDimensionalPlanetFromFile(
   //     pathtoprem, pathtotomo, npoly, lmax, maxstep, ballrad);
   // Density3D testtomo = Density3D::SphericalHomogeneousPlanet(
   //     normrad, normrho, lengthnorm, timenorm, massnorm, maxstep, ballrad);
   SimpleModels::spherical_model _sphere1D =
       SimpleModels::spherical_model::HomogeneousSphere(
           normrad, normrho, lengthnorm, timenorm, massnorm);

   // model with no map
   Density3D testtomo(_sphere1D, PlanetaryModel::TomographyZeroModel(), map0,
                      npoly, lmax, maxstep, ballrad);

   // model with map
   Density3D testtomo2(_sphere1D, PlanetaryModel::TomographyZeroModel(), map1,
                       npoly, lmax, maxstep, ballrad);
   timer1.stop("Time to make model");

   // do PSM solution
   timer1.start();
   auto stdvec_potsol =
       FindGravitationalPotential(testtomo, std::pow(10.0, -12.0));
   auto stdvec_potsol_exact =
       FindGravitationalPotential(testtomo2, std::pow(10.0, -12.0));
   timer1.stop("Time for 3D");

   // solution using spherical integration method:
   timer1.start();
   auto vec_integral_potential = GravitationalSphericalIntegral(testtomo);
   timer1.stop("Time for integral");

   //    // find perturbation force associated with BP
   //    Eigen::VectorXcd testforce = FindBoundaryPerturbationForce(testtomo,
   //    map1);

   //    // find perturbation associated with ``advection''
   //    Eigen::VectorXcd testchange =
   //        AdvectiveBoundaryPerturbation(testtomo, map1, testforce);

   // testing perturbation class
   // std::cout << "HELLO\n";

   // find perturbed solution
   MappingPerturbation testpertmap(testtomo, map1);
   auto stdvec_potsol1 = FindGravitationalPotentialPerturbation(
       testtomo, testpertmap, std::pow(10.0, -12.0));

   // variation with h:
   timer1.start();
   int lenh = 100;
   std::vector<double> vec_h(lenh);
   std::vector<double> vec_zerosol(lenh), vec_exactsol(lenh);
   for (int idx = 0; idx < lenh; ++idx) {
      vec_h[idx] = 0.2 / lenh * idx;
   }
   {
      double pi_db = 3.1415926535897932;
      double nopertphi = -2.0 * testtomo.GravitationalConstant() * pi_db *
                         normrho * std::pow(normrad, 2.0);
      for (int idx = 0; idx < lenh; ++idx) {
         inp_map mapidx(vec_h[idx], normrad);
         MappingPerturbation pertmapidx(testtomo, mapidx);
         auto stdvec_potsol_idx = FindGravitationalPotentialPerturbation(
             testtomo, pertmapidx, std::pow(10.0, -12.0));
         vec_zerosol[idx] = stdvec_potsol_idx[0][0][0].real() /
                            sqrt(4.0 * pi_db) * testtomo.PotentialNorm();
         vec_exactsol[idx] =
             nopertphi / (1.0 + vec_h[idx]) * testtomo.PotentialNorm();
      }
   }
   timer1.stop("Time for multiple solutions");

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
   std::string pathtofolder1 = "./work/Bench7/exact";
   std::string pathtofolder2 = "./work/Bench7/perturbed";
   std::string pathtofolder3 = "./work/Bench7/unperturbed";
   testtomo.PhysicalOutputAtElement(pathtofolder1, stdvec_potsol_exact);
   testtomo.PhysicalOutputAtElement(pathtofolder2, stdvec_potsol1);
   testtomo.PhysicalOutputAtElement(pathtofolder3, stdvec_potsol);

   std::string pathtofile = "./work/Bench7/zeroradius/MatrixSolution.out";
   auto file2 = std::ofstream(pathtofile);
   for (int idx = 0; idx < lenh; ++idx) {
      file2 << vec_h[idx] << ";" << vec_zerosol[idx] << ";" << vec_exactsol[idx]
            << std::endl;
   }
   file2.close();
   // testtomo.ReferentialOutputAtElement(pathtofolderlm, stdvec_potsol);
   // testtomo.OutputAtElement(pathtofolderlm, vec_integral_potential);
   //  equivalent_sphere.OutputAtElement(pathtofolder, inp_radii2,
   //  vec_exactsol2);

   //    testprem.OutputAtElement(pathtofolder, vec_exactsol);

   return 0;
}
