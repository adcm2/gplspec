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
   double ballrad = 1.2;
   int npoly = 12;
   //  int lmax = 1;
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
   // class inp_map {
   //  public:
   //    inp_map() {};
   //    inp_map(const double h, const double r) : _h{h}, _prad{r} {};
   //    auto RadialMapping(int i) const {
   //       auto lambdamap = [hmult = _h, pr = _prad](double r, double theta,
   //                                                 double phi) {
   //          double rscale = r / pr;
   //          // return hmult * std::sin(theta) * std::sin(phi) * r * (1.0 -
   //          // rscale);
   //          return hmult * std::sin(theta) * std::cos(theta) * std::sin(phi)
   //          *
   //                 r * (1.0 - rscale);
   //       };
   //       return lambdamap;
   //    }

   //  private:
   //    double _h = 0.0;
   //    double _prad = 1.0;
   // };

   h = 2;
   lmax = 20;

   // make model
   timer1.start();
   Density3D testtomo = Density3D::SphericalThreeDimensionalPlanetFromFile(
       pathtoprem, pathtotomo, npoly, lmax, maxstep, ballrad);
   timer1.stop("Time to make model");

   //////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////
   // testing preconditioner
   Eigen::SparseMatrix<std::complex<double>> testmat =
       -testtomo.SpectralElementInformation().fullmatrix<std::complex<double>>(
           0, testtomo.GSH_Grid().MaxDegree());
   testmat.makeCompressed();
   //    std::cout << "\n\n" << testmat.block(0, 26, 200, 1) << "\n\n";
   //    std::cout << "\n\n" << testmat.block(0, 0, 30, 1) << "\n\n";

   //////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////

   // do PSM solution
   timer1.start();
   auto stdvec_potsol =
       FindGravitationalPotential(testtomo, std::pow(10.0, -14.0));
   timer1.stop("Time for 3D");

   //    std::cout << "lm components at 0\n\n";
   //    int mycount = 0;
   //    for (auto idx : stdvec_potsol[0][0]) {
   //       ++mycount;
   //       if (std::abs(idx) > std::pow(10.0, -15.0)) {
   //          std::cout << "Damn: " << mycount << " " << idx << "\n";
   //       }
   //    }
   //    std::cout << "\n\n";
   // solution using spherical integration method:
   timer1.start();
   auto vec_integral_potential = GravitationalSphericalIntegral(testtomo);
   timer1.stop("Time for integral");

   // sensitivity kernel
   timer1.start();
   std::vector<int> vec_l{0, 1, 2, 3, 4}, vec_m{0, 0, 2, 3, 0};
   std::vector<double> multval{1.0, 3.0, 5.0, 7.0, 9.0};
   //    auto stdvec_senskernel = SphericalHarmonicSensitivityKernel(testtomo,
   //    2, 2);
   auto stdvec_senskernel =
       SphericalHarmonicSensitivityKernel(testtomo, vec_l, vec_m, multval);
   timer1.stop("Time for sensitivity kernel");

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
   std::string pathtofolder = "./work/Bench6";
   std::string pathtofolderlm = "./work/Bench6/lm";
   testtomo.PhysicalOutputAtElement(pathtofolder, stdvec_potsol);
   testtomo.PhysicalOutputAtElement(pathtofolder, vec_integral_potential);

   testtomo.ReferentialOutputAtElement(pathtofolderlm, stdvec_potsol);
   testtomo.OutputAtElement(pathtofolderlm, vec_integral_potential);
   testtomo.ReferentialOutputAtElement(pathtofolderlm, stdvec_senskernel, true);
   //  equivalent_sphere.OutputAtElement(pathtofolder, inp_radii2,
   //  vec_exactsol2);

   //    testprem.OutputAtElement(pathtofolder, vec_exactsol);

   return 0;
}
