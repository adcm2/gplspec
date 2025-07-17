#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <PlanetaryModel/All>

#include "Spherical_Integrator.h"
#include "Timer_Class.h"
#include <FFTWpp/Ranges>
#include <GSHTrans/All>
#include <TomographyModels/All>
#include <algorithm>
#include <chrono>
#include <complex>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>
// #include "Pseudospectral_Matrix_Wrapper.h"

double
azimuthtolongitude(double phi) {
   return 180.0 / 3.1415926535 * phi;
}
double
polartolatitude(double theta) {
   return 90.0 - (180.0 / 3.1415926535 * theta);
}

int
main() {
   using namespace std::chrono;
   using namespace Interpolation;
   using namespace GSHTrans;
   using namespace EarthModels;
   using namespace GravityFunctions;
   using Real = double;
   using Complex = std::complex<Real>;
   using Scalar = Complex;
   using MRange = All;
   using NRange = NonNegative;
   using Grid = GaussLegendreGrid<Real, MRange, NRange>;

   auto lMax = 4;                  // max L
   auto nMax = 0;                  // upper index
   auto grid = Grid(lMax, nMax);   // grid for transformations
   auto size = GSHIndices<All>(lMax, lMax, nMax).size();   // dof

   // polynomial order
   int npoly = 6;

   //////////////////////////////////////////////////////////////////

   // background 1D model
   auto myprem = PREM();

   // number of layers and scale factor
   int nlayer = myprem.NumberOfLayers();
   double rad_scale = myprem.OuterRadius();

   // find the nodes
   std::vector<double> vec_noderadii = Radial_Node(myprem, 100000.0);

   int nelem = vec_noderadii.size() - 1;   // total number of elements
   int matlen = nelem * npoly + 1;         // size of matrix

   // function:
   int N = 5;
   std::vector<double> x(N), y(N);

   // finite element matrices
   Eigen::VectorXd vecsol = Eigen::VectorXd::Zero(matlen);
   Eigen::SparseMatrix<double> matspec(matlen, matlen);

   // generate Gauss grid
   auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoly + 1);
   auto vec_allradii = GravityFunctions::All_Node(vec_noderadii, q);
   std::vector<double> vec_alldepths =
       GravityFunctions::DepthRadiusSwap(vec_allradii);

   // change units to metres
   RescaleVector(vec_allradii, myprem.OuterRadius());
   RescaleVector(vec_alldepths, myprem.OuterRadius());

   // values of G, pi and the scales relevant
   const double bigg_db = 6.6743 * std::pow(10.0, -11.0);
   const double pi_db = 3.1415926535;
   double scdiff = rad_scale / myprem.OuterRadius();
   double multfact = 4.0 * pi_db * bigg_db *
                     std::pow(rad_scale, 2.0);   // multiplication factor
   int laynum = 0;                               // initialise layer number

   //////////////////////////////////////////////////////////////////
   // Tomography model
   Timer timer1, timer2;
   timer1.start();
   std::filesystem::path cwd = std::filesystem::current_path();
   std::string fwd = cwd / "modeldata/S40RTS_dvs.nc";
   auto tomo = Tomography(fwd);
   timer1.stop("Time to initialise tomography model");

   // Density decomposition:
   auto vec_denslm = TomographyModelTransform(
       myprem, tomo, grid, CoefficientOrdering::RadialClumped, q, vec_noderadii,
       vec_alldepths);

   //////////////////////////////////////////////////////////////////
   // getting the force vector
   auto vec_fullforce = FindForce(myprem, CoefficientOrdering::RadialClumped, q,
                                  vec_noderadii, vec_denslm, lMax);

   //////////////////////////////////////////////////////////////////
   // setting up big matrix and solving for the force vector
   timer2.start();
   SphericalGeometryPoissonSolver<double, EarthModels::PREM> mybigtest =
       SphericalGeometryPoissonSolver<double, EarthModels::PREM>(
           myprem, CoefficientOrdering::RadialClumped, vec_noderadii, npoly,
           lMax);
   Eigen::VectorXcd mybigsol = mybigtest.solve(vec_fullforce);
   timer2.stop("Time taken to set up and solve using my class");

   // converting to grid, ie doing inverse transformation:
   auto vec_gridsol = SphdecompToGridSol(mybigsol, grid, nelem, npoly, size);

   auto vec_premsol = PotentialSolver1D(myprem, q, vec_noderadii);

   // derivative from 1D
   auto vecderiv = Gravity1D(q, vecsol, vec_noderadii, nelem, npoly, rad_scale);

   //////////////////////////////////////////////////////////////////
   // outputting result at a particular
   std::string pathtofile = "./work/BigSol.out";
   auto file2 = std::ofstream(pathtofile);
   for (int i = 0; i < nelem + 1; ++i) {
      file2 << std::setprecision(16) << vec_noderadii[i] << ";"
            << vec_premsol[0][i];
      for (int idxl1 = 0; idxl1 < lMax + 1; ++idxl1) {
         for (int idxl2 = 0; idxl2 < lMax + 1; ++idxl2) {
            file2 << ";" << vec_gridsol[i * npoly](idxl1, idxl2).real() << ";"
                  << vec_gridsol[i * npoly](idxl1, idxl2).imag();
         }
      }
      file2 << std::endl;
      // file2 << std::endl;
   };
}
