#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <PlanetaryModel/All>

#include "Pseudospectral_Matrix_Wrapper.h"
#include "SphericalGeometryPreconditioner.h"
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
   using NRange = All;
   using Grid = GaussLegendreGrid<Real, MRange, NRange>;
   using MATRIX = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;

   auto lMax = 10;                 // max L
   auto nMax = 1;                  // upper index
   auto grid = Grid(lMax, nMax);   // grid for transformations
   auto size = GSHIndices<All>(lMax, lMax, 0).size();   // dof
   auto _size0 = GSHIndices<All>(lMax, lMax, 0).size();
   auto _sizepm = GSHIndices<All>(lMax, lMax, 1).size();

   // polynomial order
   int npoly = 5;

   //////////////////////////////////////////////////////////////////

   // background 1D model
   auto myprem = PREM();
   auto mypertprem = PERTPREM();

   // number of layers and scale factor
   int nlayer = myprem.NumberOfLayers();
   double rad_scale = myprem.OuterRadius();

   // find the nodes
   std::vector<double> vec_noderadii = Radial_Node(myprem, 300000.0);

   int nelem = vec_noderadii.size() - 1;   // total number of elements
   int matlen = nelem * npoly + 1;         // size of matrix
   std::vector<double> vec_elemwidth(nelem);
   for (int idx = 0; idx < nelem; ++idx) {
      vec_elemwidth[idx] = vec_noderadii[idx + 1] - vec_noderadii[idx];
   }

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

   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _mat_gaussderiv;
   _mat_gaussderiv.resize(npoly + 1, npoly + 1);
   {
      auto pleg = Interpolation::LagrangePolynomial(q.Points().begin(),
                                                    q.Points().end());
      for (int idxi = 0; idxi < npoly + 1; ++idxi) {
         for (int idxj = 0; idxj < npoly + 1; ++idxj) {
            _mat_gaussderiv(idxi, idxj) = pleg.Derivative(idxi, q.X(idxj));
         }
      }
   }

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
   mybigtest.compute();
   Eigen::VectorXcd mybigsol = mybigtest.solve(vec_fullforce);
   timer2.stop("Time taken to set up and solve using my class");

   // converting to grid, ie doing inverse transformation:
   auto vec_gridsol = SphdecompToGridSol(mybigsol, grid, nelem, npoly, size);

   auto vec_premsol = PotentialSolver1D(myprem, q, vec_noderadii);

   // derivative from 1D
   auto vecderiv = Gravity1D(q, vecsol, vec_noderadii, nelem, npoly, rad_scale);

   //////////////////////////////////////////////////////////////////
   // solving with bicgstab, preconditioning and my matrix replacement class
   // first finding the value of h on grid at each of the radial nodes
   auto intsize = grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();
   std::cout << "intsize: " << intsize
             << ", (l + 1)^2: " << std::pow((lMax + 1), 2)
             << ", size0: " << _size0 << std::endl;
   std::vector<std::vector<std::complex<double>>> vec_h(matlen);
   auto mytheta = grid.CoLatitudes();
   auto myphi = grid.Longitudes();

   // for (int i = 0; i < lMax + 1; ++i) {
   //    std::cout << "i-th value " << mypoints[i] << std::endl;
   // }
   // get spatial values of h
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {

      int idxpolymax = npoly + 1;
      if (idxelem < nelem - 1) {
         idxpolymax = npoly;
      } else {
         idxpolymax = npoly + 1;
      }
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         std::vector<std::complex<double>> tmp_h(intsize);
         auto radr = GravityFunctions::StandardIntervalMap(
             q.X(idxpoly), vec_noderadii[idxelem], vec_noderadii[idxelem + 1]);
         for (int idxtheta = 0; idxtheta < grid.NumberOfCoLatitudes();
              ++idxtheta) {
            for (int idxphi = 0; idxphi < grid.NumberOfLongitudes(); ++idxphi) {
               tmp_h[idxtheta * grid.NumberOfLongitudes() + idxphi] =
                   mypertprem.RadialMap(radr, mytheta[idxtheta], myphi[idxphi]);
            }
         }
         vec_h[idxelem * npoly + idxpoly] = tmp_h;
      }
   }
   // std::cout << "Hello 1 \n";
   // convert to spherical harmonics
   using veccomp = std::vector<std::complex<double>>;
   using vecvech = std::vector<veccomp>;

   vecvech vec_hlm(matlen, veccomp(_size0, 0.0));
   // std::cout << "intsize: " << intsize
   //           << ", size: " << GSHIndices<NonNegative>(lMax, lMax, 0).size()
   //           << std::endl;
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {

      int idxpolymax = npoly + 1;
      if (idxelem < nelem - 1) {
         idxpolymax = npoly;
      } else {
         idxpolymax = npoly + 1;
      }
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         auto myidx = idxelem * npoly + idxpoly;
         grid.ForwardTransformation(lMax, 0, vec_h[myidx], vec_hlm[myidx]);
      }
   }
   // std::cout << "Hello 2 \n";
   // finding a
   // first we need to find (\nabla h)^{\alpha}
   using MATRIX3 = Eigen::Matrix<std::complex<double>, 3, 3>;
   auto mat_0 = MATRIX3::Zero();
   std::vector<std::vector<MATRIX3>> vec_a(
       nelem * (npoly + 1), std::vector<MATRIX3>(intsize, mat_0));
   // std::cout << "Hello 3 \n";
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {

      int idxpolymax = npoly + 1;
      int idxpolymin = 0;
      if (idxelem == 0) {
         idxpolymin = 1;
      }
      vecvech vec_hp(npoly + 1, veccomp(_sizepm, 0.0)),
          vec_hm(npoly + 1, veccomp(_sizepm, 0.0)),
          vec_h0(npoly + 1, veccomp(_size0, 0.0));
      double inv2 = 2.0 / vec_elemwidth[idxelem];
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {

         for (int idxl = 0; idxl < lMax + 1; ++idxl) {
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               auto idxref = idxm + idxl * (1 + idxl);
               for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                  auto idxouter = idxelem * npoly + idxn;
                  vec_h0[idxpoly][idxref] += vec_hlm[idxouter][idxref] *
                                             _mat_gaussderiv(idxn, idxpoly) *
                                             inv2;
               }
            }
         }
      }
      for (int idxpoly = idxpolymin; idxpoly < idxpolymax; ++idxpoly) {
         auto idxouter = idxelem * npoly + idxpoly;
         for (int idxl = 1; idxl < lMax + 1; ++idxl) {
            // omega_l^0
            double Omegal0 = std::sqrt(static_cast<double>(idxl) *
                                       static_cast<double>(idxl + 1) / 2.0);
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               auto idx0 = idxm + idxl * (1 + idxl);
               auto idxref = idxm + idxl * (1 + idxl) - 1;
               vec_hp[idxpoly][idxref] =
                   vec_hlm[idxouter][idxref] * Omegal0 / vec_allradii[idxouter];
               vec_hm[idxpoly][idxref] = vec_hp[idxpoly][idxref];
            }
         }
      }

      // transform back into spatial
      vecvech vec_nh0(npoly + 1, veccomp(intsize, 0.0)),
          vec_nhp(npoly + 1, veccomp(intsize, 0.0)),
          vec_nhm(npoly + 1, veccomp(intsize, 0.0));
      // std::cout << "Hello 4 \n ";
      // std::cout << "intsize: " << intsize
      //           << ", size: " << GSHIndices<NonNegative>(lMax, lMax,
      //           0).size()
      //           << std::endl;
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         grid.InverseTransformation(lMax, 0, vec_h0[idxpoly], vec_nh0[idxpoly]);
         grid.InverseTransformation(lMax, 1, vec_hp[idxpoly], vec_nhp[idxpoly]);
         grid.InverseTransformation(lMax, -1, vec_hm[idxpoly],
                                    vec_nhm[idxpoly]);
      }
      // std::cout << "Hello 5 \n";

      // find a
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         for (int idxrad = 0; idxrad < intsize; ++idxrad) {
            auto idxuse = idxelem * (npoly + 1) + idxpoly;
            auto nhm = vec_nhm[idxpoly][idxrad];
            auto nhp = vec_nhp[idxpoly][idxrad];
            auto nh0 = vec_nh0[idxpoly][idxrad];
            std::complex<double> hdivr;
            if (idxuse == 0) {
               hdivr = 0.0;
            } else {
               hdivr = vec_h[idxelem * npoly + idxpoly][idxrad] /
                       vec_allradii[idxelem * npoly + idxpoly];
            }
            vec_a[idxuse][idxrad](0, 1) = -nhm;
            vec_a[idxuse][idxrad](1, 0) = -nhm;
            vec_a[idxuse][idxrad](1, 2) = -nhp;
            vec_a[idxuse][idxrad](2, 1) = -nhp;
            vec_a[idxuse][idxrad](0, 2) = -(1.0 + nh0);
            vec_a[idxuse][idxrad](2, 0) = -(1.0 + nh0);
            vec_a[idxuse][idxrad](1, 1) =
                1.0 - nh0 + 2.0 * hdivr +
                ((nh0 - hdivr) * (nh0 - hdivr) - 2.0 * nhp * nhm) / (1.0 + nh0);
         }
      }
   }

   // check
   int ncheck1 = 0;
   int ncheck2 = 20;

   Timer timer3;
   timer3.start();
   // pseudospectral matrix replacement class
   MatrixReplacement<Complex> mymatrix(vec_noderadii, grid, q, vec_a);

   Eigen::VectorXcd myvectrial = mymatrix * mybigsol;
   // for (int idx = 0; idx < 2 * matlen; ++idx) {
   //    std::cout << "Difference: " << std::isnan(myvectrial(idx).real()) << "
   //    "
   //              << myvectrial(idx) << std::endl;
   // }

   // bicgstab
   Eigen::BiCGSTAB<MatrixReplacement<Complex>,
                   Eigen::SphericalGeometryPreconditioner<Complex>>
       solver;
   solver.compute(mymatrix);
   solver.preconditioner().addmatrix(mybigtest.specmat());
   timer3.stop("Time to initialise matrix replacement and precondition");

   solver.setTolerance(std::pow(10.0, -6.0));
   Eigen::VectorXcd testsol = solver.solve(vec_fullforce);
   std::cout << testsol(0) << std::endl;
   // Eigen::VectorXcd testsol = solver.solveWithGuess(vec_fullforce, mybigsol);
   timer3.stop("Time to solve using BiCGSTAB and matrix replacement");
   std::cout << "Number of iterations: " << solver.iterations() << std::endl;
   std::cout << "Size of matrix: " << mymatrix.rows() << std::endl;

   // converting to grid, ie doing inverse transformation:
   auto vec_gridsol2 = SphdecompToGridSol(testsol, grid, nelem, npoly, size);

   //////////////////////////////////////////////////////////////////
   // outputting result at a particular
   std::string pathtofile = "./work/BigSolAsph.out";
   auto file2 = std::ofstream(pathtofile);
   for (int i = 0; i < nelem + 1; ++i) {
      file2 << std::setprecision(16) << vec_noderadii[i] << ";"
            << vec_premsol[0][i];
      for (int idxl1 = 0; idxl1 < lMax + 1; ++idxl1) {
         for (int idxl2 = 0; idxl2 < lMax + 1; ++idxl2) {
            file2 << ";" << vec_gridsol[i * npoly](idxl1, idxl2).real() << ";"
                  << vec_gridsol[i * npoly](idxl1, idxl2).imag() << ";"
                  << vec_gridsol2[i * npoly](idxl1, idxl2).real() << ";"
                  << vec_gridsol2[i * npoly](idxl1, idxl2).imag();
         }
      }
      file2 << std::endl;
      // file2 << std::endl;
   };
}
