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
#include <ranges>
#include <vector>

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
   using MATRIX = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;
   using Scalar = Complex;
   using MRange = All;
   using NRange = All;
   using Grid = GaussLegendreGrid<Real, MRange, NRange>;

   auto lMax = 5;                  // max L
   auto nMax = 1;                  // upper index
   auto grid = Grid(lMax, nMax);   // grid for transformations
   auto size = GSHIndices<All>(lMax, lMax, 0).size();   // dof
   std::cout << "Upper index: " << grid.MinUpperIndex() << std::endl;
   // polynomial order
   int npoly = 6;

   //////////////////////////////////////////////////////////////////

   // background 1D model
   auto myprem = PREM();

   // number of layers and scale factor
   int nlayer = myprem.NumberOfLayers();
   double rad_scale = myprem.OuterRadius();

   // find the nodes
   std::vector<double> vec_noderadii = Radial_Node(myprem, 1000000.0);

   int nelem = vec_noderadii.size() - 1;   // total number of elements
   int matlen = nelem * npoly + 1;         // size of matrix
   for (auto &didx : vec_noderadii) {
      // std::cout << "vec_noderadii: " << didx << std::endl;
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
   //???????
   // what is the definition of vec_denslm
   //???????

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
   //////////////////////////////////////////////////////////////////
   // calculate the gradient of zeta
   // vector across grid at one radius
   auto vec_sol = FFTWpp::vector<Scalar>(grid.NumberOfLongitudes() *
                                         grid.NumberOfCoLatitudes());
   laynum = 0;
   for (int idxelem = 0; idxelem < nelem;
        ++idxelem) {   // finding the layer number

      if (!(vec_noderadii[idxelem] <
            myprem.UpperRadius(laynum) / myprem.OuterRadius())) {
         laynum += 1;
      };
      for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
      };
   };

   // sizes of vectors at one radius
   auto size0 = GSHIndices<All>(lMax, lMax, 0).size();     // 0 component
   auto sizep1 = GSHIndices<All>(lMax, lMax, 1).size();    //+ component
   auto sizem1 = GSHIndices<All>(lMax, lMax, -1).size();   //- component

   // vectors to hold data:
   auto zeta_lm0 = FFTWpp::vector<Complex>(size0);   // scalar field
   auto dzeta_lm0 =
       FFTWpp::vector<Complex>(size0);   // 0 component of derivative
   auto dzeta_lm1 = FFTWpp::vector<Complex>(sizep1);    //+ ""
   auto dzeta_lmm1 = FFTWpp::vector<Complex>(sizem1);   //- ""
   dzeta_lm0 = zeta_lm0;

   // Eigen matrix to hold derivative information for the quadrature
   auto pleg =
       Interpolation::LagrangePolynomial(q.Points().begin(), q.Points().end());
   Eigen::MatrixXcd mat_gaussderiv(npoly + 1, npoly + 1);
   for (int idxi = 0; idxi < npoly + 1; ++idxi) {
      for (int idxj = 0; idxj < npoly + 1; ++idxj) {
         mat_gaussderiv(idxi, idxj) = pleg.Derivative(idxi, q.X(idxj));
      }
   }
   // mat_GD(i,j) is the derivative of the ith Lagrange polynomial at the jth
   // point

   // vectors for (\nabla \zeta):
   // int matlen = nelem * npoly + 1; // size of matrix
   std::vector<Complex> nz0(size0 * matlen, 0.0), nzp1(sizep1 * matlen, 0.0),
       nzm1(sizem1 * matlen, 0.0);
   std::cout << "nz0: " << nz0.size() << " " << nzp1.size() << std::endl;

   /////////////////////////////////////////////////////////////////
   // finding (\nabla \zeta)_{lm}^0 from vector //
   /////////////////////////////////////////////////////////////////

   // looping over elements
   laynum = 0;
   auto myidxfunc = [npoly, lMax](int idxelem, int idxpoly, int idxl,
                                  int idxm) {
      return static_cast<int>(std::pow(idxl, 2) + idxl + idxm +
                              (idxpoly + npoly * idxelem) *
                                  std::pow(lMax + 1, 2));
   };
   auto myidxpm = [npoly, lMax](int idxelem, int idxpoly, int idxl, int idxm) {
      return static_cast<int>(std::pow(idxl, 2) + idxl + idxm - 1 +
                              (idxpoly + npoly * idxelem) *
                                  (std::pow(lMax + 1, 2) - 1));
   };

   // index return
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // finding the layer number
      if (!(vec_noderadii[idxelem] <
            myprem.UpperRadius(laynum) / myprem.OuterRadius())) {
         laynum += 1;
      };

      // looping over quadrature nodes
      int idxpolymax = npoly + 1;
      if (idxelem < nelem - 1) {
         idxpolymax = npoly;
      } else {
         idxpolymax = npoly + 1;
      }
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         // looping over l,m values
         for (int idxl = 0; idxl < lMax + 1; ++idxl) {

            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               // index of vector
               std::size_t mynewidx = myidxfunc(idxelem, idxpoly, idxl, idxm);

               // sum over all derivatives h_n'(r) within the element
               for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                  nz0[mynewidx] +=
                      mybigsol(myidxfunc(idxelem, idxn, idxl, idxm)) *
                      mat_gaussderiv(idxn, idxpoly) * 2.0 /
                      (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]);
               }
            }
         }
      };
   };
   laynum = 0;
   // index return
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // finding the layer number
      if (!(vec_noderadii[idxelem] <
            myprem.UpperRadius(laynum) / myprem.OuterRadius())) {
         laynum += 1;
      };

      // looping over quadrature nodes
      int idxpolymax = npoly + 1;
      if (idxelem < nelem - 1) {
         idxpolymax = npoly;
      } else {
         idxpolymax = npoly + 1;
      }
      int idxlowerbound = 0;
      if (idxelem == 0) {
         idxlowerbound = 1;
      }
      for (int idxpoly = idxlowerbound; idxpoly < idxpolymax; ++idxpoly) {
         // looping over l,m values
         if (lMax > 0) {
            for (int idxl = 1; idxl < lMax + 1; ++idxl) {
               double Omegal0 = std::sqrt(static_cast<double>(idxl) *
                                          static_cast<double>(idxl + 1) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  // index of vector
                  std::size_t mynewidx = myidxpm(idxelem, idxpoly, idxl, idxm);
                  nzp1[mynewidx] +=
                      mybigsol(myidxfunc(idxelem, idxpoly, idxl, idxm)) *
                      Omegal0;
                  nzp1[mynewidx] *=
                      std::pow(GravityFunctions::StandardIntervalMap(
                                   q.X(idxpoly), vec_noderadii[idxelem],
                                   vec_noderadii[idxelem + 1]),
                               -1.0);
                  nzm1[mynewidx] += nzp1[mynewidx];
               }
            }
         }
      };
   };

   /////////////////////////////////////////////////////////////////
   // transforming (\nabla \zeta) to spacial basis (via GSPH) //
   /////////////////////////////////////////////////////////////////
   // do transformations at each radius into the spatial domain
   std::vector<FFTWpp::vector<Scalar>> vec_mzeta, vec_pzeta, vec_zzeta;
   laynum = 0;
   // index return
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // finding the layer number
      if (!(vec_noderadii[idxelem] <
            myprem.UpperRadius(laynum) / myprem.OuterRadius())) {
         laynum += 1;
      };

      // looping over quadrature nodes
      int idxpolymax = npoly + 1;
      if (idxelem < nelem - 1) {
         idxpolymax = npoly;
      } else {
         idxpolymax = npoly + 1;
      }
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         auto wp = std::ranges::subrange(
             nzp1.begin() + sizep1 * (idxelem * npoly + idxpoly),
             nzp1.begin() + sizep1 * (idxelem * npoly + idxpoly + 1));
         auto wm = std::ranges::subrange(
             nzm1.begin() + sizem1 * (idxelem * npoly + idxpoly),
             nzm1.begin() + sizem1 * (idxelem * npoly + idxpoly + 1));
         auto w0 = std::ranges::subrange(
             nz0.begin() + size0 * (idxelem * npoly + idxpoly),
             nz0.begin() + size0 * (idxelem * npoly + idxpoly + 1));
         // if (idxelem == nelem - 1)
         // {
         //     if (idxpoly == npoly)
         //     {
         //         std::cout << "Size 0: " << size0 * (idxelem * npoly +
         //         idxpoly + 1) << " " << nz0.size() << std::endl; std::cout <<
         //         "Sizes: " << sizep1 * (idxelem * npoly + idxpoly + 1) << " "
         //         << nzp1.size() << std::endl; std::cout << "Sizes: " <<
         //         sizem1 * (idxelem * npoly + idxpoly + 1) << " " <<
         //         nzm1.size() << std::endl;
         //     }
         // }

         grid.InverseTransformation(grid.MaxDegree(), 0, w0, vec_sol);
         vec_zzeta.push_back(vec_sol);
         grid.InverseTransformation(grid.MaxDegree(), 1, wp, vec_sol);
         vec_pzeta.push_back(vec_sol);
         grid.InverseTransformation(grid.MaxDegree(), -1, wm, vec_sol);
         vec_mzeta.push_back(vec_sol);
      }
   }
   std::cout << "Size of 0: " << vec_zzeta.size() << ", 1: " << vec_pzeta.size()
             << ", -1: " << vec_mzeta.size() << std::endl;
   std::cout << lMax << " " << grid.MaxDegree() << "\n";
   /////////////////////////////////////////////////////////////////
   // finding q = a\nabla \zeta, need to multiply at each spatial point//
   /////////////////////////////////////////////////////////////////
   // now we need to multiply by the matrix a
   std::vector<std::vector<MATRIX>> vec_a(nelem * npoly + 1);
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {

      int idxpolymax = npoly + 1;
      if (idxelem < nelem - 1) {
         idxpolymax = npoly;
      } else {
         idxpolymax = npoly + 1;
      }
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         std::vector<MATRIX> single_radius_a(grid.NumberOfLongitudes() *
                                             grid.NumberOfCoLatitudes());
         for (int idxradius = 0; idxradius < grid.NumberOfLongitudes() *
                                                 grid.NumberOfCoLatitudes();
              ++idxradius) {
            MATRIX mat_a(3, 3);
            mat_a << 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0;
            single_radius_a[idxradius] = mat_a;
         }
         vec_a[idxelem * npoly + idxpoly] = single_radius_a;
      }
   }

   // now multiplying the matrix through at each point
   std::vector<FFTWpp::vector<Scalar>> vec_mq(matlen), vec_pq(matlen),
       vec_zq(matlen);
   // std::vector<Complex> vec_mq(sizem1 * matlen, 0.0), vec_pq(sizep1 * matlen,
   // 0.0), vec_zq(size0 * matlen, 0.0);
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      std::vector<MATRIX> single_radius_a(grid.NumberOfLongitudes() *
                                          grid.NumberOfCoLatitudes());

      int idxpolymax = npoly + 1;
      if (idxelem < nelem - 1) {
         idxpolymax = npoly;
      } else {
         idxpolymax = npoly + 1;
      }
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         FFTWpp::vector<Scalar> vec_mqadd(
             grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes(), 0.0),
             vec_pqadd(grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes(),
                       0.0),
             vec_zqadd(grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes(),
                       0.0);
         for (int idxr = 0;
              idxr < grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();
              ++idxr) {
            vec_mqadd[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](0, 0) *
                               vec_mzeta[idxelem * npoly + idxpoly][idxr];
            vec_mqadd[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](0, 1) *
                               vec_zzeta[idxelem * npoly + idxpoly][idxr];
            vec_mqadd[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](0, 2) *
                               vec_pzeta[idxelem * npoly + idxpoly][idxr];
            vec_zqadd[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](1, 0) *
                               vec_mzeta[idxelem * npoly + idxpoly][idxr];
            vec_zqadd[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](1, 1) *
                               vec_zzeta[idxelem * npoly + idxpoly][idxr];
            vec_zqadd[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](1, 2) *
                               vec_pzeta[idxelem * npoly + idxpoly][idxr];
            vec_pqadd[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](2, 0) *
                               vec_mzeta[idxelem * npoly + idxpoly][idxr];
            vec_pqadd[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](2, 1) *
                               vec_zzeta[idxelem * npoly + idxpoly][idxr];
            vec_pqadd[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](2, 2) *
                               vec_pzeta[idxelem * npoly + idxpoly][idxr];
         }
         vec_mq[idxelem * npoly + idxpoly] = vec_mqadd;
         vec_zq[idxelem * npoly + idxpoly] = vec_zqadd;
         vec_pq[idxelem * npoly + idxpoly] = vec_pqadd;
      }
      // vec_a[idxelem] = single_radius_a;
   }
   // std::cout << "Size of bigsol: " << mybigsol.size() << " " << size0 *
   // matlen << "\n";
   //  grid.InverseTransformation(grid.MaxDegree(), 0, dzeta_lm0, vec_sol);

   // std::copy(zeta_lm0.begin() + 1, zeta_lm0.end(), dzeta_lm1.begin());

   // std::ranges::subrange w(dzeta_lm0.begin() + 1, dzeta_lm0.end());
   /////////////////////////////////////////////////////////////////
   // finding q_{lm}^{\alpha} //
   /////////////////////////////////////////////////////////////////
   std::vector<Complex> q0(size0 * matlen, 0.0), qp1(sizep1 * matlen, 0.0),
       qm1(sizem1 * matlen, 0.0);
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      int idxpolymax = npoly + 1;
      if (idxelem < nelem - 1) {
         idxpolymax = npoly;
      } else {
         idxpolymax = npoly + 1;
      }
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         // transformation for q^{-1}
         auto wm1 =
             std::ranges::subrange(vec_mq[idxelem * npoly + idxpoly].begin(),
                                   vec_mq[idxelem * npoly + idxpoly].end());
         auto qlmm1out = std::ranges::subrange(
             qm1.begin() + (idxelem * npoly + idxpoly) * sizem1,
             qm1.begin() + (idxelem * npoly + idxpoly + 1) * sizem1);
         grid.ForwardTransformation(grid.MaxDegree(), -1, wm1, qlmm1out);

         // transformation for q^0
         auto w0 =
             std::ranges::subrange(vec_zq[idxelem * npoly + idxpoly].begin(),
                                   vec_zq[idxelem * npoly + idxpoly].end());
         auto qlmout = std::ranges::subrange(
             q0.begin() + (idxelem * npoly + idxpoly) * size0,
             q0.begin() + (idxelem * npoly + idxpoly + 1) * size0);
         grid.ForwardTransformation(grid.MaxDegree(), 0, w0, qlmout);

         // transformation for q^{+1}
         auto wp1 =
             std::ranges::subrange(vec_pq[idxelem * npoly + idxpoly].begin(),
                                   vec_pq[idxelem * npoly + idxpoly].end());
         auto qlmp1out = std::ranges::subrange(
             qp1.begin() + (idxelem * npoly + idxpoly) * sizep1,
             qp1.begin() + (idxelem * npoly + idxpoly + 1) * sizep1);
         grid.ForwardTransformation(grid.MaxDegree(), +1, wp1, qlmp1out);
      }
   }
   // the values of q_{lm}^{-1} at the first point
   // for (int i = 0; i < sizep1; ++i)
   // {
   //     std::cout << qp1[i] << qm1[i] << q0[i] << std::endl;
   // }
   /////////////////////////////////////////////////////////////////
   // integrating over 0 to b//
   /////////////////////////////////////////////////////////////////
   std::vector<Complex> vec_output(matlen * size0);
   laynum = 0;
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // finding the layer number
      if (!(vec_noderadii[idxelem] <
            myprem.UpperRadius(laynum) / myprem.OuterRadius())) {
         laynum += 1;
      };

      // looping over quadrature nodes
      int idxpolymax = npoly + 1;
      if (idxelem < nelem - 1) {
         // idxpolymax = npoly;
      } else {
         idxpolymax = npoly + 1;
      }
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         // looping over l,m values
         for (int idxl = 0; idxl < lMax + 1; ++idxl) {

            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               // index of vector
               std::size_t mynewidx = myidxfunc(idxelem, idxpoly, idxl, idxm);
               // std::cout << mynewidx << std::endl;

               // sum over all derivatives h_n'(r) within the element
               for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                  vec_output[mynewidx] +=
                      q0[myidxfunc(idxelem, idxn, idxl, idxm)] *
                      std::pow(GravityFunctions::StandardIntervalMap(
                                   q.X(idxn), vec_noderadii[idxelem],
                                   vec_noderadii[idxelem + 1]),
                               2.0) *
                      mat_gaussderiv(idxpoly, idxn) * q.W(idxn);
               }
               if (idxl > 0) {
                  double omegal0 =
                      std::sqrt(static_cast<double>(idxl) *
                                (static_cast<double>(idxl + 1)) * 0.5);
                  vec_output[mynewidx] +=
                      GravityFunctions::StandardIntervalMap(
                          q.X(idxpoly), vec_noderadii[idxelem],
                          vec_noderadii[idxelem + 1]) *
                      omegal0 * q.W(idxpoly) *
                      (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) *
                      0.5 *
                      (qm1[myidxpm(idxelem, idxpoly, idxl, idxm)] +
                       qp1[myidxpm(idxelem, idxpoly, idxl, idxm)]);
               }
               if (idxelem == nelem - 1 && idxpoly == npoly) {
                  vec_output[mynewidx] += (static_cast<double>(idxl) + 1.0) *
                                          vec_noderadii.back() *
                                          mybigsol(mynewidx);
               }
            }
         }
      };
   };

   /////////////////////////////////////////////////////////////
   // testing alternate method of doing it via loop over all nodes
   //  index return
   std::vector<Complex> vec_output2(matlen * size0, 0.0);
   FFTWpp::vector<Scalar> vec_smallmzeta(grid.NumberOfLongitudes() *
                                         grid.NumberOfCoLatitudes()),
       vec_smallpzeta(grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes()),
       vec_smallzzeta(grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes());

   // auto vec_sol = FFTWpp::vector<Scalar>(grid.NumberOfLongitudes() *
   // grid.NumberOfCoLatitudes());

   auto idxsmall = [npoly, lMax](int idxl, int idxm) {
      return static_cast<int>(std::pow(idxl, 2) + idxl + idxm);
   };
   auto idxpmsmall = [npoly, lMax](int idxl, int idxm) {
      return static_cast<int>(std::pow(idxl, 2) + idxl + idxm - 1);
   };
   int tcount = 0;
   laynum = 0;
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // finding the layer number
      if (!(vec_noderadii[idxelem] <
            myprem.UpperRadius(laynum) / myprem.OuterRadius())) {
         laynum += 1;
      };

      // looping over quadrature nodes
      int idxpolymax = npoly + 1;
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         std::vector<Complex> vec_smallnz0(size0, 0.0),
             vec_smallnzp1(sizep1, 0.0), vec_smallnzm1(sizem1, 0.0);
         std::vector<Complex> vec_smallq0(size0, 0.0),
             vec_smallqp1(sizep1, 0.0), vec_smallqm1(sizem1, 0.0);
         FFTWpp::vector<Scalar> vec_smallmq(
             grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes(), 0.0),
             vec_smallpq(grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes(),
                         0.0),
             vec_smallzq(grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes(),
                         0.0);
         ///////////////////////////////////////
         // find nabla zeta
         // looping over l,m values
         for (int idxl = 0; idxl < lMax + 1; ++idxl) {
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               // index of vector
               std::size_t mynewidx = myidxfunc(idxelem, idxpoly, idxl, idxm);
               std::size_t mynewidx2 = idxsmall(idxl, idxm);

               // sum over all derivatives h_n'(r) within the element
               for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                  vec_smallnz0[mynewidx2] +=
                      mybigsol(myidxfunc(idxelem, idxn, idxl, idxm)) *
                      mat_gaussderiv(idxn, idxpoly) * 2.0 /
                      (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]);
               }
            }
         }
         if (tcount > 0) {
            for (int idxl = 1; idxl < lMax + 1; ++idxl) {
               double Omegal0 = std::sqrt(static_cast<double>(idxl) *
                                          static_cast<double>(idxl + 1) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  // index of vector
                  std::size_t mynewidx = idxpmsmall(idxl, idxm);
                  if (idxelem == 0 && idxpoly == 0) {
                     // std::cout << mynewidx << "\n";
                  }
                  vec_smallnzp1[mynewidx] +=
                      mybigsol(myidxfunc(idxelem, idxpoly, idxl, idxm)) *
                      Omegal0;
                  vec_smallnzp1[mynewidx] *=
                      std::pow(GravityFunctions::StandardIntervalMap(
                                   q.X(idxpoly), vec_noderadii[idxelem],
                                   vec_noderadii[idxelem + 1]),
                               -1.0);
                  vec_smallnzm1[mynewidx] += vec_smallnzp1[mynewidx];
               }
            }
         }
         ///////////////////////////////////////////////////////////////////
         ///////////////////////////////////////////////////////////////////
         ///////////////////////////////////////////////////////////////////
         // debugging
         // if (idxelem == 0 && idxpoly == 0)
         // {
         //     for (int idxl = 0; idxl < lMax + 1; ++idxl)
         //     {
         //         for (int idxm = -idxl; idxm < idxl + 1; ++idxm)
         //         {
         //             std::cout << "Old: " << nz0[myidxfunc(idxelem, idxpoly,
         //             idxl, idxm)] << ", new: " << vec_smallnz0[idxsmall(idxl,
         //             idxm)] << "\n"; if (idxl > 0)
         //             {
         //                 std::cout << "Old p1: " << nzp1[myidxpm(idxelem,
         //                 idxpoly, idxl, idxm)] << ", new: " <<
         //                 vec_smallnzp1[idxpmsmall(idxl, idxm)] << "\n";
         //             }
         //         }
         //     }
         // }
         // for (int idxl = 0; idxl < lMax + 1; ++idxl)
         // {
         //     for (int idxm = -idxl; idxm < idxl + 1; ++idxm)
         //     {

         //         auto norm1 = std::abs(nz0[myidxfunc(idxelem, idxpoly, idxl,
         //         idxm)]); auto normdiff = std::abs(nz0[myidxfunc(idxelem,
         //         idxpoly, idxl, idxm)] - vec_smallnz0[idxsmall(idxl, idxm)]);
         //         if (norm1 > std::numeric_limits<double>::epsilon())
         //         {
         //             if (normdiff / norm1 >
         //             std::numeric_limits<double>::epsilon())
         //             {
         //                 std::cout << std::setprecision(16) << "Element: " <<
         //                 idxelem << ", node: " << idxpoly << ", l: " << idxl
         //                 << ", m: " << idxm << " Old: " <<
         //                 nz0[myidxfunc(idxelem, idxpoly, idxl, idxm)] << ",
         //                 new: " << vec_smallnz0[idxsmall(idxl, idxm)] <<
         //                 "\n";
         //             }
         //         }
         //     }
         // }

         ///////////////////////////////////////////////////////////////////
         ///////////////////////////////////////////////////////////////////
         ///////////////////////////////////////////////////////////////////

         auto wp =
             std::ranges::subrange(vec_smallnzp1.begin(), vec_smallnzp1.end());
         auto wm =
             std::ranges::subrange(vec_smallnzm1.begin(), vec_smallnzm1.end());
         auto w0 =
             std::ranges::subrange(vec_smallnz0.begin(), vec_smallnz0.end());

         grid.InverseTransformation(grid.MaxDegree(), 0, w0, vec_smallzzeta);
         grid.InverseTransformation(grid.MaxDegree(), 1, wp, vec_smallpzeta);
         grid.InverseTransformation(grid.MaxDegree(), -1, wm, vec_smallmzeta);
         ///////////////////////////////////////////////////////////////////
         ///////////////////////////////////////////////////////////////////
         ///////////////////////////////////////////////////////////////////
         // debugging
         // if (idxelem == 0 && idxpoly == 0)
         // {
         //     for (int idxl = 0; idxl < lMax + 1; ++idxl)
         //     {
         //         for (int idxm = -idxl; idxm < idxl + 1; ++idxm)
         //         {
         //             std::cout << "Old: " << nz0[myidxfunc(idxelem, idxpoly,
         //             idxl, idxm)] << ", new: " << vec_smallnz0[idxsmall(idxl,
         //             idxm)] << "\n"; if (idxl > 0)
         //             {
         //                 std::cout << "Old p1: " << nzp1[myidxpm(idxelem,
         //                 idxpoly, idxl, idxm)] << ", new: " <<
         //                 vec_smallnzp1[idxpmsmall(idxl, idxm)] << "\n";
         //             }
         //         }
         //     }
         // }
         for (int idxl = 0; idxl < lMax + 1; ++idxl) {
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               for (int idxr = 0; idxr < grid.NumberOfLongitudes() *
                                             grid.NumberOfCoLatitudes();
                    ++idxr) {
                  // std::cout << "Element: " << idxelem << ", npoly: " <<
                  // idxpoly << ", l: " << idxl << ", m: " << idxm << ", r: " <<
                  // idxr << ", new idx: " << myidxfunc(idxelem, idxpoly, idxl,
                  // idxm) << std::endl;
                  auto norm1 =
                      std::abs(vec_pzeta[idxelem * npoly + idxpoly][idxr]);
                  // std::cout << "l: " << idxl << ", m: " << idxm << ", r: " <<
                  // idxr << ", norm:" << norm1 << std::endl;
                  auto normdiff =
                      std::abs(vec_smallpzeta[idxr] -
                               vec_pzeta[idxelem * npoly + idxpoly][idxr]);
                  if (norm1 > std::numeric_limits<double>::epsilon()) {
                     if (normdiff / norm1 >
                             std::numeric_limits<double>::epsilon() &&
                         idxpoly != 7) {
                        // std::cout << std::setprecision(16) << "Element: " <<
                        // idxelem << ", node: " << idxpoly << ", l: " << idxl
                        // << ", m: " << idxm << " Old: " << vec_zzeta[idxelem *
                        // npoly + idxpoly][idxr] << ", new: " <<
                        // vec_smallzzeta[idxr] << "\n";
                     }
                  }
               }
            }
         }

         ///////////////////////////////////////////////////////////////////
         ///////////////////////////////////////////////////////////////////
         ///////////////////////////////////////////////////////////////////

         for (int idxr = 0;
              idxr < grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();
              ++idxr) {
            vec_smallmq[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](0, 0) *
                                 vec_smallmzeta[idxr];
            vec_smallmq[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](0, 1) *
                                 vec_smallzzeta[idxr];
            vec_smallmq[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](0, 2) *
                                 vec_smallpzeta[idxr];
            vec_smallzq[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](1, 0) *
                                 vec_smallmzeta[idxr];
            vec_smallzq[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](1, 1) *
                                 vec_smallzzeta[idxr];
            vec_smallzq[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](1, 2) *
                                 vec_smallpzeta[idxr];
            vec_smallpq[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](2, 0) *
                                 vec_smallmzeta[idxr];
            vec_smallpq[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](2, 1) *
                                 vec_smallzzeta[idxr];
            vec_smallpq[idxr] += vec_a[idxelem * npoly + idxpoly][idxr](2, 2) *
                                 vec_smallpzeta[idxr];
         }

         // transformation for q^{-1}
         auto ftwm1 =
             std::ranges::subrange(vec_smallmq.begin(), vec_smallmq.end());
         auto qlmm1out =
             std::ranges::subrange(vec_smallqm1.begin(), vec_smallqm1.end());
         grid.ForwardTransformation(grid.MaxDegree(), -1, ftwm1, qlmm1out);

         // transformation for q^0
         auto ftw0 =
             std::ranges::subrange(vec_smallzq.begin(), vec_smallzq.end());
         auto qlmout =
             std::ranges::subrange(vec_smallq0.begin(), vec_smallq0.end());
         grid.ForwardTransformation(grid.MaxDegree(), 0, ftw0, qlmout);

         // transformation for q^{+1}
         auto ftwp1 =
             std::ranges::subrange(vec_smallpq.begin(), vec_smallpq.end());
         auto qlmp1out =
             std::ranges::subrange(vec_smallqp1.begin(), vec_smallqp1.end());
         grid.ForwardTransformation(grid.MaxDegree(), +1, ftwp1, qlmp1out);

         ///////////////////////////////////////////////////////////////////
         ///////////////////////////////////////////////////////////////////
         ///////////////////////////////////////////////////////////////////
         // debugging

         // for (int idxl = 0; idxl < lMax + 1; ++idxl)
         // {
         //     for (int idxm = -idxl; idxm < idxl + 1; ++idxm)
         //     {

         //         // std::cout << "Element: " << idxelem << ", npoly: " <<
         //         idxpoly << ", l: " << idxl << ", m: " << idxm << ", r: " <<
         //         idxr << ", new idx: " << myidxfunc(idxelem, idxpoly, idxl,
         //         idxm) << std::endl; auto norm1 =
         //         std::abs(q0[myidxfunc(idxelem, idxpoly, idxl, idxm)]); auto
         //         norm2 = std::abs(qp1[myidxpm(idxelem, idxpoly, idxl,
         //         idxm)]); auto norm3 = std::abs(qm1[myidxpm(idxelem, idxpoly,
         //         idxl, idxm)]);
         //         // std::cout << "l: " << idxl << ", m: " << idxm << ", r: "
         //         << idxr << ", norm:" << norm1 << std::endl; auto normdiff =
         //         std::abs(vec_smallq0[idxsmall(idxl, idxm)] -
         //         q0[myidxfunc(idxelem, idxpoly, idxl, idxm)]); auto normdiff2
         //         = std::abs(vec_smallqp1[idxpmsmall(idxl, idxm)] -
         //         qp1[myidxpm(idxelem, idxpoly, idxl, idxm)]); auto normdiff3
         //         = std::abs(vec_smallqm1[idxpmsmall(idxl, idxm)] -
         //         qm1[myidxpm(idxelem, idxpoly, idxl, idxm)]); if (norm1 >
         //         std::numeric_limits<double>::epsilon())
         //         {
         //             if (normdiff / norm1 >
         //             std::numeric_limits<double>::epsilon() && idxpoly != 7)
         //             {
         //                 std::cout << std::setprecision(16) << "Element: " <<
         //                 idxelem << ", node: " << idxpoly << ", l: " << idxl
         //                 << ", m: " << idxm << " Old: " <<
         //                 q0[myidxfunc(idxelem, idxpoly, idxl, idxm)] << ",
         //                 new: " << vec_smallq0[idxsmall(idxl, idxm)] << "\n";
         //             }
         //         }
         //         if (norm2 > std::numeric_limits<double>::epsilon())
         //         {
         //             if (normdiff2 / norm2 >
         //             std::numeric_limits<double>::epsilon() && idxpoly != 7)
         //             {
         //                 std::cout << std::setprecision(16) << "Element: " <<
         //                 idxelem << ", node: " << idxpoly << ", l: " << idxl
         //                 << ", m: " << idxm << " Old: " <<
         //                 qp1[myidxpm(idxelem, idxpoly, idxl, idxm)] << ",
         //                 new: " << vec_smallqp1[idxpmsmall(idxl, idxm)] <<
         //                 "\n";
         //             }
         //         }
         //         if (norm3 > std::numeric_limits<double>::epsilon())
         //         {
         //             if (normdiff3 / norm3 >
         //             std::numeric_limits<double>::epsilon())
         //             {
         //                 // std::cout << std::setprecision(16) << "Element: "
         //                 << idxelem << ", node: " << idxpoly << ", l: " <<
         //                 idxl << ", m: " << idxm << " Old: " <<
         //                 qm1[myidxpm(idxelem, idxpoly, idxl, idxm)] << ",
         //                 new: " << vec_smallqm1[idxpmsmall(idxl, idxm)] <<
         //                 "\n"; std::cout << std::setprecision(16) <<
         //                 "Element: " << idxelem << ", node: " << idxpoly <<
         //                 ", l: " << idxl << ", m: " << idxm << " Ratio: " <<
         //                 normdiff3 / norm3 << "\n";
         //             }
         //         }
         //     }
         // }

         ///////////////////////////////////////////////////////////////////
         ///////////////////////////////////////////////////////////////////
         ///////////////////////////////////////////////////////////////////
         // looping over l,m values
         for (int idxl = 0; idxl < lMax + 1; ++idxl) {

            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               // index of vector
               std::size_t mynewidx = myidxfunc(idxelem, idxpoly, idxl, idxm);

               if (idxpoly == 0 || idxpoly == npoly) {
                  // std::cout << mynewidx << std::endl;
               }
               // sum over all derivatives h_n'(r) within the element
               for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                  std::size_t mynewidx2 = myidxfunc(idxelem, idxn, idxl, idxm);
                  vec_output2[mynewidx2] -=
                      vec_smallq0[idxsmall(idxl, idxm)] *
                      std::pow(GravityFunctions::StandardIntervalMap(
                                   q.X(idxpoly), vec_noderadii[idxelem],
                                   vec_noderadii[idxelem + 1]),
                               2.0) *
                      mat_gaussderiv(idxn, idxpoly) * q.W(idxpoly);
               }
               if (idxl > 0) {
                  double omegal0 =
                      std::sqrt(static_cast<double>(idxl) *
                                (static_cast<double>(idxl + 1)) * 0.5);
                  vec_output2[mynewidx] -=
                      GravityFunctions::StandardIntervalMap(
                          q.X(idxpoly), vec_noderadii[idxelem],
                          vec_noderadii[idxelem + 1]) *
                      omegal0 * q.W(idxpoly) *
                      (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) *
                      0.5 *
                      (vec_smallqm1[idxpmsmall(idxl, idxm)] +
                       vec_smallqp1[idxpmsmall(idxl, idxm)]);
               }
               if (idxelem == nelem - 1 && idxpoly == npoly) {
                  vec_output2[mynewidx] -= (static_cast<double>(idxl) + 1.0) *
                                           vec_noderadii.back() *
                                           mybigsol(mynewidx);
               }
            }
         }
         ++tcount;
      };
   };

   //////////////////////////////////////////////////////////////////
   // testing output

   // matrix replacement

   MatrixReplacement<Complex> mymatrix(vec_noderadii, grid, q, vec_a);
   // MatrixReplacement<Complex> mymatrix;
   //   int n = 10;
   // Eigen::SparseMatrix<Complex> S(mybigsol.size(), mybigsol.size());
   // mymatrix.attachMyMatrix(S);
   // dense
   timer1.start();
   Eigen::MatrixXcd densemat = Eigen::MatrixXcd(mybigtest.specmat());
   // std::cout << densemat.block(0,0,5,5)<<std::endl;
   Eigen::LDLT<Eigen::MatrixXcd, Eigen::Lower> chol_solver;
   chol_solver.compute(densemat);
   Eigen::MatrixXcd mysol2 = chol_solver.solve(vec_fullforce);
   Eigen::VectorXcd test5 = densemat.selfadjointView<Eigen::Lower>() * mysol2;
   // std::cout << mybigtest.specmat().coeff(0,0)<<std::endl;
   // std::cout <<  test5(0) << std::endl;
   timer1.stop("The time for dense calc is");

   // Eigen::VectorXcd test3 = mymatrix * mybigsol;
   Eigen::VectorXcd test3 = mymatrix * mybigsol;
   Eigen::VectorXcd test4;
   test4 = mybigtest.sphmatmult(mybigsol);
   std::cout << "Testing: " << test3(0) << std::endl;
   // testing preconditioner
   Timer timer3;
   timer3.start();
   Eigen::SphericalGeometryPreconditioner<Complex> myprecond;
   myprecond.compute(mybigtest.specmat());
   Eigen::VectorXcd testsol = myprecond.solve(vec_fullforce);

   timer3.stop("Time for preconditioner");

   timer3.start();
   { Eigen::VectorXcd test6 = mybigtest.specmat() * testsol; }
   timer3.stop("Time for matrix multiplication");

   timer3.start();
   Eigen::VectorXcd test6 = mymatrix * testsol;
   timer3.stop("Time for pseudospectral matrix multiplication");
   for (int i = 0; i < nelem; ++i) {
      //    std::cout << std::setprecision(4) << "Force: " << vec_fullforce(i *
      //    npoly * size0) << ", first: " << std::abs(vec_fullforce(i * npoly *
      //    size0) + vec_output[i * npoly * size0]) << ", second: " <<
      //    std::abs(vec_fullforce(i * npoly * size0) - vec_output2[i * npoly *
      //    size0]) << ", third: " << std::abs(vec_fullforce(i * npoly * size0)
      //    - test3(i * npoly * size0)) << std::endl;
      // std::cout << std::setprecision(4) << "Matrix: " << std::abs(test4(i *
      // npoly * size0) - vec_fullforce[i * npoly * size0]) << ", first: " <<
      // std::abs(test4(i * npoly * size0) + vec_output[i * npoly * size0]) <<
      // ", second: " << std::abs(test4(i * npoly * size0) - vec_output2[i *
      // npoly * size0]) << ", third: " << std::abs(test4(i * npoly * size0) -
      // test3(i * npoly * size0)) << std::endl; std::cout <<
      // std::setprecision(4) << "Exact: " << vec_fullforce(i * npoly * size0)
      // << ", matrix error: " << std::abs(test4(i * npoly * size0) -
      // vec_fullforce(i * npoly * size0)) << ", pseudospectral: " <<
      // std::abs(vec_fullforce(i * npoly * size0) - test3(i * npoly * size0))
      // << std::endl; std::cout << std::setprecision(4) << "Exact: " <<
      // vec_fullforce(i * npoly * size0) << ", matrix: " << test4(i * npoly *
      // size0) << ", pseudospectral: " << test3(i * npoly * size0) <<
      // std::endl; std::cout << vec_fullforce.cwiseAbs().maxCoeff() <<
      // std::endl;
      std::cout << std::setprecision(4) << "Difference sparse: "
                << std::abs((vec_fullforce(i * npoly * size0) -
                             test3(i * npoly * size0)) /
                            vec_fullforce.cwiseAbs().maxCoeff())
                << ", difference dense: "
                << std::abs((vec_fullforce(i * npoly * size0) -
                             test5(i * npoly * size0)) /
                            vec_fullforce.cwiseAbs().maxCoeff())
                << ", difference precon: "
                << std::abs((vec_fullforce(i * npoly * size0) -
                             test6(i * npoly * size0)) /
                            vec_fullforce.cwiseAbs().maxCoeff())
                << std::endl;
   }
   // timer1.stop("The time is");
   // test matrix replacement
   // MatrixReplacement mymat;
   // mymat.addData(vec_noderadii, grid, q);
   // mymat.addData(vec_noderadii, grid, q);
   std::cout << nelem << " " << npoly << " " << matlen << "\n";

   // other test
   // Eigen::SparseMatrix<Complex> sm1(3, 3);
   // {
   //     using T = Eigen::Triplet<Complex>;
   //     std::vector<T> tripletlist;

   //     tripletlist.reserve(3);
   //     tripletlist.push_back(T(0, 0, 0.5));
   //     tripletlist.push_back(T(1, 0, 0.25));
   //     tripletlist.push_back(T(2, 1, 0.9));
   //     sm1.setFromTriplets(tripletlist.begin(), tripletlist.end()); // set it
   //     sm1.makeCompressed();
   //     std::cout << sm1 << std::endl;
   // }
   // Eigen::VectorXcd myvec1, myvec2, dv3;
   // myvec1 = Eigen::VectorXcd::Zero(3);
   // myvec1(0) = 1.0;
   // myvec1(1) = 1.0;
   // myvec2 = sm1 * myvec1;
   // dv3 = sm1.selfadjointView<Eigen::Lower>() * myvec1;
   // std::cout << "Vector 1: " << myvec1 << "\n Vector 2: " << myvec2 << "\n
   // Vector 3: " << dv3 << "\n";

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
