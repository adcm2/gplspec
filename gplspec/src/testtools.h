#ifndef TEST_TOOLS_H
#define TEST_TOOLS_H

#include <Eigen/Core>
#include <Eigen/Dense>
// #include <Eigen/CholmodSupport>
#include <Eigen/Sparse>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <PlanetaryModel/All>
// #include <TomographyModels/All>
#include "SphericalGeometryPreconditioner.h"
#include "Shared_Utilities.h"
#include "Model_Input.h"
#include <FFTWpp/Ranges>
#include <GSHTrans/All>
#include <TomographyModels/All>
#include <algorithm>
#include <complex>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>

#include <chrono>
#include <random>
#include <vector>

namespace TestTools {

enum class CoefficientOrdering {
   RadialClumped,
   SphericalDecomp
};

using MATRIX3cd = Eigen::Matrix<std::complex<double>, 3, 3>;
using EARTHMATRIX3 = std::vector<std::vector<MATRIX3cd>>;
using EARTHVEC = std::vector<std::vector<std::vector<double>>>;
using RADIUSVEC = std::vector<std::vector<std::complex<double>>>;
using EARTHVECX = std::vector<RADIUSVEC>;

template <typename FLOAT, template <typename, typename> class sphericalmodel,
          template <typename, typename> class TOMPERTMODEL,
          template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange>
Eigen::VectorXcd
FindBoundaryPerturbationForce(
    const sphericalmodel<FLOAT, int> &backgroundmodel,
    const TOMPERTMODEL<FLOAT, int> &pertmodel,
    GRID<FLOAT, OrderIndexRange, IndexRange> &grid,
    const CoefficientOrdering &whichsetup,
    const GaussQuad::Quadrature1D<FLOAT> &q,
    const std::vector<FLOAT> &vec_noderadii,
    const std::vector<FFTWpp::vector<std::complex<FLOAT>>> &vec_denslm,
    const int &lMax,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
        &mat_gaussderiv) {
   using Real = FLOAT;
   using Complex = std::complex<Real>;
   std::size_t nelem = vec_noderadii.size() - 1;
   int npoly = q.N() - 1;
   int matlen = nelem * npoly + 1;   // size of matrix

   auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).Size();
   auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).Size();
   auto rphys = [&vec_noderadii](int idxelem, double x) {
      return GPLSpec::detail::StandardIntervalMap(
          x, vec_noderadii[idxelem],
          vec_noderadii[idxelem + 1]);
   };

   //    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _mat_gaussderiv;
   // finding \rho \mathbf{s} in canonical components spatially
   // firstly we need to have the background model density and the perturbation
   // from the perturbing model:
   // go through all layers, find what s is, as we currently just have radial,
   // that should be fairly simple
   int totnum = grid.NumberOfCoLatitudes() * grid.NumberOfLongitudes();

   int laynum = 0;

   auto myidxfunc = [npoly, lMax](int idxelem, int idxpoly, int idxl,
                                  int idxm) {
      return static_cast<std::size_t>(std::pow(idxl, 2) + idxl + idxm +
                                      (idxpoly + npoly * idxelem) *
                                          std::pow(lMax + 1, 2));
   };
   using VECTOR = Eigen::Vector<Complex, Eigen::Dynamic>;
   VECTOR vec_output2 = VECTOR::Zero(matlen * _size0);
   std::vector<double> gaussquadweights, gaussquadpoints;
   auto elemwidth = [&vec_noderadii](int idxelem) {
      return vec_noderadii[idxelem + 1] - vec_noderadii[idxelem];
   };
   for (int idx = 0; idx < q.Points().size(); ++idx) {
      gaussquadpoints.push_back(q.X(idx));
      gaussquadweights.push_back(q.W(idx));
   }
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      double inv2 = 2.0 / (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]);
      if (laynum < backgroundmodel.NumberOfLayers()) {
         if (!(vec_noderadii[idxelem] < backgroundmodel.UpperRadius(laynum))) {
            laynum += 1;
            std::cout << "laynum: " << laynum << ", idxelem: " << idxelem
                      << std::endl;
         };
      }

      // looping over quadrature nodes
      int idxpolymax = npoly + 1;
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {

         auto rval = GPLSpec::detail::StandardIntervalMap(
             gaussquadpoints[idxpoly], vec_noderadii[idxelem],
             vec_noderadii[idxelem + 1]);
         auto invr = 1 / rval;
         auto rval2 = rval * rval;

         /////////////////////////////////////////////////////////////

         // finding spatial q0, qp, qm
         std::vector<Complex> spatial_q0(totnum, 0.0), spatial_qm(totnum, 0.0),
             spatial_qp(totnum, 0.0);

         auto raddensity = backgroundmodel.Density(laynum)(rval);
         for (int idxt = 0; idxt < grid.NumberOfCoLatitudes(); ++idxt) {
            for (int idxp = 0; idxp < grid.NumberOfLongitudes(); ++idxp) {
               auto idxuse = idxt * grid.NumberOfLongitudes() + idxp;
               spatial_q0[idxuse] =
                   pertmodel.RadialMap(rval, grid.CoLatitudes()[idxt],
                                       grid.Longitudes()[idxp]) *
                   raddensity;
            }
         }

         /////////////////////////////////////////////////////////////
         // step 4: transforming s from spatial to spherical harmonic
         // basis
         //  transformation for q^{-1}
         std::vector<Complex> gsph_q0(_size0, 0.0), gsph_qp1(_sizepm, 0.0),
             gsph_qm1(_sizepm, 0.0);   // declaration of GSPH vectors
         {
            // transformation for q^{-1}

            grid.ForwardTransformation(lMax, -1, spatial_qm, gsph_qm1);
            grid.ForwardTransformation(lMax, 0, spatial_q0, gsph_q0);
            grid.ForwardTransformation(lMax, +1, spatial_qp, gsph_qp1);
         }

         // step 5: evaluating radial integrals to give Ax

         // 0 component term
         {
            // multiplication factor and index
            std::size_t idx2 = myidxfunc(idxelem, 0, 0, 0);
            auto mult1 = gaussquadweights[idxpoly] * rval2;

            // loop over radii
            for (int idxn = 0; idxn < npoly + 1; ++idxn) {

               // index and multiplication factor
               std::size_t idx3 = 0;
               auto multfact = mat_gaussderiv(idxn, idxpoly) * mult1;

               // loop over l and m
               for (int idxl = 0; idxl < lMax + 1; ++idxl) {
                  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                     // increment global
                     vec_output2(idx2) -= gsph_q0[idx3] * multfact;

                     // increment indices
                     ++idx2;
                     ++idx3;
                  }
               }
            }
         }

         // pm component
         {
            std::size_t idx1 = myidxfunc(idxelem, idxpoly, 1, -1);
            std::size_t idx3 = 0;
            auto mult1 =
                elemwidth(idxelem) * 0.5 * gaussquadweights[idxpoly] * rval;
            for (int idxl = 1; idxl < lMax + 1; ++idxl) {

               // multiplication factors
               double omegal0 = std::sqrt(static_cast<double>(idxl) *
                                          (static_cast<double>(idxl + 1)) *
                                          0.5);   // omega_l^0
               auto multfact = omegal0 * mult1;

               // loop over m
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {

                  vec_output2(idx1) -=
                      multfact *
                      (gsph_qm1[idx3] + gsph_qp1[idx3]);   // global increment

                  // increment indices
                  ++idx1;
                  ++idx3;
               }
            }
         }
      }
   }
   const double bigg_db = 6.6743 * std::pow(10.0, -11.0);
   const double pi_db = 3.1415926535;
   double multcomp =
       4.0 * pi_db * bigg_db * std::pow(backgroundmodel.LengthNorm(), 2.0);
   // std::size_t idxmaxval = std::pow(lMax + 1, 2) * matlen;
   for (auto idx1 = 0; idx1 < matlen * _size0; ++idx1) {
      vec_output2(idx1) *= backgroundmodel.DensityNorm() * multcomp;
   }

   return vec_output2;
};

template <typename FLOAT, template <typename, typename> class sphericalmodel,
          template <typename, typename> class TOMPERTMODEL,
          template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange>
Eigen::VectorXcd
AdvectiveBoundaryPerturbation(
    const sphericalmodel<FLOAT, int> &backgroundmodel,
    const TOMPERTMODEL<FLOAT, int> &pertmodel,
    GRID<FLOAT, OrderIndexRange, IndexRange> &grid,
    const CoefficientOrdering &whichsetup,
    const GaussQuad::Quadrature1D<FLOAT> &q,
    const std::vector<FLOAT> &vec_noderadii,
    const std::vector<FFTWpp::vector<std::complex<FLOAT>>> &vec_denslm,
    const int &lMax,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &mat_gaussderiv,
    const Eigen::Vector<std::complex<double>, Eigen::Dynamic> &vec_phi,
    const double &maxballrad) {
   using Real = FLOAT;
   using Complex = std::complex<Real>;
   std::size_t nelem = vec_noderadii.size() - 1;
   int npoly = q.N() - 1;
   int matlen = nelem * npoly + 1;   // size of matrix

   auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).Size();
   auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).Size();
   auto rphys = [&vec_noderadii](int idxelem, double x) {
      return GPLSpec::detail::StandardIntervalMap(
          x, vec_noderadii[idxelem],
          vec_noderadii[idxelem + 1]);
   };
   // finding \rho \mathbf{s} in canonical components spatially
   // firstly we need to have the background model density and the perturbation
   // from the perturbing model:
   // go through all layers, find what s is, as we currently just have radial,
   // that should be fairly simple
   int totnum = grid.NumberOfCoLatitudes() * grid.NumberOfLongitudes();

   int laynum = 0;

   auto myidxfunc = [npoly, lMax](int idxelem, int idxpoly, int idxl,
                                  int idxm) {
      return static_cast<std::size_t>(std::pow(idxl, 2) + idxl + idxm +
                                      (idxpoly + npoly * idxelem) *
                                          std::pow(lMax + 1, 2));
   };
   using VECTOR = Eigen::Vector<Complex, Eigen::Dynamic>;
   VECTOR vec_output2 = VECTOR::Zero(matlen * _size0);
   std::vector<double> gaussquadweights, gaussquadpoints;
   auto elemwidth = [&vec_noderadii](int idxelem) {
      return vec_noderadii[idxelem + 1] - vec_noderadii[idxelem];
   };
   for (int idx = 0; idx < q.Points().size(); ++idx) {
      gaussquadpoints.push_back(q.X(idx));
      gaussquadweights.push_back(q.W(idx));
   }

   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      double inv2 = 2.0 / elemwidth(idxelem);
      // std::cout << "inv2: " << inv2 << std::endl;
      // looping over quadrature nodes
      int idxpolymax = npoly + 1;
      int idxpolymin = 0;
      if (idxelem < nelem - 1) {
         idxpolymax = npoly;
      }
      if (idxelem == 0) {
         idxpolymin = 1;
      }
      for (int idxpoly = idxpolymin; idxpoly < idxpolymax; ++idxpoly) {
         // radius, inverse and squared
         auto rval = GPLSpec::detail::StandardIntervalMap(
             gaussquadpoints[idxpoly], vec_noderadii[idxelem],
             vec_noderadii[idxelem + 1]);
         auto invr = 1 / rval;
         auto rval2 = rval * rval;

         ///////////////////////////////////////
         // step 1: finding nabla zeta
         std::vector<std::complex<double>> gsph_nz0(_size0, 0.0),
             gsph_nzp1(_sizepm, 0.0), gsph_nzm1(_sizepm, 0.0);   // declaration

         // looping over l,m values
         {

            // cache friendly trial
            //  getting first index
            std::size_t idx1 = myidxfunc(idxelem, 0, 0, 0);

            // looping over radii
            for (int idxn = 0; idxn < npoly + 1; ++idxn) {
               std::size_t idx2 = 0;
               auto multfact = mat_gaussderiv(idxn, idxpoly);

               // looping over l and m
               for (int idxl = 0; idxl < lMax + 1; ++idxl) {
                  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {

                     // increment local temporary
                     gsph_nz0[idx2] += vec_phi(idx1) * multfact;

                     // increment indices
                     ++idx1;
                     ++idx2;
                  }
               }
            }
            for (auto &idx : gsph_nz0) {
               idx *= inv2;
            }
         }

         // only do the pm if not at the zero radius
         // if (tcount > 0) {
         {
            std::size_t idx1 = 0;
            std::size_t idx2 = myidxfunc(idxelem, idxpoly, 1, -1);
            for (int idxl = 1; idxl < lMax + 1; ++idxl) {
               // omega_l^0
               double Omegal0 = std::sqrt(static_cast<double>(idxl) *
                                          static_cast<double>(idxl + 1) / 2.0);
               auto multfact = Omegal0 * invr;
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {

                  auto tmp = vec_phi(idx2) * multfact;
                  gsph_nzp1[idx1] += tmp;
                  gsph_nzm1[idx1] += tmp;

                  // increment
                  ++idx1;
                  ++idx2;
               }
            }
         }
         // std::cout << "invr: " << invr << std::endl;
         // }

         ///////////////////////////////////////////////////////////////////////////
         // step 2: transforming (\nabla \zeta) into spatial domain
         FFTWpp::vector<Complex> spatial_nzm1(totnum, 0.0),
             spatial_nzp1(totnum, 0.0),
             spatial_nz0(totnum, 0.0);   // declaration

         // GSPH transforms:
         // lhs.gridused().InverseTransformation(lhs.lMax(), 0, gsph_nz0,
         //                                      spatial_nz0);   // 0th order
         grid.InverseTransformation(lMax, 0, gsph_nz0,
                                    spatial_nz0);   // 0th order
         grid.InverseTransformation(lMax, 1, gsph_nzp1,
                                    spatial_nzp1);   //+1 order
         grid.InverseTransformation(lMax, -1, gsph_nzm1,
                                    spatial_nzm1);   //-1 order

         ///////////////////////////////////////////////////////////////////////////
         // finding s in canonical components spatially
         std::vector<Complex> spatial_q0(totnum, 0.0), spatial_qm(totnum, 0.0),
             spatial_qp(totnum, 0.0);

         for (int idxt = 0; idxt < grid.NumberOfCoLatitudes(); ++idxt) {
            for (int idxp = 0; idxp < grid.NumberOfLongitudes(); ++idxp) {
               auto idxuse = idxt * grid.NumberOfLongitudes() + idxp;
               double radscalval;
               if (rval < backgroundmodel.OuterRadius()) {
                  spatial_q0[idxuse] = pertmodel.RadialMap(
                      rval, grid.CoLatitudes()[idxt], grid.Longitudes()[idxp]);
               } else {
                  spatial_q0[idxuse] =
                      pertmodel.RadialMap(backgroundmodel.OuterRadius(),
                                          grid.CoLatitudes()[idxt],
                                          grid.Longitudes()[idxp]) *
                      (maxballrad - rval) /
                      (maxballrad - backgroundmodel.OuterRadius());
               }
            }
         }

         ///////////////////////////////////////////////////////////////////////////
         // finding the multiple of the two
         std::vector<Complex> spatial_multtot(totnum, 0.0);
         for (int idxt = 0; idxt < grid.NumberOfCoLatitudes(); ++idxt) {
            for (int idxp = 0; idxp < grid.NumberOfLongitudes(); ++idxp) {
               auto idxuse = idxt * grid.NumberOfLongitudes() + idxp;
               spatial_multtot[idxuse] -=
                   spatial_qm[idxuse] * spatial_nzp1[idxuse];
               spatial_multtot[idxuse] +=
                   spatial_q0[idxuse] * spatial_nz0[idxuse];
               spatial_multtot[idxuse] -=
                   spatial_qp[idxuse] * spatial_nzm1[idxuse];
            }
         }

         ///////////////////////////////////////////////////////////////////////////
         // converting back into spherical harmonics
         std::vector<std::complex<double>> gsph_multtot(_size0, 0.0);
         grid.ForwardTransformation(lMax, 0, spatial_multtot, gsph_multtot);
         std::size_t idx1 = myidxfunc(idxelem, idxpoly, 0, 0);
         std::size_t idx3 = 0;
         for (int idxl = 0; idxl < lMax + 1; ++idxl) {
            // loop over m
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               vec_output2(idx1) = gsph_multtot[idx3];   // global increment

               // increment indices
               ++idx1;
               ++idx3;
            }
         }
      }
   }

   return vec_output2;
};


}   // namespace TestTools

// if (idxinner > 0) {
//             if (radius == layered_radii[laynum][idxinner - 1]) {
//                layered_radii.push_back({radius});
//                layered_rho.push_back({rho});
//                layered_vpv.push_back({vpv});
//                layered_vsv.push_back({vsv});
//                layered_qkappa.push_back({qkappa});
//                layered_qshear.push_back({qshear});
//                layered_vph.push_back({vph});
//                layered_vsh.push_back({vsh});
//                layered_eta.push_back({eta});
//                idxinner = 0;
//                ++laynum;
//             } else {
//                layered_radii[laynum].push_back(radius);
//                layered_rho[laynum].push_back(rho);
//                layered_vpv[laynum].push_back(vpv);
//                layered_vsv[laynum].push_back(vsv);
//                layered_qkappa[laynum].push_back(qkappa);
//                layered_qshear[laynum].push_back(qshear);
//                layered_vph[laynum].push_back(vph);
//                layered_vsh[laynum].push_back(vsh);
//                layered_eta[laynum].push_back(eta);
//             }
//          } else {
//             layered_radii[laynum].push_back(radius);
//             layered_rho[laynum].push_back(rho);
//             layered_vpv[laynum].push_back(vpv);
//             layered_vsv[laynum].push_back(vsv);
//             layered_qkappa[laynum].push_back(qkappa);
//             layered_qshear[laynum].push_back(qshear);
//             layered_vph[laynum].push_back(vph);
//             layered_vsh[laynum].push_back(vsh);
//             layered_eta[laynum].push_back(eta);
//          }

#endif
