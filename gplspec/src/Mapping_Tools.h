#ifndef MAPPING_TOOLS_H
#define MAPPING_TOOLS_H

#include <Eigen/Core>
#include <Eigen/Dense>
// #include <Eigen/CholmodSupport>
#include <Eigen/Sparse>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <PlanetaryModel/All>
// #include <TomographyModels/All>
#include "SphericalGeometryPreconditioner.h"
#include "Spherical_Integrator.h"
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

namespace MappingTools {
using namespace GravityFunctions;

using MATRIX3cd = Eigen::Matrix<std::complex<double>, 3, 3>;
using EARTHMATRIX3 = std::vector<std::vector<MATRIX3cd>>;
using EARTHVEC = std::vector<std::vector<std::vector<double>>>;
using RADIUSVEC = std::vector<std::vector<std::complex<double>>>;
using EARTHVECX = std::vector<RADIUSVEC>;

template <typename FLOAT, template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange,
          template <typename, typename> class TOMMODEL,
          template <typename, typename> class TOMPERTMODEL>
auto
tom_to_dxi(GRID<FLOAT, OrderIndexRange, IndexRange> &grid, std::size_t npoly,
           std::vector<double> &vec_noderadii,
           std::vector<double> &vec_allradii,
           const GaussQuad::Quadrature1D<FLOAT> &q,
           const TOMPERTMODEL<FLOAT, int> &mypertprem,
           const TOMMODEL<FLOAT, int> &backgroundmodel) {

   auto mytheta = grid.CoLatitudes();
   auto myphi = grid.Longitudes();
   auto intsize = grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();
   auto nelem = vec_noderadii.size() - 1;
   auto ballrad = vec_noderadii.back();

   // std::vector<std::vector<std::complex<double>>> vec_h(nelem * npoly + 1);
   std::vector<std::complex<double>> vec0(3, 0.0);
   MappingTools::EARTHVECX vec_dxi(nelem * npoly + 1,
                                   MappingTools::RADIUSVEC(intsize, vec0));

   FLOAT outerrad = backgroundmodel.OuterRadius();
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      int idxpolymax = npoly;
      if (idxelem == nelem - 1) {
         idxpolymax = npoly + 1;
      }
      auto mytheta = grid.CoLatitudes();
      auto myphi = grid.Longitudes();
      auto normval = std::sqrt(4.0 * 3.1415926535);
      // std::cout << "Hello internal\n";
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         auto radr = GravityFunctions::StandardIntervalMap(
             q.X(idxpoly), vec_noderadii[idxelem], vec_noderadii[idxelem + 1]);

         for (int idxtheta = 0; idxtheta < grid.NumberOfCoLatitudes();
              ++idxtheta) {
            for (int idxphi = 0; idxphi < grid.NumberOfLongitudes(); ++idxphi) {
               if (vec_allradii[idxelem * npoly + idxpoly] <
                   backgroundmodel.OuterRadius()) {
                  vec_dxi[idxelem * npoly + idxpoly]
                         [idxtheta * grid.NumberOfLongitudes() + idxphi][1] =
                             mypertprem.RadialMap(radr, mytheta[idxtheta],
                                                  myphi[idxphi]);
               } else {
                  double mytmp =
                      mypertprem.RadialMap(backgroundmodel.OuterRadius(),
                                           mytheta[idxtheta], myphi[idxphi]);
                  auto idx1 = idxelem * npoly + idxpoly;
                  auto idx2 = idxtheta * grid.NumberOfLongitudes() + idxphi;
                  vec_dxi[idx1][idx2][1] =
                      mytmp * (ballrad - vec_allradii[idx1]) /
                      (ballrad - backgroundmodel.OuterRadius());
               }
            }
         }
      }
   }
   return vec_dxi;
};

template <typename FLOAT, template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange>
auto
dxi_to_dxilm(
    GRID<FLOAT, OrderIndexRange, IndexRange> &grid, std::size_t npoly,
    std::vector<double> &vec_elemwidth, std::vector<double> &vec_allradii,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &_mat_gaussderiv,
    const EARTHVECX &vec_dxi) {
   // vec_dxi is in the form of a vector of vectors
   // it contains the (-,0,+) components of dxi spatially

   // various sizes and definitions
   auto lMax = grid.MaxDegree();
   auto nMax = grid.MaxUpperIndex();
   auto size =
       GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();   // dof
   auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
   auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
   auto nelem = vec_elemwidth.size();
   int matlen = nelem * npoly + 1;   // size of matrix
   auto intsize = grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();

   auto mat_0 = MATRIX3cd::Zero();
   std::vector<std::complex<double>> vec_0(3, 0.0);
   //    EARTHMATRIX3 vec_da(nelem * (npoly + 1),
   //                        std::vector<MATRIX3cd>(intsize, mat_0));
   //    EARTHMATRIX3 vec_Ddxi(nelem * (npoly + 1),
   //                          std::vector<MATRIX3cd>(intsize, mat_0));
   EARTHVECX vec_dxilm(nelem * npoly + 1, RADIUSVEC(intsize, vec_0));

   // first step is to find the gradient of dxi
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      int idxpolymax = npoly;
      if (idxelem == nelem - 1) {
         idxpolymax = npoly + 1;
      }
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         std::vector<std::complex<double>> vec_dxim(intsize, 0.0),
             vec_dxi0(intsize, 0.0), vec_dxip(intsize, 0.0);

         for (int idxint = 0; idxint < intsize; ++idxint) {
            vec_dxim[idxint] = vec_dxi[idxelem * npoly + idxpoly][idxint][0];
            vec_dxi0[idxint] = vec_dxi[idxelem * npoly + idxpoly][idxint][1];
            vec_dxip[idxint] = vec_dxi[idxelem * npoly + idxpoly][idxint][2];
         }

         // transform
         std::vector<std::complex<double>> vec_lm_dxim(_sizepm, 0.0),
             vec_lm_dxi0(_size0, 0.0), vec_lm_dxip(_sizepm, 0.0);

         grid.ForwardTransformation(lMax, -1, vec_dxim, vec_lm_dxim);
         grid.ForwardTransformation(lMax, 0, vec_dxi0, vec_lm_dxi0);
         grid.ForwardTransformation(lMax, 1, vec_dxip, vec_lm_dxip);

         auto intoutput = idxelem * npoly + idxpoly;
         vec_dxilm[intoutput][0][1] = vec_lm_dxi0[0];
         for (int idxl = 1; idxl < lMax + 1; ++idxl) {
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               auto idxlm = idxl * idxl + idxl + idxm;
               vec_dxilm[intoutput][idxlm][0] = vec_lm_dxim[idxlm - 1];
               vec_dxilm[intoutput][idxlm][1] = vec_lm_dxi0[idxlm];
               vec_dxilm[intoutput][idxlm][2] = vec_lm_dxip[idxlm - 1];
            }
         }
      }
   }

   return vec_dxilm;
}

template <typename FLOAT, template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange>
auto
dxitodf(GRID<FLOAT, OrderIndexRange, IndexRange> &grid, std::size_t npoly,
        std::vector<double> &vec_elemwidth, std::vector<double> &vec_allradii,
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &_mat_gaussderiv,
        std::vector<std::vector<std::complex<double>>> &vec_h,
        const EARTHMATRIX3 &vec_a, const EARTHVECX &vec_dxi) {
   // vec_dxi is in the form of a vector of vectors
   // it contains the (-,0,+) components of dxi in lm components

   // various sizes and definitions
   auto lMax = grid.MaxDegree();
   auto nMax = grid.MaxUpperIndex();
   auto size =
       GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();   // dof
   auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
   auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
   auto _sizepp = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 2).size();
   auto nelem = vec_elemwidth.size();
   int matlen = nelem * npoly + 1;   // size of matrix
   auto intsize = grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();

   auto mat_0 = MATRIX3cd::Zero();
   EARTHMATRIX3 vec_da(nelem * (npoly + 1),
                       std::vector<MATRIX3cd>(intsize, mat_0));
   EARTHMATRIX3 vec_df(nelem * (npoly + 1),
                       std::vector<MATRIX3cd>(intsize, mat_0));

   // first step is to find the gradient of dxi
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // find 0-component derivative
      // components of derivative, ie \nabla h:
      using veccomp = std::vector<std::complex<double>>;
      using vecvech = std::vector<veccomp>;

      // finding 0-component derivative
      for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
         veccomp vec_ddxim(_sizepm, 0.0);
         veccomp vec_ddxi0(_size0, 0.0);
         veccomp vec_ddxip(_sizepm, 0.0);
         veccomp vec_ddxip0(_sizepm, 0.0);
         veccomp vec_ddxim0(_sizepm, 0.0);
         veccomp vec_ddxipm(_size0, 0.0);
         veccomp vec_ddximp(_size0, 0.0);
         veccomp vec_ddxipp(_sizepp, 0.0);
         veccomp vec_ddximm(_sizepp, 0.0);

         double inv2 = 2.0 / vec_elemwidth[idxelem];
         auto idxoverall = idxelem * nelem + idxpoly;
         // idxoverall = 1;
         // finding \partial^0 u^{\alpha}:
         {
            // looping over radii
            for (int idxn = 0; idxn < npoly + 1; ++idxn) {
               auto multfact = _mat_gaussderiv(idxn, idxpoly) * inv2;
               auto idxouter = idxelem * nelem + idxn;

               // looping over l and m
               auto idxmax = (lMax + 1) * (lMax + 1);
               vec_ddxi0[0] += vec_dxi[idxouter][0][1] * multfact;
               for (int idx2 = 1; idx2 < idxmax; ++idx2) {
                  vec_ddxim[idx2 - 1] += vec_dxi[idxouter][idx2][0] * multfact;
                  vec_ddxi0[idx2] += vec_dxi[idxouter][idx2][1] * multfact;
                  vec_ddxip[idx2 - 1] += vec_dxi[idxouter][idx2][2] * multfact;
               }
            }
         }

         // finding \partial^{\pm}u^0:
         if (idxoverall != 0) {
            auto idxmax = (lMax + 1) * (lMax + 1);
            int idx2 = 1;
            auto idxouter = idxelem * nelem + idxpoly;
            for (int idxl = 1; idxl < lMax + 1; ++idxl) {
               auto omegal0 =
                   std::sqrt(static_cast<double>(idxl) *
                             (static_cast<double>(idxl) + 1.0) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  auto tmp1 = omegal0 * vec_dxi[idxouter][idx2][1];
                  vec_ddxim0[idx2 - 1] += (tmp1 - vec_dxi[idxouter][idx2][0]) /
                                          vec_allradii[idxouter];
                  vec_ddxip0[idx2 - 1] += (tmp1 - vec_dxi[idxouter][idx2][2]) /
                                          vec_allradii[idxouter];
                  ++idx2;
               }
            }
         }
         // finding \partial^{\pm}u^{\pm}:
         if (idxoverall != 0) {
            auto idxmax = (lMax + 1) * (lMax + 1);
            int idx1 = 0;
            int idx2 = 4;
            auto idxouter = idxelem * nelem + idxpoly;
            for (int idxl = 2; idxl < lMax + 1; ++idxl) {
               auto omegal2 =
                   std::sqrt((static_cast<double>(idxl) + 2.0) *
                             (static_cast<double>(idxl) - 1.0) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  // auto tmp1 = omegal0 * vec_dxi[idxouter][idx2][1];
                  vec_ddximm[idx1] += omegal2 * vec_dxi[idxouter][idx2][0] /
                                      vec_allradii[idxouter];
                  vec_ddxipp[idx1] += omegal2 * vec_dxi[idxouter][idx2][2] /
                                      vec_allradii[idxouter];

                  ++idx1;
                  ++idx2;
               }
            }
         }
         // finding \partial^{\pm}u^{\mp}:
         if (idxoverall != 0) {
            auto idxmax = (lMax + 1) * (lMax + 1);
            int idx2 = 0;
            auto idxouter = idxelem * nelem + idxpoly;
            for (int idxl = 0; idxl < lMax + 1; ++idxl) {
               auto omegal0 =
                   std::sqrt(static_cast<double>(idxl) *
                             (static_cast<double>(idxl) + 1.0) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  auto tmp1 = omegal0 * vec_dxi[idxouter][idx2][0];
                  vec_ddxipm[idx2] += (omegal0 * vec_dxi[idxouter][idx2][0] -
                                       vec_dxi[idxouter][idx2][1]) /
                                      vec_allradii[idxouter];
                  vec_ddximp[idx2] += (omegal0 * vec_dxi[idxouter][idx2][2] -
                                       vec_dxi[idxouter][idx2][1]) /
                                      vec_allradii[idxouter];
                  ++idx2;
               }
            }
         }

         /////////////////////////////////////////////////////////////////
         // declare spatial variables
         veccomp vec_ddxim_spatial(intsize, 0.0);
         veccomp vec_ddxi0_spatial(intsize, 0.0);
         veccomp vec_ddxip_spatial(intsize, 0.0);
         veccomp vec_ddxip0_spatial(intsize, 0.0);
         veccomp vec_ddxim0_spatial(intsize, 0.0);
         veccomp vec_ddxipm_spatial(intsize, 0.0);
         veccomp vec_ddximp_spatial(intsize, 0.0);
         veccomp vec_ddxipp_spatial(intsize, 0.0);
         veccomp vec_ddximm_spatial(intsize, 0.0);

         // transforming
         // 00
         grid.InverseTransformation(lMax, 0, vec_ddxi0, vec_ddxi0_spatial);

         // 0\pm
         grid.InverseTransformation(lMax, -1, vec_ddxim, vec_ddxim_spatial);
         grid.InverseTransformation(lMax, +1, vec_ddxip, vec_ddxip_spatial);

         //\pm 0
         grid.InverseTransformation(lMax, -1, vec_ddxim0, vec_ddxim0_spatial);
         grid.InverseTransformation(lMax, +1, vec_ddxip0, vec_ddxip0_spatial);

         //\pm \mp
         grid.InverseTransformation(lMax, 0, vec_ddxipm, vec_ddxipm_spatial);
         grid.InverseTransformation(lMax, 0, vec_ddximp, vec_ddximp_spatial);

         //\pm \pm
         grid.InverseTransformation(lMax, -2, vec_ddximm, vec_ddximm_spatial);
         grid.InverseTransformation(lMax, +2, vec_ddxipp, vec_ddxipp_spatial);

         // finding pm derivative of 0-order part
         //   vecvech vec_ddxip0(npoly + 1, veccomp(_sizepm), 0.0);
         //   vecvech vec_ddxim0(npoly + 1, veccomp(_sizepm), 0.0);
         //   vecvech vec_ddxipm(npoly + 1, veccomp(_size0), 0.0);
         //   vecvech vec_ddximp(npoly + 1, veccomp(_size0), 0.0);
         //   vecvech vec_ddxipp(npoly + 1, veccomp(_sizepp), 0.0);
         //   vecvech vec_ddximm(npoly + 1, veccomp(_sizepp), 0.0);

         /////////////////////////////////////////////////////////////////
         // filling out dxi
         {
            std::size_t idxvec = 0;
            for (int idxl1 = 0; idxl1 < grid.NumberOfCoLatitudes(); ++idxl1) {
               for (int idxl2 = 0; idxl2 < grid.NumberOfLongitudes(); ++idxl2) {
                  // going along first row (and transposing)
                  vec_df[idxelem * (npoly + 1) + idxpoly][idxvec](0, 0) =
                      vec_ddximm_spatial[idxvec];
                  vec_df[idxelem * (npoly + 1) + idxpoly][idxvec](0, 1) =
                      vec_ddxim_spatial[idxvec];
                  vec_df[idxelem * (npoly + 1) + idxpoly][idxvec](0, 2) =
                      vec_ddxipm_spatial[idxvec];

                  // going along first row (and transposing)
                  vec_df[idxelem * (npoly + 1) + idxpoly][idxvec](1, 0) =
                      vec_ddxim0_spatial[idxvec];
                  vec_df[idxelem * (npoly + 1) + idxpoly][idxvec](1, 1) =
                      vec_ddxi0_spatial[idxvec];
                  vec_df[idxelem * (npoly + 1) + idxpoly][idxvec](1, 2) =
                      vec_ddxip0_spatial[idxvec];

                  // going along first row (and transposing)
                  vec_df[idxelem * (npoly + 1) + idxpoly][idxvec](2, 0) =
                      vec_ddximp_spatial[idxvec];
                  vec_df[idxelem * (npoly + 1) + idxpoly][idxvec](2, 1) =
                      vec_ddxip_spatial[idxvec];
                  vec_df[idxelem * (npoly + 1) + idxpoly][idxvec](2, 2) =
                      vec_ddxipp_spatial[idxvec];

                  ++idxvec;
               }
            }
         }
      }
   }
   return vec_df;
}

template <typename FLOAT, template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange>
auto
h_to_f(GRID<FLOAT, OrderIndexRange, IndexRange> &grid, std::size_t npoly,
       std::vector<double> &vec_elemwidth, std::vector<double> &vec_allradii,
       Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &_mat_gaussderiv,
       std::vector<std::vector<std::complex<double>>> &vec_h) {

   // various sizes and definitions
   auto lMax = grid.MaxDegree();
   auto nMax = grid.MaxUpperIndex();
   auto size =
       GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();   // dof
   auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
   auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
   auto nelem = vec_elemwidth.size();
   int matlen = nelem * npoly + 1;   // size of matrix
   auto intsize = grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();

   // declaring typenames
   using veccomp = std::vector<std::complex<double>>;
   using vecvech = std::vector<veccomp>;
   using MATRIX3 = Eigen::Matrix<std::complex<double>, 3, 3>;

   // declaring variables used in calculation
   // vecvech vec_hlm(matlen, veccomp(_size0, 0.0));
   auto vec_hlm = vecgrid_to_sphdecomp(grid, vec_h, npoly, nelem);
   auto vec_hcheck = sphdecomp_to_vecgrid(grid, vec_hlm, npoly, nelem);

   // std::cout << std::endl;
   auto mat_0 = MATRIX3::Zero();
   std::vector<std::vector<MATRIX3>> vec_f(
       nelem * (npoly + 1), std::vector<MATRIX3>(intsize, mat_0));

   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // max for idxpoly loop
      int idxpolymax = npoly + 1;
      int idxpolymin = 0;
      if (idxelem == 0) {
         idxpolymin = 1;
      }

      // components of derivative, ie \nabla h:
      vecvech vec_hp(npoly + 1, veccomp(_sizepm, 0.0)),
          vec_hm(npoly + 1, veccomp(_sizepm, 0.0));
      double inv2 = 2.0 / vec_elemwidth[idxelem];

      // finding h0 component
      vecvech vec_h02(npoly + 1, veccomp(intsize, 0.0));

      // finding h0 component
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         for (int idxrad = 0; idxrad < intsize; ++idxrad) {
            for (int idxn = 0; idxn < npoly + 1; ++idxn) {
               auto idxouter = idxelem * npoly + idxn;
               vec_h02[idxpoly][idxrad] += vec_h[idxouter][idxrad] *
                                           _mat_gaussderiv(idxn, idxpoly) *
                                           inv2;
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
                   vec_hlm[idxouter][idx0] * Omegal0 / vec_allradii[idxouter];
               vec_hm[idxpoly][idxref] = vec_hp[idxpoly][idxref];
            }
         }
      }

      // transform back into spatial
      vecvech vec_nhp(npoly + 1, veccomp(intsize, 0.0)),
          vec_nhm(npoly + 1, veccomp(intsize, 0.0));

      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         grid.InverseTransformation(lMax, 1, vec_hp[idxpoly], vec_nhp[idxpoly]);
         grid.InverseTransformation(lMax, -1, vec_hm[idxpoly],
                                    vec_nhm[idxpoly]);
      }

      // find a
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         for (int idxrad = 0; idxrad < intsize; ++idxrad) {
            auto idxuse = idxelem * (npoly + 1) + idxpoly;
            auto nhm = vec_nhm[idxpoly][idxrad];
            auto nhp = vec_nhp[idxpoly][idxrad];
            // auto nh0 = vec_nh0[idxpoly][idxrad];
            auto nh0 = vec_h02[idxpoly][idxrad];
            std::complex<double> hdivr;

            if (idxuse == 0) {
               hdivr = 0.0;
            } else {
               hdivr = vec_h[idxelem * npoly + idxpoly][idxrad] /
                       vec_allradii[idxelem * npoly + idxpoly];
            }

            auto tmp1 = 1.0 / (1.0 + hdivr);
            auto tmp2 = tmp1 / (1.0 + nh0);

            if (std::isnan(std::abs(tmp1))) {
               std::cout << "idxe: " << idxelem << ", idxp: " << idxpoly
                         << std::endl;
            }

            vec_f[idxuse][idxrad](0, 2) = -tmp1;
            vec_f[idxuse][idxrad](1, 0) = -nhm * tmp2;
            vec_f[idxuse][idxrad](1, 1) = tmp1 - (nh0 - hdivr) * tmp2;
            vec_f[idxuse][idxrad](1, 2) = -nhp * tmp2;
            vec_f[idxuse][idxrad](2, 0) = -tmp1;

            // if (idxelem == 0 && idxpoly == 0 && idxrad == 0) {
            //    std::cout << "Hello\n\n";
            //    std::cout << tmp1 << " " << tmp2 << " " << nhm << " " << nhp
            //              << " " << hdivr << std::endl;
            //    std::cout << vec_f[0][0] << std::endl;
            //    std::cout << "Hello 2\n\n";
            // }
         }
      }
   }
   // std::cout << vec_f[0][0] << std::endl;
   return vec_f;
};

auto
df_to_da(const EARTHMATRIX3 &vec_a, const EARTHMATRIX3 &vec_invf,
         const EARTHMATRIX3 &vec_df) {
   auto numnodes = vec_a.size();
   auto numspatial = vec_a[0].size();

   EARTHMATRIX3 vec_da(numnodes,
                       std::vector<MATRIX3cd>(numspatial, MATRIX3cd::Zero()));
   MATRIX3cd mat_metric;
   mat_metric(0, 2) = -1.0;
   mat_metric(1, 1) = 1.0;
   mat_metric(2, 0) = -1.0;

   // loop through
   for (int idxouter = 0; idxouter < numnodes; ++idxouter) {
      for (int idxinner = 0; idxinner < numspatial; ++idxinner) {
         MATRIX3cd tmp1 = vec_invf[idxouter][idxinner] * mat_metric *
                          vec_df[idxouter][idxinner];
         vec_da[idxouter][idxinner] += vec_a[idxouter][idxinner] *
                                       (-tmp1(2, 0) + tmp1(1, 1) - tmp1(0, 2));
         MATRIX3cd tmp2 = tmp1 * mat_metric * vec_a[idxouter][idxinner];
         vec_da[idxouter][idxinner] -= tmp2;
         vec_da[idxouter][idxinner] -= tmp2.transpose();
      }
   }
   return vec_da;
};

// // function for finding the force vector with perturbation to mapping
// template <typename FLOAT, template <typename, typename> class sphericalmodel>
// Eigen::VectorXcd
// FindForce(
//     const sphericalmodel<FLOAT, int> &backgroundmodel,
//     const CoefficientOrdering &whichsetup,
//     const GaussQuad::Quadrature1D<FLOAT> &q, std::vector<double>
//     &vec_elemwidth, const std::vector<FLOAT> &vec_noderadii, const int &lMax,
//     const Eigen::VectorXcd &vec_sol0, const EARTHMATRIX3 &vec_da,
//     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &_mat_gaussderiv) {
//    // auto intsize = grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();
//    // lambda to use in mapping [-1,1] to [idxelem]
//    auto rscale = [&vec_noderadii](int idxelem, double x) {
//       return (x + (vec_noderadii[idxelem + 1] + vec_noderadii[idxelem]) /
//                       (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]));
//    };

//    std::size_t nelem = vec_noderadii.size() - 1;
//    int npoly = q.N() - 1;
//    auto myidxfunc = [&npoly, &lMax](int idxelem, int idxpoly, int idxl,
//                                     int idxm) {
//       return static_cast<int>(idxl * idxl + idxl + idxm +
//                               (idxpoly + npoly * idxelem) *
//                                   std::pow(lMax + 1, 2));
//    };
//    int matlen = nelem * npoly + 1;   // size of matrix
//    Eigen::VectorXcd vec_fullforce =
//        Eigen::VectorXcd::Zero(std::pow(lMax + 1, 2) * matlen);

//    int laynum = 0;   // layer number
//    // constants
//    const double bigg_db = 6.6743 * std::pow(10.0, -11.0);
//    const double pi_db = 3.1415926535;
//    double multfact =
//        4.0 * pi_db * bigg_db * std::pow(backgroundmodel.LengthNorm(), 2.0);

//    auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
//    auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();

//    for (int idxelem = 0; idxelem < nelem; ++idxelem) {
//       double inv2 = 2.0 / vec_elemwidth[idxelem];
//       // looping over quadrature nodes
//       int idxpolymax = npoly + 1;
//       for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
//          // radius, inverse and squared
//          auto rval = GravityFunctions::StandardIntervalMap(
//              q.X(idxpoly), vec_noderadii[idxelem], vec_noderadii[idxelem +
//              1]);
//          auto invr = 1 / rval;
//          auto rval2 = rval * rval;

//          ///////////////////////////////////////
//          // step 1: finding nabla zeta
//          std::vector<std::complex<double>> gsph_nz0(_size0, 0.0),
//              gsph_nzp1(_size0, 0.0), gsph_nzm1(_sizepm, 0.0);   //
//              declaration

//          // looping over l,m values
//          {

//             // cache friendly trial
//             //  getting first index
//             std::size_t idx1 = myidxfunc(idxelem, 0, 0, 0);

//             // looping over radii
//             for (int idxn = 0; idxn < npoly + 1; ++idxn) {
//                std::size_t idx2 = 0;
//                auto multfact = _mat_gaussderiv(idxn, idxpoly);

//                // looping over l and m
//                for (int idxl = 0; idxl < lMax + 1; ++idxl) {
//                   for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {

//                      // increment local temporary
//                      gsph_nz0[idx2] += vec_sol0(idx1) * multfact;

//                      // increment indices
//                      ++idx1;
//                      ++idx2;
//                   }
//                }
//             }
//             // for (auto &idx : gsph_nz0) {
//             //    idx *= inv2;
//             // }
//          }

//          // only do the pm if not at the zero radius
//          // if (tcount > 0) {
//          {
//             std::size_t idx1 = 0;
//             std::size_t idx2 = myidxfunc(idxelem, idxpoly, 1, -1);
//             for (int idxl = 1; idxl < lMax + 1; ++idxl) {
//                // omega_l^0
//                double Omegal0 = std::sqrt(static_cast<double>(idxl) *
//                                           static_cast<double>(idxl + 1)
//                                           / 2.0);
//                auto multfact = Omegal0 * invr;
//                for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {

//                   auto tmp = vec_sol0(idx2) * multfact;
//                   gsph_nzp1[idx1] += tmp;
//                   gsph_nzm1[idx1] += tmp;

//                   // increment
//                   ++idx1;
//                   ++idx2;
//                }
//             }
//          }
//       }
//    }

//    // looping over elements
//    // for (int idxelem = 0; idxelem < nelem; ++idxelem) {
//    //    // finding the layer number
//    //    // if (laynum < backgroundmodel.NumberOfLayers()) {
//    //    //    if (!(vec_noderadii[idxelem] <
//    //    backgroundmodel.UpperRadius(laynum)))
//    //    //    {
//    //    //       laynum += 1;
//    //    //    };
//    //    // }

//    //    // looping over quadrature nodes
//    //    for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
//    //       // looping over l,m values
//    //       for (int idxl = 0; idxl < lMax + 1; ++idxl) {
//    //          for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {

//    //             std::size_t mynewidx =
//    //                 std::pow(idxl, 2) + idxl + idxm +
//    //                 (idxpoly + npoly * idxelem) * std::pow(lMax + 1, 2);

//    //             vec_fullforce(mynewidx) +=
//    //                 multfact * q.W(idxpoly) *
//    //                 std::pow(0.5 * (vec_noderadii[idxelem + 1] -
//    //                                 vec_noderadii[idxelem]),
//    //                          3.0) *
//    //                 std::pow(rscale(idxelem, q.X(idxpoly)), 2.0) *
//    //                 vec_denslm[idxelem * (npoly + 1) + idxpoly]
//    //                           [std::pow(idxl, 2) + idxl + idxm];
//    //          }
//    //       }
//    //    };
//    // };
//    // for (int idxl = 1; idxl < lMax + 1; ++idxl) {
//    //    for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//    //       std::size_t mynewidx = std::pow(idxl, 2) + idxl + idxm;
//    //       vec_fullforce(mynewidx) *= 0;
//    //    }
//    // }
//    // std::size_t idxmaxval = std::pow(lMax + 1, 2) * matlen;
//    // for (auto idx1 = 0; idx1 < idxmaxval; ++idx1) {
//    //    vec_fullforce(idx1) *= backgroundmodel.DensityNorm();
//    // }

//    return vec_fullforce;
// };

}   // namespace MappingTools

#endif