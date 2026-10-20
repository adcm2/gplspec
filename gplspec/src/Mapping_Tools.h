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
#include "Shared_Utilities.h"
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

using MATRIX3cd = Eigen::Matrix<std::complex<double>, 3, 3>;
using EARTHMATRIX3 = std::vector<std::vector<MATRIX3cd>>;
using EARTHVEC = std::vector<std::vector<std::vector<double>>>;
using RADIUSVEC = std::vector<std::vector<std::complex<double>>>;
using EARTHVECX = std::vector<RADIUSVEC>;

template <typename FLOAT, template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange>
inline auto
h_to_f(GRID<FLOAT, OrderIndexRange, IndexRange> &grid, std::size_t npoly,
       std::vector<double> &vec_elemwidth, std::vector<double> &vec_allradii,
       Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &_mat_gaussderiv,
       std::vector<std::vector<std::complex<double>>> &vec_h) {

   // various sizes and definitions
   auto lMax = grid.MaxDegree();
   auto nMax = grid.MaxUpperIndex();
   auto size =
       GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).Size();   // dof
   auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).Size();
   auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).Size();
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

inline auto
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
