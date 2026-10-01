#ifndef GPLSPEC_MAPPING_PERTURBATION_UTILITIES_H
#define GPLSPEC_MAPPING_PERTURBATION_UTILITIES_H

#include "Earth_Density_Models_3D.h"
#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

namespace GeneralEarthModels::MappingPerturbationDetail {

// Construct the canonical spatial df tensor from the existing spectral
// displacement gradients. This is shared by the radial-map and file-coefficient
// constructors; callers retain their own displacement source and storage setup.
inline void ConstructMappingVectorGradient(
    const Density3D &inp_model,
    const std::vector<std::vector<std::vector<Eigen::Vector3cd>>> &_vec_dxilm,
    std::vector<std::vector<std::vector<Eigen::Matrix3cd>>> &_vec_df,
    std::size_t _num_layers, std::size_t spatialsize, std::size_t nnode,
    std::size_t lMax) {
   auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).Size();
   auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).Size();
   auto _sizepp = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 2).Size();

   // first step is to find the gradient of dxi
   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
      // find 0-component derivative
      // components of derivative, ie \nabla h:
      using veccomp = std::vector<std::complex<double>>;
      using vecvech = std::vector<veccomp>;

      // finding 0-component derivative
      for (int idxpoly = 0; idxpoly < nnode; ++idxpoly) {
         veccomp vec_ddxim(_sizepm, 0.0);
         veccomp vec_ddxi0(_size0, 0.0);
         veccomp vec_ddxip(_sizepm, 0.0);
         veccomp vec_ddxip0(_sizepm, 0.0);
         veccomp vec_ddxim0(_sizepm, 0.0);
         veccomp vec_ddxipm(_size0, 0.0);
         veccomp vec_ddximp(_size0, 0.0);
         veccomp vec_ddxipp(_sizepp, 0.0);
         veccomp vec_ddximm(_sizepp, 0.0);
         double radr =
             inp_model.Node_InformationP().NodeRadius(idxelem, idxpoly);
         double inv2 =
             2.0 / inp_model.Node_InformationP().ElementWidth(idxelem);
         auto idxoverall = idxelem * _num_layers + idxpoly;
         // idxoverall = 1;
         // finding \partial^0 u^{\alpha}:
         {
            // looping over radii
            for (int idxn = 0; idxn < nnode; ++idxn) {
               auto multfact =
                   inp_model.GaussDerivative(idxn, idxpoly) * inv2;
               auto idxouter = idxelem * _num_layers + idxn;

               // looping over l and m
               auto idxmax = (lMax + 1) * (lMax + 1);
               vec_ddxi0[0] += _vec_dxilm[idxelem][idxn][0](1) * multfact;
               for (int idx2 = 1; idx2 < idxmax; ++idx2) {
                  vec_ddxim[idx2 - 1] +=
                      _vec_dxilm[idxelem][idxn][idx2](0) * multfact;
                  vec_ddxi0[idx2] +=
                      _vec_dxilm[idxelem][idxn][idx2](1) * multfact;
                  vec_ddxip[idx2 - 1] +=
                      _vec_dxilm[idxelem][idxn][idx2](2) * multfact;
               }
            }
         }

         // finding \partial^{\pm}u^0:
         if (idxoverall != 0) {
            auto idxmax = (lMax + 1) * (lMax + 1);
            int idx2 = 1;
            auto idxouter = idxelem * _num_layers + idxpoly;

            for (int idxl = 1; idxl < lMax + 1; ++idxl) {
               auto omegal0 =
                   std::sqrt(static_cast<double>(idxl) *
                             (static_cast<double>(idxl) + 1.0) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  auto tmp1 =
                      omegal0 * _vec_dxilm[idxelem][idxpoly][idx2][1];
                  vec_ddxim0[idx2 - 1] +=
                      (tmp1 - _vec_dxilm[idxelem][idxpoly][idx2](0)) / radr;
                  vec_ddxip0[idx2 - 1] +=
                      (tmp1 - _vec_dxilm[idxelem][idxpoly][idx2](2)) / radr;
                  ++idx2;
               }
            }
         }
         // finding \partial^{\pm}u^{\pm}:
         if (idxoverall != 0) {
            auto idxmax = (lMax + 1) * (lMax + 1);
            int idx1 = 0;
            int idx2 = 4;
            //    auto idxouter = idxelem * _num_layers + idxpoly;
            for (int idxl = 2; idxl < lMax + 1; ++idxl) {
               auto omegal2 =
                   std::sqrt((static_cast<double>(idxl) + 2.0) *
                             (static_cast<double>(idxl) - 1.0) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  // auto tmp1 = omegal0 * vec_dxi[idxouter][idx2][1];
                  vec_ddximm[idx1] +=
                      omegal2 * _vec_dxilm[idxelem][idxpoly][idx2](0) / radr;
                  vec_ddxipp[idx1] +=
                      omegal2 * _vec_dxilm[idxelem][idxpoly][idx2](2) / radr;

                  ++idx1;
                  ++idx2;
               }
            }
         }
         // finding \partial^{\pm}u^{\mp}:
         if (idxoverall != 0) {
            auto idxmax = (lMax + 1) * (lMax + 1);
            int idx2 = 0;
            //    auto idxouter = idxelem * _num_layers + idxpoly;
            for (int idxl = 0; idxl < lMax + 1; ++idxl) {
               auto omegal0 =
                   std::sqrt(static_cast<double>(idxl) *
                             (static_cast<double>(idxl) + 1.0) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  auto tmp1 =
                      omegal0 * _vec_dxilm[idxelem][idxpoly][idx2](0);
                  vec_ddxipm[idx2] +=
                      (omegal0 * _vec_dxilm[idxelem][idxpoly][idx2](0) -
                       _vec_dxilm[idxelem][idxpoly][idx2](1)) /
                      radr;
                  vec_ddximp[idx2] +=
                      (omegal0 * _vec_dxilm[idxelem][idxpoly][idx2](2) -
                       _vec_dxilm[idxelem][idxpoly][idx2](1)) /
                      radr;
                  ++idx2;
               }
            }
         }

         /////////////////////////////////////////////////////////////////
         // declare spatial variables
         veccomp vec_ddxim_spatial(spatialsize, 0.0);
         veccomp vec_ddxi0_spatial(spatialsize, 0.0);
         veccomp vec_ddxip_spatial(spatialsize, 0.0);
         veccomp vec_ddxip0_spatial(spatialsize, 0.0);
         veccomp vec_ddxim0_spatial(spatialsize, 0.0);
         veccomp vec_ddxipm_spatial(spatialsize, 0.0);
         veccomp vec_ddximp_spatial(spatialsize, 0.0);
         veccomp vec_ddxipp_spatial(spatialsize, 0.0);
         veccomp vec_ddximm_spatial(spatialsize, 0.0);

         // transforming
         // 00
         inp_model.GSH_GridP().InverseTransformation(lMax, 0, vec_ddxi0,
                                                     vec_ddxi0_spatial);

         // 0\pm
         inp_model.GSH_GridP().InverseTransformation(lMax, -1, vec_ddxim,
                                                     vec_ddxim_spatial);
         inp_model.GSH_GridP().InverseTransformation(lMax, +1, vec_ddxip,
                                                     vec_ddxip_spatial);

         //\pm 0
         inp_model.GSH_GridP().InverseTransformation(lMax, -1, vec_ddxim0,
                                                     vec_ddxim0_spatial);
         inp_model.GSH_GridP().InverseTransformation(lMax, +1, vec_ddxip0,
                                                     vec_ddxip0_spatial);

         //\pm \mp
         inp_model.GSH_GridP().InverseTransformation(lMax, 0, vec_ddxipm,
                                                     vec_ddxipm_spatial);
         inp_model.GSH_GridP().InverseTransformation(lMax, 0, vec_ddximp,
                                                     vec_ddximp_spatial);

         //\pm \pm
         inp_model.GSH_GridP().InverseTransformation(lMax, -2, vec_ddximm,
                                                     vec_ddximm_spatial);
         inp_model.GSH_GridP().InverseTransformation(lMax, +2, vec_ddxipp,
                                                     vec_ddxipp_spatial);

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
            for (auto it : inp_model.GSH_GridP().CoLatitudes()) {
               for (auto ip : inp_model.GSH_GridP().Longitudes()) {
                  // going along first row (and transposing)
                  _vec_df[idxelem][idxpoly][idxvec](0, 0) =
                      vec_ddximm_spatial[idxvec];
                  _vec_df[idxelem][idxpoly][idxvec](0, 1) =
                      vec_ddxim_spatial[idxvec];
                  _vec_df[idxelem][idxpoly][idxvec](0, 2) =
                      vec_ddxipm_spatial[idxvec];

                  // going along first row (and transposing)
                  _vec_df[idxelem][idxpoly][idxvec](1, 0) =
                      vec_ddxim0_spatial[idxvec];
                  _vec_df[idxelem][idxpoly][idxvec](1, 1) =
                      vec_ddxi0_spatial[idxvec];
                  _vec_df[idxelem][idxpoly][idxvec](1, 2) =
                      vec_ddxip0_spatial[idxvec];

                  // going along first row (and transposing)
                  _vec_df[idxelem][idxpoly][idxvec](2, 0) =
                      vec_ddximp_spatial[idxvec];
                  _vec_df[idxelem][idxpoly][idxvec](2, 1) =
                      vec_ddxip_spatial[idxvec];
                  _vec_df[idxelem][idxpoly][idxvec](2, 2) =
                      vec_ddxipp_spatial[idxvec];

                  ++idxvec;
               }
            }
         }
      }
   }
}

}  // namespace GeneralEarthModels::MappingPerturbationDetail

#endif
