#ifndef PERTURBATION_MODEL_3D_H
#define PERTURBATION_MODEL_3D_H

// #include "../SimpleModels"
#include "../Radial_Tools.h"
// #include "Simple_Models.h"
#include "../Spectral_Element_Tools.h"
#include "../testtools.h"
#include "Earth_Density_Models_3D.h"
#include "Density_Model_Constructors.h"
#include "Density_Model_Return.h"
#include <GSHTrans/All>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <PlanetaryModel/All>
#include <TomographyModels/All>

namespace GeneralEarthModels {
class MappingPerturbation {

 private:
   using vecdb = std::vector<double>;
   using vvecdb = std::vector<vecdb>;
   using vvvecdb = std::vector<vvecdb>;
   using veceigvec = std::vector<Eigen::Vector3cd>;
   using vveceigvec = std::vector<veceigvec>;
   using vvveceigvec = std::vector<vveceigvec>;
   using veceig = std::vector<Eigen::Matrix3cd>;
   using vveceig = std::vector<veceig>;
   using vvveceig = std::vector<vveceig>;
   vvveceigvec _vec_dxi, _vec_dxilm;
   vvveceig _vec_df, _vec_da;
   const Density3D &_ref_model;

 public:
   // no perturbation
   MappingPerturbation(const Density3D &);

   // radial map only
   template <class mapclass>
      requires PlanetaryModel::RadialMappingClass<mapclass>
   MappingPerturbation(const Density3D &, const mapclass &);

   // radial map with file input
   MappingPerturbation(const Density3D &, const std::string &, const int,
                       const int);

   auto dxi() const { return _vec_dxi; };
   auto dxilm() const { return _vec_dxilm; };
   auto df() const { return _vec_df; };
   auto da() const { return _vec_da; };
   const vvveceig &ref_da() const { return _vec_da; };
   const Density3D &ref_inp_model() const { return _ref_model; };
};

// no map constructor
MappingPerturbation::MappingPerturbation(const Density3D &inp_model)
    : _ref_model{inp_model} {
   std::size_t _num_layers = inp_model.Num_Elements();
   std::size_t spatialsize = inp_model.GSH_Grid().NumberOfCoLatitudes() *
                             inp_model.GSH_Grid().NumberOfLongitudes();
   std::size_t lMax = inp_model.GSH_Grid().MaxDegree();

   // construct dxi = 0
   _vec_dxi = vvveceigvec(
       _num_layers,
       vveceigvec(inp_model.q().N(),
                  veceigvec(spatialsize, Eigen::Vector3cd::Zero())));
   _vec_dxilm = vvveceigvec(
       _num_layers,
       vveceigvec(inp_model.q().N(), veceigvec((lMax + 1) * (lMax + 1),
                                               Eigen::Vector3cd::Zero())));
   // construct da = 0
   Eigen::Matrix3cd mat_a0;
   mat_a0 = Eigen::Matrix3cd::Zero(3, 3);
   _vec_da = vvveceig(_num_layers,
                      vveceig(inp_model.q().N(), veceig(spatialsize, mat_a0)));
};

// radial map only
template <class mapclass>
   requires PlanetaryModel::RadialMappingClass<mapclass>
MappingPerturbation::MappingPerturbation(const Density3D &inp_model,
                                         const mapclass &inp_map)
    : _ref_model{inp_model} {
   std::size_t _num_layers = inp_model.Num_Elements();
   std::size_t spatialsize = inp_model.GSH_Grid().NumberOfCoLatitudes() *
                             inp_model.GSH_Grid().NumberOfLongitudes();
   std::size_t nnode = inp_model.q().N();
   std::size_t lMax = inp_model.GSH_Grid().MaxDegree();

   // construct dxi = 0
   _vec_dxi = vvveceigvec(
       _num_layers, vveceigvec(inp_model.q().N(), veceigvec(spatialsize)));

   //    std::cout << "\nCheck 1\n";
   // now fill out dxi:
   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
      int laynum = inp_model.Node_Information().LayerNumber(idxelem);
      for (int idxnode = 0; idxnode < nnode; ++idxnode) {

         double radr =
             inp_model.Node_Information().NodeRadius(idxelem, idxnode);
         //  auto radr = node_data.NodeRadius(idxelem, idxnode);
         auto multfact = 1.0;
         auto raduse = radr;

         // check if within planet
         if (radr > inp_model.Node_Information().PlanetRadius()) {
            double planetrad = inp_model.Node_Information().PlanetRadius();
            double outerrad = inp_model.Node_Information().OuterRadius();
            raduse = planetrad;
            multfact = (outerrad - radr) / (outerrad - planetrad);
         }

         // fill out dxi
         std::size_t idxspatial = 0;
         for (auto it : inp_model.GSH_Grid().CoLatitudes()) {
            for (auto ip : inp_model.GSH_Grid().Longitudes()) {
               Eigen::Vector3cd tmpeig = Eigen::Vector3cd::Zero();
               tmpeig(1) =
                   inp_map.RadialMapping(laynum)(raduse, it, ip) * multfact;
               _vec_dxi[idxelem][idxnode][idxspatial++] = tmpeig;
            }
         }
      }
   }
   //    std::cout << "\nCheck 2\n";
   // construct dxi in spherical harmonic decomposition
   //  first step is to find the gradient of dxi
   _vec_dxilm =
       vvveceigvec(_num_layers, vveceigvec(inp_model.q().N(),
                                           veceigvec((lMax + 1) * (lMax + 1))));
   {
      auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
      auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();

      //   std::cout << "\nCheck 2.1\n";
      for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
         //   int idxpolymax = n;
         //   if (idxelem == nelem - 1) {
         //      idxpolymax = npoly + 1;
         //   }
         for (int idxpoly = 0; idxpoly < nnode; ++idxpoly) {
            std::vector<std::complex<double>> vec_dxim(spatialsize, 0.0),
                vec_dxi0(spatialsize, 0.0), vec_dxip(spatialsize, 0.0);

            for (int idxint = 0; idxint < spatialsize; ++idxint) {
               vec_dxim[idxint] = _vec_dxi[idxelem][idxpoly][idxint](0);
               vec_dxi0[idxint] = _vec_dxi[idxelem][idxpoly][idxint](1);
               vec_dxip[idxint] = _vec_dxi[idxelem][idxpoly][idxint](2);
            }

            // transform
            std::vector<std::complex<double>> vec_lm_dxim(_sizepm, 0.0),
                vec_lm_dxi0(_size0, 0.0), vec_lm_dxip(_sizepm, 0.0);

            inp_model.GSH_Grid().ForwardTransformation(lMax, -1, vec_dxim,
                                                       vec_lm_dxim);
            inp_model.GSH_Grid().ForwardTransformation(lMax, 0, vec_dxi0,
                                                       vec_lm_dxi0);
            inp_model.GSH_Grid().ForwardTransformation(lMax, 1, vec_dxip,
                                                       vec_lm_dxip);

            // std::cout << "\nCheck 2.2\n";
            // fill out "storage"
            _vec_dxilm[idxelem][idxpoly][0](1) = vec_lm_dxi0[0];
            // std::cout << "\nCheck 2.3\n";
            for (int idxl = 1; idxl < lMax + 1; ++idxl) {
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  auto idxlm = idxl * idxl + idxl + idxm;
                  _vec_dxilm[idxelem][idxpoly][idxlm](0) =
                      vec_lm_dxim[idxlm - 1];
                  _vec_dxilm[idxelem][idxpoly][idxlm](1) = vec_lm_dxi0[idxlm];
                  _vec_dxilm[idxelem][idxpoly][idxlm](2) =
                      vec_lm_dxip[idxlm - 1];
               }
            }
            // std::cout << "\nCheck 2.4\n";
         }
      }
   }
   //    std::cout << "\nCheck 3\n";
   {
      // construct da = 0
      Eigen::Matrix3cd mat_f0;
      mat_f0 = Eigen::Matrix3cd::Zero(3, 3);
      _vec_df = vvveceig(
          _num_layers, vveceig(inp_model.q().N(), veceig(spatialsize, mat_f0)));
      auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
      auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
      auto _sizepp = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 2).size();
      //   EARTHMATRIX3 vec_df(nelem * (npoly + 1),
      //                       std::vector<MATRIX3cd>(intsize, mat_0));

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
                inp_model.Node_Information().NodeRadius(idxelem, idxpoly);
            double inv2 =
                2.0 / inp_model.Node_Information().ElementWidth(idxelem);
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
            inp_model.GSH_Grid().InverseTransformation(lMax, 0, vec_ddxi0,
                                                       vec_ddxi0_spatial);

            // 0\pm
            inp_model.GSH_Grid().InverseTransformation(lMax, -1, vec_ddxim,
                                                       vec_ddxim_spatial);
            inp_model.GSH_Grid().InverseTransformation(lMax, +1, vec_ddxip,
                                                       vec_ddxip_spatial);

            //\pm 0
            inp_model.GSH_Grid().InverseTransformation(lMax, -1, vec_ddxim0,
                                                       vec_ddxim0_spatial);
            inp_model.GSH_Grid().InverseTransformation(lMax, +1, vec_ddxip0,
                                                       vec_ddxip0_spatial);

            //\pm \mp
            inp_model.GSH_Grid().InverseTransformation(lMax, 0, vec_ddxipm,
                                                       vec_ddxipm_spatial);
            inp_model.GSH_Grid().InverseTransformation(lMax, 0, vec_ddximp,
                                                       vec_ddximp_spatial);

            //\pm \pm
            inp_model.GSH_Grid().InverseTransformation(lMax, -2, vec_ddximm,
                                                       vec_ddximm_spatial);
            inp_model.GSH_Grid().InverseTransformation(lMax, +2, vec_ddxipp,
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
               for (auto it : inp_model.GSH_Grid().CoLatitudes()) {
                  for (auto ip : inp_model.GSH_Grid().Longitudes()) {
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

   //    std::cout << "\nCheck 4\n";
   // construct da = 0
   Eigen::Matrix3cd mat_a0;
   mat_a0 = Eigen::Matrix3cd::Zero(3, 3);
   _vec_da =
       vvveceig(_num_layers, vveceig(inp_model.q().N(), veceig(spatialsize)));
   {
      Eigen::Matrix3cd mat_metric = Eigen::Matrix3cd::Zero();
      mat_metric(0, 2) = -1.0;
      mat_metric(1, 1) = 1.0;
      mat_metric(2, 0) = -1.0;

      // loop through
      for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

         // finding 0-component derivative
         for (int idxpoly = 0; idxpoly < nnode; ++idxpoly) {
            for (int idxinner = 0; idxinner < spatialsize; ++idxinner) {
               Eigen::Matrix3cd tmp1 =
                   inp_model.InverseF_Point(idxelem, idxpoly, idxinner) *
                   mat_metric * _vec_df[idxelem][idxpoly][idxinner];
               _vec_da[idxelem][idxpoly][idxinner] +=
                   inp_model.LaplaceTensor_Point(idxelem, idxpoly, idxinner) *
                   (-tmp1(2, 0) + tmp1(1, 1) - tmp1(0, 2));
               Eigen::Matrix3cd tmp2 =
                   tmp1 * mat_metric *
                   inp_model.LaplaceTensor_Point(idxelem, idxpoly, idxinner);
               _vec_da[idxelem][idxpoly][idxinner] -= tmp2;
               _vec_da[idxelem][idxpoly][idxinner] -= tmp2.transpose();
            }
         }
      }
   }
};

// radial map only
MappingPerturbation::MappingPerturbation(const Density3D &inp_model,
                                         const std::string &pathtofile,
                                         const int lmin, const int lmodmax)
    : _ref_model{inp_model} {
   std::size_t _num_layers = inp_model.Num_Elements();
   std::size_t spatialsize = inp_model.GSH_Grid().NumberOfCoLatitudes() *
                             inp_model.GSH_Grid().NumberOfLongitudes();
   std::size_t nnode = inp_model.q().N();
   std::size_t lMax = inp_model.GSH_Grid().MaxDegree();

   // now this constructor deals with reading in Phobos data in particular
   // find length of file
   std::size_t num_lines = TestTools::filelength(pathtofile);
   std::size_t model_lmax = (std::sqrt(1 + 8 * num_lines) - 3) / 2;
   std::vector<int> vec_l(num_lines), vec_m(num_lines);
   std::vector<double> vec_A(num_lines), vec_SA(num_lines), vec_B(num_lines),
       vec_SB(num_lines);
   std::fstream modelfile;
   modelfile.open(pathtofile, std::ios::in);
   // getting information out of file
   if (modelfile.is_open()) {
      // test first line
      // getline(modelfile, modeltitle);
      for (int idx = 0; idx < num_lines; ++idx) {
         modelfile >> vec_l[idx] >> vec_m[idx] >> vec_A[idx] >> vec_SA[idx] >>
             vec_B[idx] >> vec_SB[idx];
         // move to next line
         modelfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }

      modelfile.close();
   } else {
      std::cout << "Couldn't open!" << "\n";
      assert("Model not found!");
   }
   // appending to vec_A if model lmax smaller than lmax:
   auto coefficientnumber =
       GSHTrans::GSHIndices<GSHTrans::NonNegative>(lMax, lMax, 0).size();
   auto coefficientnumberall =
       GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
   std::vector<std::complex<double>> vec_lm_rad(coefficientnumber, 0.0);
   //    // std::vector<std::complex<double>>
   //    vec_lm_rad_full(coefficientnumberall, 0.0);
   //    int maxmodell = std::min(lmodmax, lMax);
   {
      std::complex<double> imag1(0.0, 1.0);
      std::size_t idxoverall = 0;
      std::size_t maxl = 0, tmp_maxl = 0;
      if (lMax > model_lmax) {
         tmp_maxl = model_lmax;
      } else {
         tmp_maxl = lMax;
      }
      if (lmodmax > tmp_maxl) {
         maxl = tmp_maxl;
      } else {
         maxl = lmodmax;
      }
      // maxl = std::min(lModMax, tmp_maxl);
      for (int l = 0; l < lmin; ++l) {
         for (int m = 0; m < l + 1; ++m) {
            ++idxoverall;
         }
      }
      for (int l = lmin; l < maxl + 1; ++l) {
         for (int m = 0; m < l + 1; ++m) {
            double multfact = 0.5;
            // multfact = 1.0;
            if (m == 0) {
               multfact = 1.0;
            }
            vec_lm_rad[idxoverall] =
                multfact * (vec_A[idxoverall] - imag1 * vec_B[idxoverall]);
            ++idxoverall;
         }
      }
   }

   //    auto spatialsize = inp_model.GSH_Grid().Longitudes().size() *
   //                       inp_model.GSH_Grid().CoLatitudes().size();
   // std::cout << "Got spatialsize\n";
   //    _vec_h = vvvecdb(_num_layers, vvecdb(inp_model.q().N(),
   //    vecdb(spatialsize, 0.0)));
   std::vector<double> vec_outerradius(spatialsize, 0.0);
   // std::vector<std::complex<double>> vec_outerrad(spatialsize, 0.0);
   // perform SH transform to get physical outer radius
   // std::cout << "Pre-transformation\n";

   inp_model.GSH_Grid().InverseTransformation(lMax, 0, vec_lm_rad,
                                              vec_outerradius);
   // _grid.InverseTransformation(lMax, 0, vec_lm_rad_full, vec_outerrad);
   //    std::vector<std::complex<double>> vec_lm_check(coefficientnumberall,
   //    0.0);

   // for (int idx = 0; idx < spatialsize; ++idx) {
   //    vec_outerradius[idx] = vec_outerrad[idx].real();
   // }
   // max value
   //    double maxval =
   //        *std::max_element(vec_outerradius.begin(), vec_outerradius.end());
   //    double minval =
   //        *std::min_element(vec_outerradius.begin(), vec_outerradius.end());
   //    std::cout << "Max radius: " << maxval << ", min element: " << minval <<
   //    "\n";
   // for (auto &idx : vec_outerradius) {
   //    idx *= 1.0 / length_norm;
   // }
   //    double physicalaverageradius = vec_A[0] / std::sqrt(4.0 * pi_db);
   double referentialradius = inp_model.Node_Information().PlanetRadius();
   std::vector<double> vec_houter(spatialsize, 0.0);
   for (int idx = 0; idx < spatialsize; ++idx) {
      vec_houter[idx] = vec_outerradius[idx] / inp_model.LengthNorm();
   }
   // _grid.ForwardTransformation(lMax,0,vec_houter,vec_lm_check);
   // std::cout << vec_lm_check[0]/sqrt(3.1415926535*4.0)<<"\n";
   // std::cout << vec_outerradius[0] << " " << vec_outerradius[100] << "\n";
   //    auto _num_layers = inp_model.Node_Information().NumberOfElements();
   // construct dxi = 0
   _vec_dxi = vvveceigvec(
       _num_layers, vveceigvec(inp_model.q().N(), veceigvec(spatialsize)));

   // fill out h from mapping
   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

      int laynum = inp_model.Node_Information().LayerNumber(idxelem);

      // looping through nodes
      for (int idxnode = 0; idxnode < inp_model.q().N(); ++idxnode) {
         auto radr = inp_model.Node_Information().NodeRadius(idxelem, idxnode);
         auto multfact = 1.0;

         // check if within planetem = 0; idxelem < _num_layers; ++idxelem) {
         //       int laynum =
         //       inp_model.Node_Information().LayerNumber(idxelem); for (int
         //       idxnode = 0; idxnode < nnode; ++idxnode) {

         //          double radr =
         //              inp_model.Node_Information().NodeRadius(idxelem,
         //              idxnode);
         //          //  auto radr = node_data.NodeRadius(idxelem, idxnode);
         //          auto multfact = 1.0;
         //          auto raduse = radr;

         //          // check if within planet
         //          if (radr > inp_model.Node_Information().PlanetRadius()) {
         //             double planetrad =
         //             inp_model.Node_Information().PlanetRadius(); double
         //             outerrad = inp_model.Node_Information().OuterRadius();
         //             raduse = planetrad;
         //             multfact = (outerrad - radr) / (outerrad - planetrad);
         //          }

         //          // fill out dxi
         //          std::size_t idxspatial = 0;
         //          for (auto it : inp_model.GSH_Grid().CoLatitudes()) {
         //             for (auto ip : inp_model.GSH_Grid().Longitudes()) {
         //                Eigen::Vector3cd tmpeig = Eigen::Vector3cd::Zero();
         //                tmpeig(1) = 0.0;
         //                _vec_dxi[idxelem][idxnode][idxspatial++] = tmpeig;
         //             }
         //          }
         //       }
         //    }
         if (radr > inp_model.Node_Information().PlanetRadius()) {
            multfact = (inp_model.Node_Information().OuterRadius() - radr) /
                       (inp_model.Node_Information().OuterRadius() -
                        inp_model.Node_Information().PlanetRadius());
         } else {
            multfact = radr / inp_model.Node_Information().PlanetRadius();
         }

         // fill out h
         int idxspatial = 0;
         for (auto it : inp_model.GSH_Grid().CoLatitudes()) {
            for (auto ip : inp_model.GSH_Grid().Longitudes()) {
               Eigen::Vector3cd tmpeig = Eigen::Vector3cd::Zero();
               tmpeig(1) = vec_houter[idxspatial] * multfact;
               _vec_dxi[idxelem][idxnode][idxspatial] = tmpeig;
               ++idxspatial;
            }
         }
      }
   }

   //    std::cout << "\nCheck 1\n";
   // now fill out dxi:
   //    for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
   //       int laynum = inp_model.Node_Information().LayerNumber(idxelem);
   //       for (int idxnode = 0; idxnode < nnode; ++idxnode) {

   //          double radr =
   //              inp_model.Node_Information().NodeRadius(idxelem, idxnode);
   //          //  auto radr = node_data.NodeRadius(idxelem, idxnode);
   //          auto multfact = 1.0;
   //          auto raduse = radr;

   //          // check if within planet
   //          if (radr > inp_model.Node_Information().PlanetRadius()) {
   //             double planetrad =
   //             inp_model.Node_Information().PlanetRadius(); double outerrad =
   //             inp_model.Node_Information().OuterRadius(); raduse =
   //             planetrad; multfact = (outerrad - radr) / (outerrad -
   //             planetrad);
   //          }

   //          // fill out dxi
   //          std::size_t idxspatial = 0;
   //          for (auto it : inp_model.GSH_Grid().CoLatitudes()) {
   //             for (auto ip : inp_model.GSH_Grid().Longitudes()) {
   //                Eigen::Vector3cd tmpeig = Eigen::Vector3cd::Zero();
   //                tmpeig(1) = 0.0;
   //                _vec_dxi[idxelem][idxnode][idxspatial++] = tmpeig;
   //             }
   //          }
   //       }
   //    }
   //    std::cout << "\nCheck 2\n";
   // construct dxi in spherical harmonic decomposition
   //  first step is to find the gradient of dxi
   _vec_dxilm =
       vvveceigvec(_num_layers, vveceigvec(inp_model.q().N(),
                                           veceigvec((lMax + 1) * (lMax + 1))));
   {
      auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
      auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();

      //   std::cout << "\nCheck 2.1\n";
      for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
         //   int idxpolymax = n;
         //   if (idxelem == nelem - 1) {
         //      idxpolymax = npoly + 1;
         //   }
         for (int idxpoly = 0; idxpoly < nnode; ++idxpoly) {
            std::vector<std::complex<double>> vec_dxim(spatialsize, 0.0),
                vec_dxi0(spatialsize, 0.0), vec_dxip(spatialsize, 0.0);

            for (int idxint = 0; idxint < spatialsize; ++idxint) {
               vec_dxim[idxint] = _vec_dxi[idxelem][idxpoly][idxint](0);
               vec_dxi0[idxint] = _vec_dxi[idxelem][idxpoly][idxint](1);
               vec_dxip[idxint] = _vec_dxi[idxelem][idxpoly][idxint](2);
            }

            // transform
            std::vector<std::complex<double>> vec_lm_dxim(_sizepm, 0.0),
                vec_lm_dxi0(_size0, 0.0), vec_lm_dxip(_sizepm, 0.0);

            inp_model.GSH_Grid().ForwardTransformation(lMax, -1, vec_dxim,
                                                       vec_lm_dxim);
            inp_model.GSH_Grid().ForwardTransformation(lMax, 0, vec_dxi0,
                                                       vec_lm_dxi0);
            inp_model.GSH_Grid().ForwardTransformation(lMax, 1, vec_dxip,
                                                       vec_lm_dxip);

            // std::cout << "\nCheck 2.2\n";
            // fill out "storage"
            _vec_dxilm[idxelem][idxpoly][0](1) = vec_lm_dxi0[0];
            // std::cout << "\nCheck 2.3\n";
            for (int idxl = 1; idxl < lMax + 1; ++idxl) {
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  auto idxlm = idxl * idxl + idxl + idxm;
                  _vec_dxilm[idxelem][idxpoly][idxlm](0) =
                      vec_lm_dxim[idxlm - 1];
                  _vec_dxilm[idxelem][idxpoly][idxlm](1) = vec_lm_dxi0[idxlm];
                  _vec_dxilm[idxelem][idxpoly][idxlm](2) =
                      vec_lm_dxip[idxlm - 1];
               }
            }
            // std::cout << "\nCheck 2.4\n";
         }
      }
   }
   //    std::cout << "\nCheck 3\n";
   {
      // construct da = 0
      Eigen::Matrix3cd mat_f0;
      mat_f0 = Eigen::Matrix3cd::Zero(3, 3);
      _vec_df = vvveceig(
          _num_layers, vveceig(inp_model.q().N(), veceig(spatialsize, mat_f0)));
      auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
      auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
      auto _sizepp = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 2).size();
      //   EARTHMATRIX3 vec_df(nelem * (npoly + 1),
      //                       std::vector<MATRIX3cd>(intsize, mat_0));

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
                inp_model.Node_Information().NodeRadius(idxelem, idxpoly);
            double inv2 =
                2.0 / inp_model.Node_Information().ElementWidth(idxelem);
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
            inp_model.GSH_Grid().InverseTransformation(lMax, 0, vec_ddxi0,
                                                       vec_ddxi0_spatial);

            // 0\pm
            inp_model.GSH_Grid().InverseTransformation(lMax, -1, vec_ddxim,
                                                       vec_ddxim_spatial);
            inp_model.GSH_Grid().InverseTransformation(lMax, +1, vec_ddxip,
                                                       vec_ddxip_spatial);

            //\pm 0
            inp_model.GSH_Grid().InverseTransformation(lMax, -1, vec_ddxim0,
                                                       vec_ddxim0_spatial);
            inp_model.GSH_Grid().InverseTransformation(lMax, +1, vec_ddxip0,
                                                       vec_ddxip0_spatial);

            //\pm \mp
            inp_model.GSH_Grid().InverseTransformation(lMax, 0, vec_ddxipm,
                                                       vec_ddxipm_spatial);
            inp_model.GSH_Grid().InverseTransformation(lMax, 0, vec_ddximp,
                                                       vec_ddximp_spatial);

            //\pm \pm
            inp_model.GSH_Grid().InverseTransformation(lMax, -2, vec_ddximm,
                                                       vec_ddximm_spatial);
            inp_model.GSH_Grid().InverseTransformation(lMax, +2, vec_ddxipp,
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
               for (auto it : inp_model.GSH_Grid().CoLatitudes()) {
                  for (auto ip : inp_model.GSH_Grid().Longitudes()) {
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

   //    std::cout << "\nCheck 4\n";
   // construct da = 0
   Eigen::Matrix3cd mat_a0;
   mat_a0 = Eigen::Matrix3cd::Zero(3, 3);
   _vec_da =
       vvveceig(_num_layers, vveceig(inp_model.q().N(), veceig(spatialsize)));
   {
      Eigen::Matrix3cd mat_metric = Eigen::Matrix3cd::Zero();
      mat_metric(0, 2) = -1.0;
      mat_metric(1, 1) = 1.0;
      mat_metric(2, 0) = -1.0;

      // loop through
      for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

         // finding 0-component derivative
         for (int idxpoly = 0; idxpoly < nnode; ++idxpoly) {
            for (int idxinner = 0; idxinner < spatialsize; ++idxinner) {
               Eigen::Matrix3cd tmp1 =
                   inp_model.InverseF_Point(idxelem, idxpoly, idxinner) *
                   mat_metric * _vec_df[idxelem][idxpoly][idxinner];
               _vec_da[idxelem][idxpoly][idxinner] +=
                   inp_model.LaplaceTensor_Point(idxelem, idxpoly, idxinner) *
                   (-tmp1(2, 0) + tmp1(1, 1) - tmp1(0, 2));
               Eigen::Matrix3cd tmp2 =
                   tmp1 * mat_metric *
                   inp_model.LaplaceTensor_Point(idxelem, idxpoly, idxinner);
               _vec_da[idxelem][idxpoly][idxinner] -= tmp2;
               _vec_da[idxelem][idxpoly][idxinner] -= tmp2.transpose();
            }
         }
      }
   }
};

}   // namespace GeneralEarthModels

#endif