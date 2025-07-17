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

namespace TestTools {

using namespace GravityFunctions;

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

   auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
   auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
   auto rphys = [&vec_noderadii](int idxelem, double x) {
      return ((vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * x +
              (vec_noderadii[idxelem + 1] + vec_noderadii[idxelem])) *
             0.5;
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

         auto rval = GravityFunctions::StandardIntervalMap(
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

   auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
   auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
   auto rphys = [&vec_noderadii](int idxelem, double x) {
      return ((vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * x +
              (vec_noderadii[idxelem + 1] + vec_noderadii[idxelem])) *
             0.5;
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
         auto rval = GravityFunctions::StandardIntervalMap(
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

auto
filelength(const std::string &filename) {
   std::ifstream inputFile;
   std::string line;
   int filelength = 0;
   inputFile.open(filename, std::ios::in);
   if (!inputFile.is_open()) {
      std::cout << "Failed to open file\n";
   } else {
      while (std::getline(inputFile, line)) {
         ++filelength;
      }
      inputFile.close();
   }
   return filelength;
};

template <typename FLOAT = double, typename INTEGRAL = int> class EarthModel {
 public:
   using size_type = INTEGRAL;
   using value_type = FLOAT;
   // using InterpA = Interpolation::Akima<std::vector<double>::iterator,
   //                                      std::vector<double>::iterator>;
   // template <typename NFLOAT = FLOAT>
   // using InterpA = Interpolation::CubicSpline<std::vector<double>::iterator,
   //                                            std::vector<double>::iterator>;
   using myvector = std::vector<FLOAT>;
   using myiter = myvector::iterator;
   using InterpA = Interpolation::CubicSpline<myiter, myiter>;
   // using InterpA = Interpolation::CubicSpline<std::vector<FLOAT>::iterator,
   //                                            std::vector<FLOAT>::iterator>;
   // constructors
   EarthModel() {};
   EarthModel(const std::string &);
   template <template <typename> class ParameterModel>
   EarthModel(const std::string &, const ParameterModel<FLOAT> &);

   // norms
   FLOAT LengthNorm() const { return length_norm; };
   FLOAT MassNorm() const { return mass_norm; };
   FLOAT TimeNorm() const { return time_norm; }
   FLOAT DensityNorm() const { return density_norm; };
   FLOAT InertiaNorm() const { return inertia_norm; };
   FLOAT VelocityNorm() const { return velocity_norm; };
   FLOAT AccelerationNorm() const { return acceleration_norm; };
   FLOAT ForceNorm() const { return force_norm; };
   FLOAT StressNorm() const { return stress_norm; };
   FLOAT GravitationalConstant() const { return gravitational_constant; };

   // Geometry of model
   int NumberOfLayers() const { return _numlayers; };
   int LayerLowerIndex(int i) const { return _vec_indices[i][0]; };
   int LayerUpperIndex(int i) const { return _vec_indices[i][1]; };
   int LayerIndexDifference(int i) const {
      return _vec_indices[i][1] - _vec_indices[i][0];
   };
   auto LayerRadii() const { return layered_radii; };
   auto LayerRadii(int i) const { return layered_radii[i]; };
   auto LayerRadiiNumber(int i) const { return layered_radii[i].size(); };

   FLOAT LowerRadius(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside the number of layers in the model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return vec_layers[i];
   }
   FLOAT UpperRadius(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside the number of layers in the model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return vec_layers[i + 1];
   }
   FLOAT OuterRadius() const { return vec_layers[_numlayers]; }

   // Isotropy/fluid/solid etc
   bool IsIsotropic() const { return _isisotropic; };

   // Solid or fluid
   bool IsSolid(INTEGRAL i) const { return _issolid[i]; }
   bool IsFluid(INTEGRAL i) const { return !IsSolid(i); }

   // Density
   InterpA Density(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_rho[i];
      // return func_rhoc[i];
   };

   // Velocities
   // InterpA VP(INTEGRAL i) const { return func_vpv[i]; };
   InterpA VPV(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };

      return func_vpv[i];
   };
   InterpA VPH(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_vph[i];
   };
   // InterpA VS(INTEGRAL i) const { return vec_s_velocity[i]; };
   InterpA VSV(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_vsv[i];
   };
   InterpA VSH(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_vsh[i];
   };

   ///////////////////////////////////////////////////////////////
   ////////////////// !!!!!!!!!!!!!!!!!!!!!!!!!!! ////////////////
   ///////////////////////////////////////////////////////////////
   InterpA VS(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_vsv[i];
   };
   InterpA VP(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_vpv[i];
   };
   ///////////////////////////////////////////////////////////////
   ////////////////// !!!!!!!!!!!!!!!!!!!!!!!!!!! ////////////////
   ///////////////////////////////////////////////////////////////

   // Returning eta:
   auto Eta(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_eta[i];
   }

   // returning A, C, N, L, kappa, mu
   auto A(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return Density(i)(x) * VPH(i)(x) * VPH(i)(x);
      };
      return aret;
   };
   auto C(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return Density(i)(x) * VPV(i)(x) * VPV(i)(x);
      };
      return aret;
   };
   auto N(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return Density(i)(x) * VSH(i)(x) * VSH(i)(x);
      };
      return aret;
   };
   auto L(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return Density(i)(x) * VSV(i)(x) * VSV(i)(x);
      };
      return aret;
   };
   auto F(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return Eta(i)(x) * (A(i)(x) - 2 * L(i)(x));
      };
      return aret;
   };
   auto Kappa(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return (C(i)(x) + 4.0 * (A(i)(x) - N(i)(x) + F(i)(x))) / 9.0;
      };
      return aret;
   };
   auto Mu(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return (C(i)(x) + A(i)(x) + 6.0 * L(i)(x) + 5.0 * N(i)(x) -
                 2.0 * F(i)(x)) /
                15.0;
      };
      return aret;
   };

 private:
   using vecdb = std::vector<double>;
   using vvecdb = std::vector<vecdb>;
   vecdb vec_radius, vec_rho, vec_vpv, vec_vsv, vec_qkappa, vec_qshear, vec_vph,
       vec_vsh, vec_eta, vec_layers;
   vvecdb layered_radii = vvecdb(1, vecdb());
   vvecdb layered_rho = vvecdb(1, vecdb());
   vvecdb layered_vpv = vvecdb(1, vecdb());
   vvecdb layered_vsv = vvecdb(1, vecdb());
   vvecdb layered_qkappa = vvecdb(1, vecdb());
   vvecdb layered_qshear = vvecdb(1, vecdb());
   vvecdb layered_vph = vvecdb(1, vecdb());
   vvecdb layered_vsh = vvecdb(1, vecdb());
   vvecdb layered_eta = vvecdb(1, vecdb());

   bool _isisotropic = true;
   std::vector<bool> _issolid = std::vector<bool>(1, true);
   std::vector<std::vector<int>> _vec_indices;
   // density
   // Interpolation::CubicSpline<std::vector<double>::iterator,
   //                            std::vector<double>::iterator>
   //     checkval;
   std::vector<InterpA> func_rho, func_vpv, func_vsv, func_qkappa, func_qshear,
       func_vph, func_vsh, func_eta;
   // std::vector<InterpC> func_rhoc;
   // vvecdb layered_radii = vvecdb(1, vecdb());
   // std::vector<std::vector<double>> layered_vpv, layered_vsv,
   //     layered_qkappa, layered_qshear, layered_vph, layered_vsh, layered_eta;

   // model information
   std::string modeltitle;
   int ifanis, ifdeck, numnodes, nic, noc, _numlayers;
   double tref;

   FLOAT length_norm, mass_norm, time_norm, density_norm, inertia_norm,
       velocity_norm, acceleration_norm, force_norm, stress_norm,
       gravitational_constant;

   ///////////////////////////////////////////////
   /////////////// ??????????????? ///////////////
   ///////////////////////////////////////////////

   // find layers of model
   std::vector<double> findlayers(const std::vector<double> &vec_sorted) {
      std::vector<double> vec_bounds;
      auto i1 = vec_sorted.begin();
      while (i1 != vec_sorted.end()) {
         vec_bounds.push_back(*i1);
         i1 = std::adjacent_find(++i1, vec_sorted.end());
      }
      vec_bounds.push_back(*(--i1));
      return vec_bounds;
   };

   // find indices
   std::vector<std::size_t>
   layerindices(const std::vector<double> &vec_sorted) {
      std::vector<std::size_t> vec_indices;
      auto i1 = vec_sorted.begin();
      while (i1 != vec_sorted.end()) {
         vec_indices.push_back(std::distance(vec_sorted.begin(), i1));
         i1 = std::adjacent_find(++i1, vec_sorted.end());
      }
      vec_indices.push_back(
          std::distance(vec_sorted.begin(), vec_sorted.end()) - 1);
      return vec_indices;
   }
};

// "default" constructor
template <typename FLOAT, typename INTEGRAL>
EarthModel<FLOAT, INTEGRAL>::EarthModel(const std::string &pathtofile)
    : EarthModel(pathtofile, EarthModels::EarthConstants<FLOAT>()){};

// full constructor
template <typename FLOAT, typename INTEGRAL>
template <template <typename> class ParameterModel>
EarthModel<FLOAT, INTEGRAL>::EarthModel(
    const std::string &pathtofile, const ParameterModel<FLOAT> &ModelConstants)
    : length_norm(ModelConstants.LengthNorm()),
      mass_norm(ModelConstants.MassNorm()),
      time_norm(ModelConstants.TimeNorm()),
      density_norm(ModelConstants.DensityNorm()),
      velocity_norm(ModelConstants.VelocityNorm()),
      acceleration_norm(ModelConstants.AccelerationNorm()),
      force_norm(ModelConstants.ForceNorm()),
      stress_norm(ModelConstants.StressNorm()),
      inertia_norm(ModelConstants.InertiaNorm()),
      gravitational_constant(ModelConstants.GravitationalConstant()) {

   // std::cout << this->density_norm << "\n";
   // opening file
   std::fstream modelfile;
   modelfile.open(pathtofile, std::ios::in);

   // getting information out of file
   if (modelfile.is_open()) {
      // get first line (title)
      getline(modelfile, modeltitle);

      // extract information from second line and move to next line
      modelfile >> ifanis >> tref >> ifdeck;
      modelfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      // extract information from third line and move to next line
      modelfile >> numnodes >> nic >> noc;
      modelfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      // loop through the deck
      int laynum = 0;
      int idxouter = 0;
      int idxinner = 0;

      while (idxouter < numnodes) {

         std::vector<double> tmp_radius;
         bool samelayer = true;

         // while (samelayer) {
         // double radius;
         double radius, rho, vpv, vsv, qkappa, qshear, vph, vsh, eta;
         modelfile >> radius >> rho >> vpv >> vsv >> qkappa >> qshear >> vph >>
             vsh >> eta;
         if (idxinner > 0 && (radius / this->length_norm ==
                              layered_radii[laynum][idxinner - 1])) {
            // move to next layer
            layered_radii.push_back({radius / this->length_norm});
            layered_rho.push_back({rho / this->density_norm});
            layered_vpv.push_back({vpv / this->velocity_norm});
            layered_vsv.push_back({vsv / this->velocity_norm});
            layered_qkappa.push_back({qkappa});
            layered_qshear.push_back({qshear});
            layered_vph.push_back({vph / this->velocity_norm});
            layered_vsh.push_back({vsh / this->velocity_norm});
            layered_eta.push_back({eta});

            // isotropy
            if (_isisotropic && vpv != vsv) {
               _isisotropic = false;
            }

            // fluid/solid:
            if (vsv == 0.0 && vsh == 0) {
               _issolid.push_back(false);
            } else {
               _issolid.push_back(true);
            }

            // set idxinner back to zero
            idxinner = 0;
            ++laynum;
         } else {
            // put next value in current layer in
            layered_radii[laynum].push_back(radius / this->length_norm);
            layered_rho[laynum].push_back(rho / this->density_norm);
            layered_vpv[laynum].push_back(vpv / this->velocity_norm);
            layered_vsv[laynum].push_back(vsv / this->velocity_norm);
            layered_qkappa[laynum].push_back(qkappa);
            layered_qshear[laynum].push_back(qshear);
            layered_vph[laynum].push_back(vph / this->velocity_norm);
            layered_vsh[laynum].push_back(vsh / this->velocity_norm);
            layered_eta[laynum].push_back(eta);

            // check whether solid or fluid
            if (idxouter == 0) {
               if (vsv == 0.0 && vsh == 0.0) {
                  _issolid[laynum] = false;
               }
            }
         }

         // move to next line
         modelfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

         // increment counters
         ++idxinner;
         ++idxouter;
         // }
      }
      // for (auto &idxouter : layered_radii) {
      //    // std::cout << "HELLO\n";
      //    std::cout << idxouter.front() << " " << idxouter.back() << "\n";
      // }
      // std::cout << "Seg test 1\n";
      _numlayers = layered_radii.size();
      vec_layers.reserve(_numlayers + 1);
      vec_layers.push_back(0.0);
      // std::generate(
      //     vec_layers.begin() + 1, vec_layers.end(),
      //     [n = 0, this]() mutable { return layered_radii[n++].back(); });
      // std::cout << "Seg test 2\n";
      _vec_indices =
          std::vector<std::vector<int>>(_numlayers, std::vector<int>(2, 0));
      _vec_indices[0][0] = 0;
      // std::cout << "Seg test 3\n";
      for (int idx = 0; idx < _numlayers; ++idx) {
         vec_layers.push_back(layered_radii[idx].back());
         if (idx != 0) {
            _vec_indices[idx][0] = _vec_indices[idx - 1][1] + 1;
         }
         _vec_indices[idx][1] =
             _vec_indices[idx][0] + layered_radii[idx].size() - 1;
         // std::cout << "Seg test " << idx + 3 << "\n";
      }
      // std::cout << "Size: " << vec_layers.size() << "\n";
      // for (auto &idx : vec_layers) {
      //    std::cout << idx << "\n";
      // }

      for (int idx = 0; idx < _numlayers; ++idx) {
         // iterators to start and end of layer of radius
         auto it1 = layered_radii[idx].begin();
         auto it2 = layered_radii[idx].end();

         // iterators to beginning of this layer for all data
         auto it_rho = layered_rho[idx].begin();
         auto it_vpv = layered_vpv[idx].begin();
         auto it_vsv = layered_vsv[idx].begin();
         auto it_qkappa = layered_qkappa[idx].begin();
         auto it_qshear = layered_qshear[idx].begin();
         auto it_vph = layered_vph[idx].begin();
         auto it_vsh = layered_vsh[idx].begin();
         auto it_eta = layered_eta[idx].begin();

         // pushback
         func_rho.push_back(InterpA(it1, it2, it_rho));
         func_vpv.push_back(InterpA(it1, it2, it_vpv));
         func_vsv.push_back(InterpA(it1, it2, it_vsv));
         func_qkappa.push_back(InterpA(it1, it2, it_qkappa));
         func_qshear.push_back(InterpA(it1, it2, it_qshear));
         func_vph.push_back(InterpA(it1, it2, it_vph));
         func_vsh.push_back(InterpA(it1, it2, it_vsh));
         func_eta.push_back(InterpA(it1, it2, it_eta));
      }

      // for (int idx = 0; idx < N; ++idx) {
      //    double radius, rho, vpv, vsv, qkappa, qshear, vph, vsh, eta;
      //    modelfile >> radius >> rho >> vpv >> vsv >> qkappa >> qshear >> vph
      //    >>
      //        vsh >> eta;
      //    vec_radius.push_back(radius);
      //    // vec_rho.push_back(rho);
      //    // vec_vpv.push_back(vpv);
      //    // vec_vsv.push_back(vsv);
      //    // vec_qkappa.push_back(qkappa);
      //    // vec_qshear.push_back(qshear);
      //    // vec_vph.push_back(vph);
      //    // vec_vsh.push_back(vsh);
      //    // vec_eta.push_back(eta);

      //    // move to next line
      //    modelfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      // }

      modelfile.close();
   } else {
      assert("Model not found!");
   }

   // finding layers
   // this->vec_layers = this->findlayers(this->vec_radius);

   // auto vec_indices = this->layerindices(this->vec_radius);
   // for (auto &idx : vec_indices) {
   //    std::cout << idx << "\n";
   // }
   // {
   //    int idxouter = 0;
   //    for (int idxlayers = 0; idxlayers < _numlayers; ++idxlayers) {
   //       std::vector<double> tmp;
   //       while (vec_radius[idxouter] != vec_layers[idxlayers + 1]) {
   //          tmp.push_back(vec_radius[idxouter]);
   //       }
   //    }
   // }
};



// full constructor

void
PhobosRead(const std::string &pathtofile) {

   // std::cout << this->density_norm << "\n";
   // opening file

   // find length of file
   std::size_t num_lines = filelength(pathtofile);

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