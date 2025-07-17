#ifndef GRAVITY_SPHERICAL_INTEGRATOR_H
#define GRAVITY_SPHERICAL_INTEGRATOR_H

#include <Eigen/Core>
#include <Eigen/Dense>
// #include <Eigen/CholmodSupport>
#include <Eigen/Sparse>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <PlanetaryModel/All>
// #include <TomographyModels/All>
#include "SphericalGeometryPreconditioner.h"
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

namespace GravityFunctions {
// return physical radius for mappings between [-1,1] and [x1,x2]
template <typename FLOAT>
FLOAT
StandardIntervalMap(const FLOAT &x, const FLOAT &x1, const FLOAT &x2) {
   return ((x2 - x1) * x + (x1 + x2)) * 0.5;
};

// return scaled x for if needed
template <typename FLOAT>
FLOAT
ScaledIntervalMap(const FLOAT &x, const FLOAT &x1, const FLOAT &x2) {
   return x + (x1 + x2) / (x2 - x1);
};

// returns depth from radius or vice versa
template <typename FLOAT>
std::vector<FLOAT>
DepthRadiusSwap(const std::vector<FLOAT> &vec_init) {
   std::vector<FLOAT> vec_final(vec_init.size());
   for (int idx = 0; idx < vec_init.size(); ++idx) {
      vec_final[idx] = 1.0 - vec_init[idx];
   }
   return vec_final;
};

// rescale:
template <typename FLOAT>
void
RescaleVector(std::vector<FLOAT> &vec_init, const FLOAT &rescaleparameter) {
   // std::vector<FLOAT> vec_final = vec_init;
   for (int idx = 0; idx < vec_init.size(); ++idx) {
      vec_init[idx] *= rescaleparameter;
   }
};

// determines radii of nodes, if certain number of element nodes specified for
// each layer
template <typename FLOAT, template <typename, typename> class sphericalmodel>
   requires PlanetaryModel::SphericalElasticModel<sphericalmodel<FLOAT, int>>
std::vector<FLOAT>
Layer_Node(const sphericalmodel<FLOAT, int> &mymodel,
           const std::vector<int> &nnodes) {
   // total number of nodes
   int numnodes = std::reduce(nnodes.cbegin(), nnodes.cend());
   numnodes -= mymodel.NumberOfLayers() - 1;
   std::vector<FLOAT> vec_noderadii(numnodes);   // node radii

   vec_noderadii[0] = 0.0;
   int idxnode = 0;   // counter
   for (int idx = 0; idx < mymodel.NumberOfLayers(); ++idx) {
      FLOAT dx = (mymodel.UpperRadius(idx) - mymodel.LowerRadius(idx)) /
                 (static_cast<double>(nnodes[idx] - 1) * mymodel.LengthNorm());
      FLOAT radius_scaled = mymodel.LowerRadius(idx) / mymodel.LengthNorm();
      for (int idx2 = 1; idx2 < nnodes[idx]; ++idx2) {
         idxnode += 1;
         vec_noderadii[idxnode] =
             radius_scaled + dx * static_cast<double>(idx2);
      }
   };
   return vec_noderadii;
};

// determines radii of edges of elements, with a maximum element size specified
template <typename FLOAT, template <typename, typename> class sphericalmodel>
   requires PlanetaryModel::SphericalElasticModel<sphericalmodel<FLOAT, int>>
std::vector<FLOAT>
Radial_Node(const sphericalmodel<FLOAT, int> &mymodel,
            const FLOAT max_element_size) {

   // IMPORTANT: MAX ELEMENT SIZE MUST BE IN SAME UNITS AS MODEL//

   // declare vec_noderadii
   std::vector<FLOAT> vec_noderadii;   // node radii

   vec_noderadii.push_back(0.0);

   // we want to have normalised radii. Thus define normalisation:
   // FLOAT mynorm = 1.0 / mymodel.OuterRadius();
   // FLOAT mynorm = 1.0 / mymodel.LengthNorm();
   // FLOAT normmaxsize = mynorm * max_element_size;

   for (int idx = 0; idx < mymodel.NumberOfLayers(); ++idx) {
      // Determining the spacing of nodes within each layer
      FLOAT laydepth = mymodel.UpperRadius(idx) - mymodel.LowerRadius(idx);
      int numelements = std::ceil(laydepth / max_element_size);
      FLOAT nodespacing = laydepth / static_cast<FLOAT>(numelements);

      // filling out the node radii vector
      for (int idxelem = 0; idxelem < numelements; ++idxelem) {
         if (idxelem < numelements - 1) {
            FLOAT currentval = vec_noderadii.back();
            vec_noderadii.push_back(currentval + nodespacing);
         } else {
            vec_noderadii.push_back(mymodel.UpperRadius(idx));
         }
      }
   };
   return vec_noderadii;
};

// determines radii of edges of elements, with a maximum element size specified
template <typename FLOAT, template <typename, typename> class sphericalmodel>
   requires PlanetaryModel::SphericalElasticModel<sphericalmodel<FLOAT, int>>
std::vector<FLOAT>
Radial_Node(const sphericalmodel<FLOAT, int> &mymodel,
            const FLOAT max_element_size, const FLOAT ballrad) {

   // IMPORTANT: MAX ELEMENT SIZE MUST BE IN SAME UNITS AS MODEL//

   // declare vec_noderadii
   std::vector<FLOAT> vec_noderadii;   // node radii

   vec_noderadii.push_back(0.0);

   // we want to have normalised radii. Thus define normalisation:
   // FLOAT mynorm = 1.0 / mymodel.LengthNorm();
   FLOAT normmaxsize = max_element_size;

   for (int idx = 0; idx < mymodel.NumberOfLayers(); ++idx) {
      // Determining the spacing of nodes within each layer
      FLOAT laydepth = mymodel.UpperRadius(idx) - mymodel.LowerRadius(idx);
      int numelements = std::ceil(laydepth / normmaxsize);
      FLOAT nodespacing = laydepth / static_cast<FLOAT>(numelements);

      // filling out the node radii vector
      for (int idxelem = 0; idxelem < numelements; ++idxelem) {
         if (idxelem < numelements - 1) {
            FLOAT currentval = vec_noderadii.back();
            vec_noderadii.push_back(currentval + nodespacing);
         } else {
            vec_noderadii.push_back(mymodel.UpperRadius(idx));
         }
      }
   };

   // adding in the outer elements between the edge of the planet and the ball
   // FLOAT multfact = 1.5;
   FLOAT ballrad2 = ballrad;
   // std::cout << "ball radius: " << ballrad2 << std::endl;
   {
      // get size and number of elements
      FLOAT laydepth = ballrad2 - mymodel.OuterRadius();
      int numelem = std::ceil(laydepth / normmaxsize);
      FLOAT nodespacing = laydepth / static_cast<FLOAT>(numelem);

      // finish filling out vector
      for (int idxelem = 0; idxelem < numelem; ++idxelem) {
         if (idxelem < numelem - 1) {
            FLOAT currentval = vec_noderadii.back();
            vec_noderadii.push_back(currentval + nodespacing);
         } else {
            vec_noderadii.push_back(ballrad2);
         }
      }
   }
   return vec_noderadii;
};

// determines radii of all nodes within the decomposition
template <typename FLOAT>
std::vector<FLOAT>
All_Node(const std::vector<FLOAT> &vec_noderadii,
         GaussQuad::Quadrature1D<FLOAT> &q) {

   // IMPORTANT: MAX ELEMENT SIZE MUST BE IN SAME UNITS AS MODEL//

   // sizes of various vectors:
   std::size_t nelem = vec_noderadii.size() - 1;
   int npoly = q.N() - 1;

   std::size_t matlen = nelem * npoly + 1;    // size of vector
   std::vector<FLOAT> vec_allradii(matlen);   // declare vector
   for (int idx = 0; idx < nelem; ++idx) {
      for (int idx2 = 0; idx2 < npoly; ++idx2) {
         vec_allradii[idx * npoly + idx2] = StandardIntervalMap(
             q.X(idx2), vec_noderadii[idx], vec_noderadii[idx + 1]);
      }
   }
   vec_allradii[0] = 0.0;
   vec_allradii[matlen - 1] = vec_noderadii[nelem];
   return vec_allradii;
};

template <typename FLOAT, template <typename, typename> class sphericalmodel>
   requires PlanetaryModel::SphericalElasticModel<sphericalmodel<FLOAT, int>>
std::vector<std::vector<FLOAT>>
PotentialSolver1D(const sphericalmodel<FLOAT, int> &mymodel,
                  GaussQuad::Quadrature1D<FLOAT> &q,
                  const std::vector<FLOAT> &vec_noderadii) {
   std::size_t nelem = vec_noderadii.size() - 1;
   int npoly = q.N() - 1;

   double phi_0;
   phi_0 = 0.0;
   int laynum = 0;
   double rad_scale = mymodel.LengthNorm();
   // double scdiff = rad_scale / mymodel.OuterRadius();

   auto rphys = [&vec_noderadii](int idxelem, double x) {
      return ((vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * x +
              (vec_noderadii[idxelem + 1] + vec_noderadii[idxelem])) *
             0.5;
   };
   const double bigg_db = 6.6743 * std::pow(10.0, -11.0);
   const double pi_db = 3.1415926535;

   // finding the phi_0 part of the potential
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // determining the layer number
      if (laynum < mymodel.NumberOfLayers()) {
         if (!(vec_noderadii[idxelem] < mymodel.UpperRadius(laynum))) {
            laynum += 1;
         };
      }

      auto funcdens = [&mymodel, &rphys, laynum, idxelem](double x) {
         return rphys(idxelem, x) * mymodel.Density(laynum)(rphys(idxelem, x));
         // return rphys(idxelem, x) * 1000.0;
      };
      phi_0 -= 2.0 * pi_db * bigg_db *
               (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) *
               std::pow(rad_scale, 2.0) * q.Integrate(funcdens);
   }

   // initialising the potential
   std::vector<double> vec_exactprempotential(nelem + 1, 0.0),
       vec_exactpremgravity(nelem + 1, 0.0);
   vec_exactprempotential[0] += phi_0;

   std::vector<double> vec_massr(nelem + 1, 0.0);
   laynum = 0;
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {

      if (laynum < mymodel.NumberOfLayers()) {
         if (!(vec_noderadii[idxelem] < mymodel.UpperRadius(laynum))) {
            laynum += 1;
         };
      }

      auto funcdens = [&mymodel, &rphys, laynum, idxelem,
                       &vec_noderadii](double x) {
         // return rphys(idxelem, x) * myprem.Density(laynum)(rphys(idxelem,
         // x) * scdiff) * (1 - rphys(idxelem, x) / vec_noderadii[idxelem]);
         return mymodel.Density(laynum)(rphys(idxelem, x)) * rphys(idxelem, x) *
                rphys(idxelem, x);
         // return 1000.0 * rphys(idxelem - 1, x) * rphys(idxelem - 1, x);
      };
      vec_massr[idxelem + 1] =
          vec_massr[idxelem] +
          2.0 * pi_db * (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) *
              q.Integrate(funcdens) * std::pow(rad_scale, 3.0);
   };

   std::vector<double> vec_rhorint(nelem + 1, 0.0);
   laynum = 0;
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      if (laynum < mymodel.NumberOfLayers()) {
         if (!(vec_noderadii[idxelem] < mymodel.UpperRadius(laynum))) {
            laynum += 1;
         };
      }

      auto funcdens = [&mymodel, &rphys, laynum, idxelem,
                       &vec_noderadii](double x) {
         return rphys(idxelem, x) * mymodel.Density(laynum)(rphys(idxelem, x));
         // return 1000.0 * rphys(idxelem - 1, x);
      };
      vec_rhorint[idxelem + 1] =
          vec_rhorint[idxelem] +
          2.0 * pi_db * bigg_db *
              (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) *
              q.Integrate(funcdens) * std::pow(rad_scale, 2.0);
   };

   // auto stop = high_resolution_clock::now();
   // auto duration = duration_cast<microseconds>(stop - start);
   // std::cout << "Time taken by exact integral solution: "
   //           << duration.count() / 1000000.0 << " seconds" << std::endl;

   ////////////////////////////////////////////////////
   ////////////////////////////////////////////////////
   ////////////////////////////////////////////////////
   ////////////////////////////////////////////////////
   // output to files

   for (int idxelem = 1; idxelem < nelem + 1; ++idxelem) {
      vec_exactprempotential[idxelem] +=
          phi_0 -
          bigg_db * vec_massr[idxelem] /
              (vec_noderadii[idxelem] * mymodel.LengthNorm()) +
          vec_rhorint[idxelem];
      vec_exactprempotential[idxelem] *= mymodel.DensityNorm();
      vec_exactpremgravity[idxelem] =
          bigg_db * vec_massr[idxelem] /
          (std::pow(vec_noderadii[idxelem] * mymodel.LengthNorm(), 2.0));
      vec_exactpremgravity[idxelem] *= mymodel.DensityNorm();
   };
   vec_exactprempotential[0] *= mymodel.DensityNorm();
   vec_exactpremgravity[0] *= mymodel.DensityNorm();

   std::vector<std::vector<FLOAT>> vec_return;
   vec_return.push_back(vec_exactprempotential);
   vec_return.push_back(vec_exactpremgravity);
   return vec_return;
}

template <typename FLOAT>
auto
PotentialIntegrator3D(
    const std::vector<FFTWpp::vector<std::complex<double>>> &vec_denslm,
    GaussQuad::Quadrature1D<FLOAT> &q, const std::vector<FLOAT> &vec_noderadii,
    const std::vector<FLOAT> &vec_allradii, const int &lMax) {
   std::size_t nelem = vec_noderadii.size() - 1;
   int npoly = q.N() - 1;

   double phi_0;
   phi_0 = 0.0;
   // int laynum = 0;
   // double rad_scale = mymodel.LengthNorm();
   // double scdiff = rad_scale / mymodel.OuterRadius();
   // vector that will contain results
   std::vector<std::vector<std::complex<FLOAT>>> vec_output(
       vec_noderadii.size(),
       std::vector<std::complex<FLOAT>>((lMax + 1) * (lMax + 1), 0.0));

   auto rphys = [&vec_noderadii](int idxelem, double x) {
      return ((vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * x +
              (vec_noderadii[idxelem + 1] + vec_noderadii[idxelem])) *
             0.5;
   };
   auto rscale = [&vec_noderadii](int idxelem, double x) {
      return ((1.0 - vec_noderadii[idxelem] / vec_noderadii[idxelem + 1]) * x +
              (1.0 + vec_noderadii[idxelem] / vec_noderadii[idxelem + 1])) *
             0.5;
   };
   const double bigg_db = 6.6743 * std::pow(10.0, -11.0);
   const double pi_db = 3.1415926535;

   // std::vector<std::complex<double>> vec_g(nelem + 1, 0.0);
   for (int idxl = 0; idxl < lMax + 1; ++idxl) {
      for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
         std::vector<std::complex<double>> vec_f(nelem + 1, 0.0),
             vec_g(nelem + 1, 0.0);
         auto idxlm = idxl * idxl + idxl + idxm;
         // relabelling g as r^{-(l+1)} \int_{0}^{r} dr' r'^{2+l} \rho_{lm}(r')
         // use Gaussian quadrature to find the integral from
         // vec_noderadii[idxelem-1] to vec_noderadii[idxelem]
         for (int idxelem = 0; idxelem < nelem; ++idxelem) {
            auto myratio = vec_noderadii[idxelem] / vec_noderadii[idxelem + 1];
            auto rscaleint = [&myratio](double x) {
               return ((1.0 - myratio) * x + (1.0 + myratio)) * 0.5;
            };
            for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
               vec_f[idxelem + 1] +=
                   q.W(idxpoly) *
                   vec_denslm[idxelem * (npoly + 1) + idxpoly][idxlm] *
                   std::pow(rscaleint(q.X(idxpoly)),
                            2.0 + static_cast<double>(idxl));
            }
            auto intfact = 0.5 * (1.0 - myratio);
            vec_f[idxelem + 1] *= intfact;
            vec_f[idxelem + 1] *=
                vec_noderadii[idxelem + 1] * vec_noderadii[idxelem + 1];

            vec_f[idxelem + 1] +=
                vec_f[idxelem] *
                std::pow(myratio, 1.0 + static_cast<double>(idxl));
         }
         // in this for loop the lower bound is 1 as at 0 it is zero except for
         // l=0
         //  auto elemlowbound ;
         //  if(idxl ==0){elemlowbound=-1;} else{elemlowbound=0;}
         for (int idxelem = nelem - 1; idxelem > 0; --idxelem) {

            auto myratio = vec_noderadii[idxelem] / vec_noderadii[idxelem + 1];
            auto myratio2 = 1.0 / myratio;
            auto rscaleint = [&myratio2](double x) {
               return ((myratio2 - 1.0) * x + (myratio2 + 1.0)) * 0.5;
            };

            for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
               vec_g[idxelem] +=
                   q.W(idxpoly) *
                   vec_denslm[idxelem * (npoly + 1) + idxpoly][idxlm] *
                   std::pow(rscaleint(q.X(idxpoly)),
                            1.0 - static_cast<double>(idxl));
            }
            auto intfact = 0.5 * (myratio2 - 1.0);
            vec_g[idxelem] *= intfact;
            vec_g[idxelem] *= vec_noderadii[idxelem] * vec_noderadii[idxelem];
            vec_g[idxelem] += vec_g[idxelem + 1] *
                              std::pow(myratio, static_cast<double>(idxl));
         }
         if (idxl == 0) {
            auto rscaleint = [&vec_noderadii](double x) {
               return ((vec_noderadii[1] - vec_noderadii[0]) * x +
                       (vec_noderadii[1] + vec_noderadii[0])) *
                      0.5;
            };
            for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
               vec_g[0] += q.W(idxpoly) * vec_denslm[idxpoly][0] *
                           rscaleint(q.X(idxpoly));
            }
            auto intfact = 0.5 * (vec_noderadii[1] - vec_noderadii[0]);
            vec_g[0] *= intfact;
            vec_g[0] += vec_g[1];
         }

         // store result
         for (int idxelem = 0; idxelem < nelem + 1; ++idxelem) {
            auto idxuse = idxl * idxl + idxl + idxm;
            vec_output[idxelem][idxuse] = vec_f[idxelem] + vec_g[idxelem];
            vec_output[idxelem][idxuse] *=
                -4.0 * pi_db * bigg_db /
                (2.0 * static_cast<double>(idxl) + 1.0);
         }
      }
   }

   return vec_output;
}

// enum class for how the coefficients are order
enum class CoefficientOrdering {
   RadialClumped,
   SphericalDecomp

   // RadialClumped means at radius r, list (0,0), (1,-1), (1,0), (1,1) ...
   // SphericalDecomp means list (0,0)(0), (0,0)(r_1), ..., (0,0)(surface),
   // (1,-1)(0), etc.
};

// function that finds the decomposition into spherical harmonics at each radius
template <typename FLOAT, template <typename, typename> class sphericalmodel,
          template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange>
std::vector<FFTWpp::vector<std::complex<FLOAT>>>
TomographyModelTransform(
    const sphericalmodel<FLOAT, int> &backgroundmodel,
    GRID<FLOAT, OrderIndexRange, IndexRange> &grid,
    const CoefficientOrdering &whichsetup,
    const GaussQuad::Quadrature1D<FLOAT> &q,
    const std::vector<FLOAT> &vec_noderadii,
    const std::vector<FLOAT> &vec_alldepths,
    const std::vector<std::vector<std::complex<double>>> &vec_j) {
   using namespace GSHTrans;
   using Scalar = std::complex<FLOAT>;
   auto lMax = grid.MaxDegree();
   auto nMax = grid.MaxUpperIndex();
   // std::cout << "lmax: " << lMax << ", nMax: " << nMax << std::endl;

   auto rphys = [&vec_noderadii](int idxelem, double x) {
      return ((vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * x +
              (vec_noderadii[idxelem + 1] + vec_noderadii[idxelem])) *
             0.5;
   };

   // Make a random coefficient.
   auto getSize = [](auto lMax, auto n) {
      if constexpr (GSHTrans::RealFloatingPoint<Scalar>) {
         return GSHIndices<NonNegative>(lMax, lMax, n).size();
      } else {
         return GSHIndices<All>(lMax, lMax, n).size();
      }
   };

   auto size = getSize(lMax, 0);

   // std::vector<FFTWpp::vector<std::complex<FLOAT>>> vec_denslm;
   // sizes of various vectors:
   std::size_t nelem = vec_noderadii.size() - 1;
   int npoly = q.N() - 1;

   std::vector<FFTWpp::vector<std::complex<FLOAT>>> vec_denslm;
   vec_denslm.reserve(nelem * (npoly + 1));

   int laynum = 0;
   for (int idxelem = 0; idxelem < nelem;
        ++idxelem) {   // finding the layer number

      if (laynum < backgroundmodel.NumberOfLayers()) {
         if (!(vec_noderadii[idxelem] < backgroundmodel.UpperRadius(laynum))) {
            laynum += 1;
         };
      }

      for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
         // vector of values which will be transformed
         auto vec_dvs = FFTWpp::vector<Scalar>(
             grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes(), 1.0);
         auto dvslm = FFTWpp::vector<std::complex<FLOAT>>(size, 0.0);

         //////////////////////////////////////////////////////////////////
         //////////////////////////////////////////////////////////////////
         //////////////////////////////////////////////////////////////////
         // case of spherically symmetric background model
         for (int idxt = 0; idxt < grid.NumberOfCoLatitudes(); ++idxt) {
            for (int idxp = 0; idxp < grid.NumberOfLongitudes(); ++idxp) {

               auto idxuse = idxp + idxt * grid.NumberOfLongitudes();
               vec_dvs[idxuse] *=
                   backgroundmodel.Density(laynum)(
                       rphys(idxelem, q.X(idxpoly))) *
                   vec_j[idxelem * (npoly + 1) + idxpoly][idxuse];
            }
         }

         // transform to spherical harmonics
         grid.ForwardTransformation(lMax, 0, vec_dvs, dvslm);

         vec_denslm.push_back(dvslm);
      }
   }

   return vec_denslm;
};

// function that finds the decomposition into spherical harmonics at each radius
template <typename FLOAT, template <typename, typename> class sphericalmodel,
          template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange>
std::vector<FFTWpp::vector<std::complex<FLOAT>>>
TomographyModelTransformReferential(
    const sphericalmodel<FLOAT, int> &backgroundmodel,
    GRID<FLOAT, OrderIndexRange, IndexRange> &grid,
    const CoefficientOrdering &whichsetup,
    const GaussQuad::Quadrature1D<FLOAT> &q,
    const std::vector<FLOAT> &vec_noderadii,
    const std::vector<FLOAT> &vec_alldepths) {
   using namespace GSHTrans;
   using Scalar = std::complex<FLOAT>;
   auto lMax = grid.MaxDegree();
   auto nMax = grid.MaxUpperIndex();
   // std::cout << "lmax: " << lMax << ", nMax: " << nMax << std::endl;

   auto rphys = [&vec_noderadii](int idxelem, double x) {
      return ((vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) * x +
              (vec_noderadii[idxelem + 1] + vec_noderadii[idxelem])) *
             0.5;
   };

   // Make a random coefficient.
   auto getSize = [](auto lMax, auto n) {
      if constexpr (GSHTrans::RealFloatingPoint<Scalar>) {
         return GSHIndices<NonNegative>(lMax, lMax, n).size();
      } else {
         return GSHIndices<All>(lMax, lMax, n).size();
      }
   };

   auto size = getSize(lMax, 0);

   // std::vector<FFTWpp::vector<std::complex<FLOAT>>> vec_denslm;
   // sizes of various vectors:
   std::size_t nelem = vec_noderadii.size() - 1;
   int npoly = q.N() - 1;

   std::vector<FFTWpp::vector<std::complex<FLOAT>>> vec_denslm;
   vec_denslm.reserve(nelem * (npoly + 1));

   int laynum = 0;
   for (int idxelem = 0; idxelem < nelem;
        ++idxelem) {   // finding the layer number

      if (laynum < backgroundmodel.NumberOfLayers()) {
         if (!(vec_noderadii[idxelem] < backgroundmodel.UpperRadius(laynum))) {
            laynum += 1;
         };
      }

      for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
         // vector of values which will be transformed
         auto vec_dvs = FFTWpp::vector<Scalar>(
             grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes(), 1.0);
         auto dvslm = FFTWpp::vector<std::complex<FLOAT>>(size, 0.0);

         //////////////////////////////////////////////////////////////////
         //////////////////////////////////////////////////////////////////
         //////////////////////////////////////////////////////////////////
         // case of spherically symmetric background model
         for (int idxt = 0; idxt < grid.NumberOfCoLatitudes(); ++idxt) {
            for (int idxp = 0; idxp < grid.NumberOfLongitudes(); ++idxp) {

               auto idxuse = idxp + idxt * grid.NumberOfLongitudes();
               vec_dvs[idxuse] *= backgroundmodel.Density(laynum)(
                   rphys(idxelem, q.X(idxpoly)));
            }
         }

         // transform to spherical harmonics
         grid.ForwardTransformation(lMax, 0, vec_dvs, dvslm);

         vec_denslm.push_back(dvslm);
      }
   }

   return vec_denslm;
};

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
    const int &lMax) {
   std::size_t nelem = vec_noderadii.size() - 1;
   int npoly = q.N() - 1;
   int matlen = nelem * npoly + 1;   // size of matrix
   Eigen::VectorXcd vec_fullforce =
       Eigen::VectorXcd::Zero(std::pow(lMax + 1, 2) * matlen);

   // finding \rho \mathbf{s} in canonical components spatially
   // firstly we need to have the background model density and the perturbation
   // from the perturbing model:
   // go through all layers, find what s is, as we currently just have radial,
   // that should be fairly simple
   int totnum = grid.NumberOfCoLatitudes() * grid.NumberOfLongitudes();
   std::vector<std::vector<std::vector<std::complex<double>>>> vec_s(
       nelem * (npoly + 1),
       std::vector<std::vector<std::complex<double>>>(
           totnum, std::vector<std::complex<double>>(3, 0.0)));
   int laynum = 0;
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // finding the layer number
      if (laynum < backgroundmodel.NumberOfLayers()) {
         if (!(vec_noderadii[idxelem] < backgroundmodel.UpperRadius(laynum))) {
            laynum += 1;
         };
      }
      for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
         auto rval = rphys(idxelem, q.X(idxpoly));
         auto raddensity = backgroundmodel.Density(laynum)(rval);
         for (int idxt = 0; idxt < grid.NumberOfCoLatitudes(); ++idxt) {
            for (int idxp = 0; idxp < grid.NumberOfLongitudes(); ++idxp) {
               auto idxuse = idxt * grid.NumberOfLongitudes() + idxp;
               vec_s[idxelem * (npoly + 1) + idxpoly][idxuse][1] =
                   pertmodel.RadialMap(rval, grid.CoLatitudes()[idxt],
                                       grid.Longitudes()[idxp]) *
                   raddensity;
            }
         }
      }
   }

   return vec_fullforce;
};

// function for finding the force vector:
template <typename FLOAT, template <typename, typename> class sphericalmodel>
Eigen::VectorXcd
FindForce(const sphericalmodel<FLOAT, int> &backgroundmodel,
          const CoefficientOrdering &whichsetup,
          const GaussQuad::Quadrature1D<FLOAT> &q,
          const std::vector<FLOAT> &vec_noderadii,
          const std::vector<FFTWpp::vector<std::complex<FLOAT>>> &vec_denslm,
          const int &lMax) {
   // auto intsize = grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();
   // lambda to use in mapping [-1,1] to [idxelem]
   auto rscale = [&vec_noderadii](int idxelem, double x) {
      return (x + (vec_noderadii[idxelem + 1] + vec_noderadii[idxelem]) /
                      (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]));
   };
   std::size_t nelem = vec_noderadii.size() - 1;
   int npoly = q.N() - 1;
   int matlen = nelem * npoly + 1;   // size of matrix
   Eigen::VectorXcd vec_fullforce =
       Eigen::VectorXcd::Zero(std::pow(lMax + 1, 2) * matlen);

   int laynum = 0;   // layer number
   // constants
   const double bigg_db = 6.6743 * std::pow(10.0, -11.0);
   const double pi_db = 3.1415926535;
   double multfact =
       4.0 * pi_db * bigg_db * std::pow(backgroundmodel.LengthNorm(), 2.0);

   // looping over elements
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // finding the layer number
      // if (laynum < backgroundmodel.NumberOfLayers()) {
      //    if (!(vec_noderadii[idxelem] < backgroundmodel.UpperRadius(laynum)))
      //    {
      //       laynum += 1;
      //    };
      // }

      // looping over quadrature nodes
      for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
         // looping over l,m values
         for (int idxl = 0; idxl < lMax + 1; ++idxl) {
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {

               std::size_t mynewidx =
                   std::pow(idxl, 2) + idxl + idxm +
                   (idxpoly + npoly * idxelem) * std::pow(lMax + 1, 2);

               vec_fullforce(mynewidx) +=
                   multfact * q.W(idxpoly) *
                   std::pow(0.5 * (vec_noderadii[idxelem + 1] -
                                   vec_noderadii[idxelem]),
                            3.0) *
                   std::pow(rscale(idxelem, q.X(idxpoly)), 2.0) *
                   vec_denslm[idxelem * (npoly + 1) + idxpoly]
                             [std::pow(idxl, 2) + idxl + idxm];
            }
         }
      };
   };

   // correcting at zero radius
   for (int idxl = 1; idxl < lMax + 1; ++idxl) {
      for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
         std::size_t mynewidx = std::pow(idxl, 2) + idxl + idxm;
         vec_fullforce(mynewidx) *= 0;
      }
   }
   std::size_t idxmaxval = std::pow(lMax + 1, 2) * matlen;
   for (auto idx1 = 0; idx1 < idxmaxval; ++idx1) {
      vec_fullforce(idx1) *= backgroundmodel.DensityNorm();
   }

   return vec_fullforce;
};

template <typename FLOAT, template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange>
std::vector<Eigen::MatrixXcd>
SphdecompToGridSol(const Eigen::VectorXcd &mybigsol,
                   GRID<FLOAT, OrderIndexRange, IndexRange> &grid,
                   const int &nelem, const int &npoly,
                   const std::size_t &size) {
   using Complex = std::complex<FLOAT>;
   using Scalar = Complex;
   std::vector<Eigen::MatrixXcd> vec_gridsol;
   {

      for (int idxelem = 0; idxelem < nelem; ++idxelem) {
         // finding the maximum index for the inner loop
         int maxpolyidx;
         if (idxelem < nelem - 1) {
            maxpolyidx = npoly;
         } else if (idxelem == nelem - 1) {
            maxpolyidx = npoly + 1;
         };

         // loop over the nodes inside each element
         for (int idxpoly = 0; idxpoly < maxpolyidx; ++idxpoly) {
            // temporary matrix to contain the grid values
            Eigen::MatrixXcd tmpmat;
            tmpmat.resize(grid.NumberOfCoLatitudes(),
                          grid.NumberOfLongitudes());

            // temporary vector to contain lm coefficients and its values
            auto dvslm = FFTWpp::vector<Complex>(size, 0.0);
            for (int idx = 0; idx < size; ++idx) {
               int myidx = (idxelem * npoly + idxpoly) * size + idx;
               dvslm[idx] = mybigsol(myidx);
            }

            // temporary vector to transform to
            auto vec_sol = FFTWpp::vector<Scalar>(grid.NumberOfLongitudes() *
                                                  grid.NumberOfCoLatitudes());
            grid.InverseTransformation(grid.MaxDegree(), 0, dvslm, vec_sol);

            // filling out the temporary matrix which will be inserted
            for (int idxt = 0; idxt < grid.NumberOfCoLatitudes(); ++idxt) {
               for (int idxp = 0; idxp < grid.NumberOfLongitudes(); ++idxp) {
                  tmpmat(idxt, idxp) =
                      vec_sol[idxt * grid.NumberOfLongitudes() + idxp];
               };
            };
            vec_gridsol.push_back(tmpmat);
         };
      };
   }

   return vec_gridsol;
};

template <typename FLOAT, template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange,
          template <typename, typename> class TOMPERTMODEL>
auto
tom_to_h(GRID<FLOAT, OrderIndexRange, IndexRange> &grid, std::size_t npoly,
         std::vector<double> &vec_noderadii,
         const GaussQuad::Quadrature1D<FLOAT> &q,
         const TOMPERTMODEL<FLOAT, int> &mypertprem) {
   auto mytheta = grid.CoLatitudes();
   auto myphi = grid.Longitudes();
   auto intsize = grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();
   auto nelem = vec_noderadii.size() - 1;
   std::vector<std::vector<std::complex<double>>> vec_h(nelem * npoly + 1);

   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // getting max idxpoly
      int idxpolymax = npoly + 1;
      if (idxelem < nelem - 1) {
         idxpolymax = npoly;
      } else {
         idxpolymax = npoly + 1;
      }

      // looping through index poly
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
   return vec_h;
};

template <typename FLOAT, template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange,
          template <typename, typename> class TOMMODEL,
          template <typename, typename> class TOMPERTMODEL>
auto
tom_to_h(GRID<FLOAT, OrderIndexRange, IndexRange> &grid, std::size_t npoly,
         std::vector<double> &vec_noderadii,
         const GaussQuad::Quadrature1D<FLOAT> &q,
         const TOMPERTMODEL<FLOAT, int> &mypertprem,
         const TOMMODEL<FLOAT, int> &backgroundmodel) {
   auto mytheta = grid.CoLatitudes();
   auto myphi = grid.Longitudes();
   auto intsize = grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();
   auto nelem = vec_noderadii.size() - 1;
   auto ballrad = vec_noderadii.back();
   std::vector<std::vector<std::complex<double>>> vec_h(nelem * npoly + 1);

   FLOAT outerrad = backgroundmodel.OuterRadius();
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // getting max idxpoly
      int idxpolymax = npoly + 1;
      if (idxelem < nelem - 1) {
         idxpolymax = npoly;
      } else {
         idxpolymax = npoly + 1;
      }

      // looping through index poly
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         std::vector<std::complex<double>> tmp_h(intsize);
         auto radr = GravityFunctions::StandardIntervalMap(
             q.X(idxpoly), vec_noderadii[idxelem], vec_noderadii[idxelem + 1]);
         if (radr <= outerrad) {
            for (int idxtheta = 0; idxtheta < grid.NumberOfCoLatitudes();
                 ++idxtheta) {
               for (int idxphi = 0; idxphi < grid.NumberOfLongitudes();
                    ++idxphi) {
                  tmp_h[idxtheta * grid.NumberOfLongitudes() + idxphi] =
                      mypertprem.RadialMap(radr, mytheta[idxtheta],
                                           myphi[idxphi]);
               }
            }
         } else {
            for (int idxtheta = 0; idxtheta < grid.NumberOfCoLatitudes();
                 ++idxtheta) {
               for (int idxphi = 0; idxphi < grid.NumberOfLongitudes();
                    ++idxphi) {
                  tmp_h[idxtheta * grid.NumberOfLongitudes() + idxphi] =
                      mypertprem.RadialMap(outerrad, mytheta[idxtheta],
                                           myphi[idxphi]) *
                      (ballrad - radr) / (ballrad - outerrad);
               }
            }
         }

         vec_h[idxelem * npoly + idxpoly] = tmp_h;
      }
   }
   return vec_h;
};

template <typename FLOAT, template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange>
auto
vecgrid_to_sphdecomp(
    GRID<FLOAT, OrderIndexRange, IndexRange> &grid,
    const std::vector<std::vector<std::complex<FLOAT>>> &vec_spatial,
    std::size_t npoly, std::size_t nelem) {
   auto lMax = grid.MaxDegree();
   auto nMax = grid.MaxUpperIndex();
   auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
   int matlen = nelem * npoly + 1;   // size of matrix
   auto intsize = grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();
   using veccomp = std::vector<std::complex<FLOAT>>;
   using vecvech = std::vector<veccomp>;

   // declaring variables used in calculation
   vecvech vec_sphharm(matlen, veccomp(_size0, 0.0));
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // max for idxpoly loop
      int idxpolymax;
      if (idxelem < nelem - 1) {
         idxpolymax = npoly;
      } else {
         idxpolymax = npoly + 1;
      }

      // finding h0 component
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {

         grid.ForwardTransformation(lMax, 0,
                                    vec_spatial[idxelem * npoly + idxpoly],
                                    vec_sphharm[idxelem * npoly + idxpoly]);
      }
   };
   return vec_sphharm;
};

template <typename FLOAT, template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange>
auto
sphdecomp_to_vecgrid(
    GRID<FLOAT, OrderIndexRange, IndexRange> &grid,
    const std::vector<std::vector<std::complex<FLOAT>>> &vec_sphdecomp,
    std::size_t npoly, std::size_t nelem) {
   auto lMax = grid.MaxDegree();
   auto nMax = grid.MaxUpperIndex();
   auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
   int matlen = nelem * npoly + 1;   // size of matrix
   auto intsize = grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();
   using veccomp = std::vector<std::complex<FLOAT>>;
   using vecvech = std::vector<veccomp>;

   // declaring variables used in calculation
   vecvech vec_grid(matlen, veccomp(intsize, 0.0));
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // max for idxpoly loop
      int idxpolymax;
      if (idxelem < nelem - 1) {
         idxpolymax = npoly;
      } else {
         idxpolymax = npoly + 1;
      }
      // finding h0 component
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {

         grid.InverseTransformation(lMax, 0,
                                    vec_sphdecomp[idxelem * npoly + idxpoly],
                                    vec_grid[idxelem * npoly + idxpoly]);
      }
   };
   return vec_grid;
};

template <typename FLOAT, template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange>
auto
h_to_j(GRID<FLOAT, OrderIndexRange, IndexRange> &grid, std::size_t npoly,
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
   // auto vec_hlm = vecgrid_to_sphdecomp(grid, vec_h, npoly, nelem);
   auto mat_0 = MATRIX3::Zero();
   std::vector<std::vector<std::complex<double>>> vec_j(
       nelem * (npoly + 1), std::vector<std::complex<double>>(intsize, 0.0));

   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      // max for idxpoly loop
      int idxpolymax = npoly + 1;
      int idxpolymin = 0;
      if (idxelem == 0) {
         idxpolymin = 1;
      }

      // components of derivative, ie \nabla h:
      vecvech vec_h0(npoly + 1, veccomp(intsize, 0.0));
      double inv2 = 2.0 / vec_elemwidth[idxelem];

      // finding h0 component
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         for (int idxrad = 0; idxrad < intsize; ++idxrad) {
            for (int idxn = 0; idxn < npoly + 1; ++idxn) {
               auto idxouter = idxelem * npoly + idxn;
               vec_h0[idxpoly][idxrad] += vec_h[idxouter][idxrad] *
                                          _mat_gaussderiv(idxn, idxpoly) * inv2;
            }
         }
      }

      // std::cout << "Hello 5 \n";

      // find a
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
         for (int idxrad = 0; idxrad < intsize; ++idxrad) {
            auto idxuse = idxelem * (npoly + 1) + idxpoly;
            std::complex<double> hdivr;
            if (idxuse == 0) {
               hdivr = 0.0;
            } else {
               hdivr = vec_h[idxelem * npoly + idxpoly][idxrad] /
                       vec_allradii[idxelem * npoly + idxpoly];
            }
            vec_j[idxuse][idxrad] =
                (1.0 + hdivr) * (1.0 + hdivr) * (1.0 + vec_h0[idxpoly][idxrad]);
         }
      }
   }

   return vec_j;
};

template <typename FLOAT, template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange, class sphericalmodel>
auto
h_to_a(GRID<FLOAT, OrderIndexRange, IndexRange> &grid,
       const sphericalmodel &inp_model) {
   return 1;
};

template <typename FLOAT, template <typename, typename, typename> class GRID,
          typename OrderIndexRange, typename IndexRange>
auto
h_to_a(GRID<FLOAT, OrderIndexRange, IndexRange> &grid, std::size_t npoly,
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
   // for (int idx1 = 0; idx1 < npoly + 1; ++idx1) {
   //    // std::cout << "Difference: " << vec_h[1][idx1] - vec_hcheck[1][idx1]
   //    //  << std::endl;
   //    std::cout << "h: " << vec_h[idx1][0] << ", hlm: "
   //              << vec_hlm[idx1][0] / std::pow(4.0 * 3.1415926535, 0.5)
   //              << std::endl;
   // }

   // std::cout << std::endl;
   auto mat_0 = MATRIX3::Zero();
   std::vector<std::vector<MATRIX3>> vec_a(
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
      //   vec_h0(npoly + 1, veccomp(_size0, 0.0));
      double inv2 = 2.0 / vec_elemwidth[idxelem];

      // finding h0 component
      // for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {

      //    for (int idxl = 0; idxl < lMax + 1; ++idxl) {
      //       for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
      //          auto idxref = idxm + idxl * (1 + idxl);

      //          for (int idxn = 0; idxn < npoly + 1; ++idxn) {
      //             auto idxouter = idxelem * npoly + idxn;
      //             vec_h0[idxpoly][idxref] += vec_hlm[idxouter][idxref] *
      //                                        _mat_gaussderiv(idxn, idxpoly) *
      //                                        inv2;
      //          }
      //       }
      //    }
      // }
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
         // grid.InverseTransformation(lMax, 0, vec_h0[idxpoly],
         // vec_nh0[idxpoly]);
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
            // auto nh0 = vec_nh0[idxpoly][idxrad];
            auto nh0 = vec_h02[idxpoly][idxrad];
            std::complex<double> hdivr;
            if (idxuse == 0) {
               hdivr = 0.0;
            } else {
               hdivr = vec_h[idxelem * npoly + idxpoly][idxrad] /
                       vec_allradii[idxelem * npoly + idxpoly];
            }
            // if (idxuse == 1 && idxrad < 5) {
            //    std::cout << hdivr << std::endl;
            // }
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

   return vec_a;
};

template <typename FLOAT>
Eigen::Matrix<FLOAT, Eigen::Dynamic, 1>
Gravity1D(const GaussQuad::Quadrature1D<FLOAT> &q,
          const Eigen::Matrix<FLOAT, Eigen::Dynamic, 1> &vecsol,
          const std::vector<FLOAT> &vec_noderadii, const int &nelem,
          const int &npoly, const FLOAT &rad_scale) {

   auto pleg =
       Interpolation::LagrangePolynomial(q.Points().begin(), q.Points().end());
   Eigen::Matrix<FLOAT, Eigen::Dynamic, 1> vecderiv =
       Eigen::Matrix<FLOAT, Eigen::Dynamic, 1>::Zero(nelem + 1);
   for (int idxelem = 0; idxelem < nelem + 1; ++idxelem) {
      if (idxelem < nelem) {
         for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
            vecderiv(idxelem) +=
                vecsol[idxelem * npoly + idxpoly] * 2.0 /
                (vec_noderadii[idxelem + 1] - vec_noderadii[idxelem]) *
                pleg.Derivative(idxpoly, -1.0);
         };
      } else {
         for (int idxpoly = 0; idxpoly < npoly + 1; ++idxpoly) {
            vecderiv(idxelem) +=
                vecsol[(idxelem - 1) * npoly + idxpoly] * 2.0 /
                (vec_noderadii[idxelem] - vec_noderadii[idxelem - 1]) *
                pleg.Derivative(idxpoly, 1.0);
         };
      };
   };
   vecderiv *= 1.0 / rad_scale;
   return vecderiv;
};

// class for the spherical integrator, with the sparse matrix and associated
// routines
template <typename FLOAT, template <typename, typename> class sphericalmodel>
   requires PlanetaryModel::SphericalElasticModel<sphericalmodel<FLOAT, int>>
class Spherical_Integrator {
 private:
   using CFLOAT = std::complex<FLOAT>;
   // using CFLOAT = FLOAT;

   using SPHMAT = Eigen::Matrix<FLOAT, Eigen::Dynamic, Eigen::Dynamic>;
   using SPHVEC = Eigen::Matrix<CFLOAT, Eigen::Dynamic, 1>;
   using SPARSEMAT = Eigen::SparseMatrix<CFLOAT>;
   using CIT_FLOAT = std::vector<FLOAT>::const_iterator;
   using IT_FLOAT = std::vector<FLOAT>::iterator;

 public:
   // constructor

   Spherical_Integrator(const sphericalmodel<FLOAT, int> &mymodel,
                        const std::vector<int> &nnodes, const int &ord_poly,
                        const int &l_val)
       : vec_nnodes{nnodes}, polyord{ord_poly}, lval{l_val},
         nlayer{mymodel.NumberOfLayers()}, rad_scale{mymodel.LengthNorm()},
         m_isInitialized{true} {
      // quadrature
      q = GaussQuad::GaussLobattoLegendreQuadrature1D<FLOAT>(polyord + 1);

      using namespace Interpolation;

      // finding radial nodes
      vec_noderadii = Layer_Node(mymodel, nnodes);

      nelem = vec_noderadii.size() - 1;   // total number of elements
      std::vector<FLOAT> nodespacing(nelem), invnodespacingtwo(nelem);
      for (int idxelem = 0; idxelem < nelem; ++idxelem) {
         nodespacing[idxelem] =
             vec_noderadii[idxelem + 1] - vec_noderadii[idxelem];
         invnodespacingtwo[idxelem] = 2.0 / nodespacing[idxelem];
      };
      matlen = nelem * polyord + 1;   // size of matrix

      using T = Eigen::Triplet<CFLOAT>;
      auto pleg = Interpolation::LagrangePolynomial(q.Points().begin(),
                                                    q.Points().end());
      vec_derivatives.reserve(std::pow((polyord + 1), 2) * (polyord + 2) / 2);
      for (int idxi = 0; idxi < polyord + 1; ++idxi) {
         for (int idxj = 0; idxj < idxi + 1; ++idxj) {
            for (int idxk = 0; idxk < polyord + 1; ++idxk) {
               vec_derivatives.push_back(pleg.Derivative(idxi, q.X(idxk)) *
                                         pleg.Derivative(idxj, q.X(idxk)));
            };
         };
      };

      // lambda to find value of r^2 l_i'(x_k) l_j'(x_k), with x_k as the
      // Gaussian quadrature point
      auto funcnew = [this, &invnodespacingtwo](int idxelem, int idxi, int idxj,
                                                int idxk) {
         return invnodespacingtwo[idxelem] *
                std::pow(StandardIntervalMap(q.X(idxk), vec_noderadii[idxelem],
                                             vec_noderadii[idxelem + 1]),
                         2.0) *
                this->vec_derivatives[this->idxderiv(idxi, idxj, idxk)];
      };

      std::vector<T> tripletlist;
      // tripletlist.reserve(nelem * (polyord + 1) * (polyord + 1) + 1);
      tripletlist.reserve(nelem * (((polyord + 1) * (polyord + 2)) / 2 - 1) +
                          1);

      {

         FLOAT zetasq = static_cast<FLOAT>(lval) * static_cast<FLOAT>(lval + 1);

         // finding the integrals using the inner product and the vector with
         // the derivative values at the interpolation points
         for (int idxelem = 0; idxelem < nelem; ++idxelem) {
            FLOAT tmp;
            for (int idxi = 0; idxi < polyord + 1; ++idxi) {
               for (int idxj = 0; idxj < idxi + 1; ++idxj) {
                  tmp = 0.0;
                  for (int idxk = 0; idxk < polyord + 1; ++idxk) {
                     tmp -= funcnew(idxelem, idxi, idxj, idxk) * q.W(idxk);
                  }

                  // push back into triplet list
                  if (idxj < idxi) {
                     tripletlist.push_back(T(idxelem * polyord + idxi,
                                             idxelem * polyord + idxj, tmp));
                  } else {
                     tripletlist.push_back(T(idxelem * polyord + idxi,
                                             idxelem * polyord + idxj,
                                             tmp - zetasq * q.W(idxi)));
                  }
               };
            };
         };

         // final point
         tripletlist.push_back(
             T(matlen - 1, matlen - 1,
               -static_cast<CFLOAT>(lval + 1) *
                   static_cast<CFLOAT>(vec_noderadii[nelem])));

         // construct sparse matrix
         specelem.resize(matlen, matlen);   // resize
         specelem.setFromTriplets(tripletlist.begin(),
                                  tripletlist.end());   // set it
         specelem.makeCompressed();                     // compress

         // Cholesky decomposition
         chol_solver.compute(specelem);
      };
   };

   SPHVEC solve(const SPHVEC &vec_force) {
      assert(m_isInitialized && "Not initialized");
      // std::cout << "Rows: " << vec_force.rows() << ". Columns: " <<
      // vec_force.cols() << std::endl; std::cout << "Rows: " <<
      // specelem.rows() << ". Columns: " << specelem.cols() << std::endl;
      auto vecsol = chol_solver.solve(vec_force);

      return vecsol;
   };

   // Spherical_Integrator(sphericalmodel mymodel, ){};

 private:
   SPARSEMAT specelem;   // vector of spectral element matrices
   int polyord, lval, nlayer, nelem, matlen;
   std::vector<int> vec_nnodes;
   FLOAT rad_scale;
   // node values:
   std::vector<FLOAT> vec_noderadii;
   std::vector<FLOAT> vec_derivatives;

   // Gauss quadrature values
   GaussQuad::Quadrature1D<FLOAT> q;
   Eigen::SimplicialLDLT<SPARSEMAT, Eigen::Lower, Eigen::AMDOrdering<int>>
       chol_solver;
   // Eigen::CholmodSupernodalLLT<SPARSEMAT, Eigen::Lower> chol_solver;
   bool m_isInitialized;

   // return the index in the vector of derivatives for function l_i'(x_k)
   // l_j'(x_k). Stored such that contiguous in memory over k:
   int idxderiv(const int &i, const int &j, const int &k) {
      assert(m_isInitialized && "Not initialized");
      int n = this->polyord + 1;
      // return k + n * (j + i * (n - 1) - (i * (i - 1)) / 2);
      return k + n * (j + (i * (i + 1)) / 2);
   };
};   // end spherical integrator class

template <typename FLOAT, template <typename, typename> class sphericalmodel>
   requires PlanetaryModel::SphericalElasticModel<sphericalmodel<FLOAT, int>>
class SphericalGeometryPoissonSolver {
 private:
   using CFLOAT = std::complex<FLOAT>;
   // using CFLOAT = FLOAT;

   using SPHMAT = Eigen::Matrix<FLOAT, Eigen::Dynamic, Eigen::Dynamic>;
   using SPHVEC = Eigen::Matrix<CFLOAT, Eigen::Dynamic, 1>;
   using SPARSEMAT = Eigen::SparseMatrix<CFLOAT>;
   using CIT_FLOAT = std::vector<FLOAT>::const_iterator;
   using IT_FLOAT = std::vector<FLOAT>::iterator;

 public:
   SphericalGeometryPoissonSolver(const sphericalmodel<FLOAT, int> &mymodel,
                                  const CoefficientOrdering &whichsetup,
                                  const std::vector<FLOAT> &noderadii,
                                  const int &ord_poly, const int &l_max)
       : vec_noderadii{noderadii}, polyord{ord_poly}, lmax{l_max},
         nlayer{mymodel.NumberOfLayers()}, rad_scale{mymodel.LengthNorm()},
         sparse_ordering{whichsetup}, m_isInitialized{true} {
      // quadrature
      q = GaussQuad::GaussLobattoLegendreQuadrature1D<FLOAT>(polyord + 1);

      using namespace Interpolation;

      // finding radial nodes
      // vec_noderadii = Layer_Node(mymodel, nnodes);

      nelem = vec_noderadii.size() - 1;   // total number of elements
      std::vector<FLOAT> nodespacing(nelem), invnodespacingtwo(nelem);
      for (int idxelem = 0; idxelem < nelem; ++idxelem) {
         nodespacing[idxelem] =
             vec_noderadii[idxelem + 1] - vec_noderadii[idxelem];
         invnodespacingtwo[idxelem] = 2.0 / nodespacing[idxelem];
      };
      matlen = nelem * polyord + 1;                  // size of single l matrix
      matfulllen = matlen * std::pow(lmax + 1, 2);   // total size of matrix

      // fill out derivative vector
      auto pleg = Interpolation::LagrangePolynomial(q.Points().begin(),
                                                    q.Points().end());
      vec_derivatives.reserve(std::pow((polyord + 1), 2) * (polyord + 2) / 2);
      for (int idxi = 0; idxi < polyord + 1; ++idxi) {
         for (int idxj = 0; idxj < idxi + 1; ++idxj) {
            for (int idxk = 0; idxk < polyord + 1; ++idxk) {
               vec_derivatives.push_back(pleg.Derivative(idxi, q.X(idxk)) *
                                         pleg.Derivative(idxj, q.X(idxk)));
            };
         };
      };

      // lambda to find value of r^2 l_i'(x_k) l_j'(x_k), with x_k as the
      // Gaussian quadrature point
      auto funcnew = [this, &invnodespacingtwo](int idxelem, int idxi, int idxj,
                                                int idxk) {
         return invnodespacingtwo[idxelem] *
                std::pow(StandardIntervalMap(q.X(idxk), vec_noderadii[idxelem],
                                             vec_noderadii[idxelem + 1]),
                         2.0) *
                this->vec_derivatives[this->idxderiv(idxi, idxj, idxk)];
      };

      // triplet list for constructing the sparse matrix
      using T = Eigen::Triplet<CFLOAT>;
      std::vector<T> tripletlist;

      tripletlist.reserve(
          (nelem * (((polyord + 1) * (polyord + 2)) / 2 - 1) + 1) *
          std::pow(lmax + 1, 2));

      for (int idxl = 0; idxl < lmax + 1; ++idxl) {
         // zeta square
         FLOAT zetasq = static_cast<FLOAT>(idxl) * static_cast<FLOAT>(idxl + 1);

         // finding the integrals using the inner product and the vector with
         // the derivative values at the interpolation points
         for (int idxelem = 0; idxelem < nelem; ++idxelem) {
            FLOAT tmp;
            for (int idxi = 0; idxi < polyord + 1; ++idxi) {
               for (int idxj = 0; idxj < idxi + 1; ++idxj)
               // for (int idxj = 0; idxj < polyord + 1; ++idxj)
               {
                  tmp = 0.0;
                  for (int idxk = 0; idxk < polyord + 1; ++idxk) {
                     tmp -= funcnew(idxelem, idxi, idxj, idxk) * q.W(idxk);
                  }
                  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                     std::size_t idx_x, idx_y;
                     idx_x = xidx(idxelem, idxi, idxj, idxl, idxm);
                     idx_y = yidx(idxelem, idxi, idxj, idxl, idxm);
                     // push back into triplet list

                     if (idxj < idxi) {
                        tripletlist.push_back(T(idx_x, idx_y, tmp));
                     } else {
                        tripletlist.push_back(
                            T(idx_x, idx_y,
                              tmp - zetasq * q.W(idxi) * nodespacing[idxelem] *
                                        0.5));
                        // std::cout << "Hello" << std::abs(tmp) << " "
                        // << std::abs(zetasq * q.W(idxi));
                     }
                     ///////////////////////////////////////////
                  };
               };
            };
         };

         for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
            std::size_t idx_x, idx_y;
            idx_x = xidx(nelem - 1, polyord, polyord, idxl, idxm);
            idx_y = yidx(nelem - 1, polyord, polyord, idxl, idxm);

            // push back final point into triplet list
            tripletlist.push_back(
                T(idx_x, idx_y,
                  -static_cast<CFLOAT>(idxl + 1) *
                      static_cast<CFLOAT>(vec_noderadii[nelem])));
         };
      };
      // construct sparse matrix
      sphspecmat.resize(matfulllen, matfulllen);   // resize
      sphspecmat.setFromTriplets(tripletlist.begin(),
                                 tripletlist.end());   // set it
      sphspecmat.makeCompressed();                     // compress

      // this->compute(sphspecmat);
      // vectest.push_back(chol_solver);
   };

   SphericalGeometryPoissonSolver &compute() {
      assert(m_isInitialized && "Not initialized");
      chol_solver.compute(this->sphspecmat);
      m_cholisInitialized = true;
      return *this;
   };
   std::size_t matrix_length() {
      assert(m_isInitialized && "Not initialized");
      return matfulllen;
   };

   SPARSEMAT const specmat() {
      assert(m_isInitialized && "Not initialized");
      return sphspecmat;
   };

   SPHVEC solve(const SPHVEC &vec_force) {
      assert(m_cholisInitialized && "Solver not initialized");
      // std::cout << "Rows: " << vec_force.rows() << ". Columns: " <<
      // vec_force.cols() << std::endl; std::cout << "Rows: " <<
      // specelem.rows() << ". Columns: " << specelem.cols() << std::endl;
      // SPHVEC vecsol = chol_solver.solve(vec_force);
      SPHVEC vecsol = chol_solver.solve(vec_force);

      return vecsol;
   };

   // template <typename FLOATVEC>
   SPHVEC sphmatmult(const SPHVEC &lhs) {
      assert(m_isInitialized && "Not initialized");
      SPHVEC vec_out;
      vec_out = sphspecmat.template selfadjointView<Eigen::Lower>() * lhs;
      // vec_out = sphspecmat * lhs;
      return vec_out;
   };

 private:
   SPARSEMAT sphspecmat;   // vector of spectral element matrices
   int polyord, lmax, nlayer, nelem, matlen;
   std::size_t matfulllen;
   // std::vector<int> vec_nnodes;
   FLOAT rad_scale;
   // node values:
   std::vector<FLOAT> vec_noderadii;
   std::vector<FLOAT> vec_derivatives;

   // Gauss quadrature values
   GaussQuad::Quadrature1D<FLOAT> q;
   Eigen::SimplicialLDLT<SPARSEMAT, Eigen::Lower> chol_solver;
   // Eigen::SparseLU<SPARSEMAT> chol_solver;
   // std::vector<Eigen::SimplicialLDLT<SPARSEMAT, Eigen::Lower,
   // Eigen::AMDOrdering<int>>> vectest;
   bool m_isInitialized = false;
   bool m_cholisInitialized = false;
   CoefficientOrdering sparse_ordering;

   // index of vector
   int idxderiv(const int &i, const int &j, const int &k) {
      assert(m_isInitialized && "Not initialized");
      int n = this->polyord + 1;
      // return k + n * (j + i * (n - 1) - (i * (i - 1)) / 2);
      return k + n * (j + (i * (i + 1)) / 2);
   };

   // index for sparse matrix
   std::size_t xidx(const int &idxelem, const int &idxi, const int &idxj,
                    const int &idxl, const int &idxm) {
      assert(m_isInitialized && "Not initialized");
      if (sparse_ordering == CoefficientOrdering::SphericalDecomp) {
         return (std::pow(idxl, 2) + idxl + idxm) * this->matlen + idxi +
                idxelem * this->polyord;
      } else if (sparse_ordering == CoefficientOrdering::RadialClumped) {
         return std::pow(idxl, 2) + idxl + idxm +
                (idxi + this->polyord * idxelem) * std::pow(this->lmax + 1, 2);
      } else {
         return 0;
      };
   };

   std::size_t yidx(const int &idxelem, const int &idxi, const int &idxj,
                    const int &idxl, const int &idxm) {
      assert(m_isInitialized && "Not initialized");
      if (sparse_ordering == CoefficientOrdering::SphericalDecomp) {
         return (std::pow(idxl, 2) + idxl + idxm) * this->matlen + idxj +
                idxelem * this->polyord;
      } else if (sparse_ordering == CoefficientOrdering::RadialClumped) {
         return std::pow(idxl, 2) + idxl + idxm +
                (idxj + this->polyord * idxelem) * std::pow(this->lmax + 1, 2);
      } else {
         return 0;
      };
   };
};

// class for the spherical integrator, with the sparse matrix and associated
// routines
template <typename FLOAT, template <typename, typename> class sphericalmodel>
   requires PlanetaryModel::SphericalElasticModel<sphericalmodel<FLOAT, int>>
class PoissonSphericalHarmonic {
 private:
   using CFLOAT = std::complex<FLOAT>;
   // using CFLOAT = FLOAT;

   using SPHMAT = Eigen::Matrix<FLOAT, Eigen::Dynamic, Eigen::Dynamic>;
   using SPHVEC = Eigen::Matrix<CFLOAT, Eigen::Dynamic, 1>;
   using SPARSEMAT = Eigen::SparseMatrix<CFLOAT>;
   using EIGCHOL =
       Eigen::SimplicialLDLT<SPARSEMAT, Eigen::Lower, Eigen::AMDOrdering<int>>;
   using CIT_FLOAT = std::vector<FLOAT>::const_iterator;
   using IT_FLOAT = std::vector<FLOAT>::iterator;

 public:
   // constructor

   PoissonSphericalHarmonic(const sphericalmodel<FLOAT, int> &mymodel,
                            const std::vector<int> &nnodes, const int &ord_poly,
                            const int &l_max)
       : vec_nnodes{nnodes}, polyord{ord_poly}, lmax{l_max},
         nlayer{mymodel.NumberOfLayers()}, rad_scale{mymodel.LengthNorm()},
         m_isInitialized{true} {
      // quadrature
      q = GaussQuad::GaussLobattoLegendreQuadrature1D<FLOAT>(polyord + 1);

      using namespace Interpolation;

      // finding radial nodes
      vec_noderadii = Layer_Node(mymodel, nnodes);

      nelem = vec_noderadii.size() - 1;   // total number of elements
      std::vector<FLOAT> nodespacing(nelem), invnodespacingtwo(nelem);
      for (int idxelem = 0; idxelem < nelem; ++idxelem) {
         nodespacing[idxelem] =
             vec_noderadii[idxelem + 1] - vec_noderadii[idxelem];
         invnodespacingtwo[idxelem] = 2.0 / nodespacing[idxelem];
      };
      matlen = nelem * polyord + 1;   // size of matrix

      using T = Eigen::Triplet<CFLOAT>;
      auto pleg = Interpolation::LagrangePolynomial(q.Points().begin(),
                                                    q.Points().end());
      vec_derivatives.reserve(std::pow((polyord + 1), 2) * (polyord + 2) / 2);
      for (int idxi = 0; idxi < polyord + 1; ++idxi) {
         for (int idxj = 0; idxj < idxi + 1; ++idxj) {
            for (int idxk = 0; idxk < polyord + 1; ++idxk) {
               vec_derivatives.push_back(pleg.Derivative(idxi, q.X(idxk)) *
                                         pleg.Derivative(idxj, q.X(idxk)));
            };
         };
      };

      // lambda to find value of r^2 l_i'(x_k) l_j'(x_k), with x_k as the
      // Gaussian quadrature point
      auto funcnew = [this, &invnodespacingtwo](int idxelem, int idxi, int idxj,
                                                int idxk) {
         return invnodespacingtwo[idxelem] *
                std::pow(StandardIntervalMap(q.X(idxk), vec_noderadii[idxelem],
                                             vec_noderadii[idxelem + 1]),
                         2.0) *
                this->vec_derivatives[this->idxderiv(idxi, idxj, idxk)];
      };

      // reserve space
      vec_specelem.resize(lmax + 1);
      // vec_chol_solver.reserve(lmax + 1);

      // go through vector and fill out values
      for (int idxl = 0; idxl < lmax + 1; ++idxl) {
         std::vector<T> tripletlist;
         // tripletlist.reserve(nelem * (polyord + 1) * (polyord + 1) + 1);
         tripletlist.reserve(nelem * (((polyord + 1) * (polyord + 2)) / 2 - 1) +
                             1);

         {

            FLOAT zetasq = static_cast<FLOAT>(idxl) * static_cast<FLOAT>(idxl);

            // finding the integrals using the inner product and the vector
            // with the derivative values at the interpolation points
            for (int idxelem = 0; idxelem < nelem; ++idxelem) {
               FLOAT tmp;
               for (int idxi = 0; idxi < polyord + 1; ++idxi) {
                  for (int idxj = 0; idxj < idxi + 1; ++idxj) {
                     tmp = 0.0;
                     for (int idxk = 0; idxk < polyord + 1; ++idxk) {
                        tmp -= funcnew(idxelem, idxi, idxj, idxk) * q.W(idxk);
                     }

                     // push back into triplet list
                     if (idxj < idxi) {
                        tripletlist.push_back(T(idxelem * polyord + idxi,
                                                idxelem * polyord + idxj, tmp));
                     } else {
                        tripletlist.push_back(T(idxelem * polyord + idxi,
                                                idxelem * polyord + idxj,
                                                tmp + zetasq * q.W(idxi)));
                     }
                  };
               };
            };

            // final point
            tripletlist.push_back(
                T(matlen - 1, matlen - 1,
                  -static_cast<CFLOAT>(idxl + 1) *
                      static_cast<CFLOAT>(vec_noderadii[nelem])));

            // construct sparse matrix
            vec_specelem[idxl].resize(matlen, matlen);   // resize
            vec_specelem[idxl].setFromTriplets(tripletlist.begin(),
                                               tripletlist.end());   // set it
            vec_specelem[idxl].makeCompressed();                     // compress

            // Cholesky decomposition
            // EIGCHOL chol_tmp;
            // chol_tmp.compute(vec_specelem[idxl]);
            // vec_chol_solver[idxl].compute(vec_specelem[idxl]);
            // vec_chol_solver.push_back(EIGCHOL(vec_specelem[idxl]));
            // vec_chol_solver[idxl].compute(vec_specelem[idxl]);
            // vec_chol_solver.emplace_back(vec_specelem[idxl]);
            chol_solver.compute(vec_specelem[idxl]);
            tripletlist.clear();
         };
      }
   };

   SPHVEC solve(const int &lval, const SPHVEC &vec_force) {
      assert(m_isInitialized && "Not initialized");
      assert(((lval < lmax + 1) && (lval > -1)) && "Incorrect l");
      std::cout << "Rows: " << vec_force.rows()
                << ". Columns: " << vec_force.cols() << std::endl;
      std::cout << "Rows: " << vec_specelem[0].rows()
                << ". Columns: " << vec_specelem[0].cols() << std::endl;
      auto vecsol = chol_solver.solve(vec_force);

      return vecsol;
   };

   // Spherical_Integrator(sphericalmodel mymodel, ){};

 private:
   // SPARSEMAT specelem; // vector of spectral element matrices
   std::vector<SPARSEMAT> vec_specelem;
   int polyord, lmax, nlayer, nelem, matlen;
   std::vector<int> vec_nnodes;
   FLOAT rad_scale;
   // node values:
   std::vector<FLOAT> vec_noderadii;
   std::vector<FLOAT> vec_derivatives;

   // Gauss quadrature values
   GaussQuad::Quadrature1D<FLOAT> q;
   Eigen::SimplicialLDLT<SPARSEMAT, Eigen::Lower, Eigen::AMDOrdering<int>>
       chol_solver;
   // std::vector<Eigen::SimplicialLDLT<SPARSEMAT, Eigen::Lower,
   // Eigen::AMDOrdering<int>>> vec_chol_solver;
   // Eigen::CholmodSupernodalLLT<SPARSEMAT, Eigen::Lower> chol_solver;
   bool m_isInitialized;

   // return the index in the vector of derivatives for function l_i'(x_k)
   // l_j'(x_k). Stored such that contiguous in memory over k:
   int idxderiv(const int &i, const int &j, const int &k) {
      assert(m_isInitialized && "Not initialized");
      int n = this->polyord + 1;
      // return k + n * (j + i * (n - 1) - (i * (i - 1)) / 2);
      return k + n * (j + (i * (i + 1)) / 2);
   };
};   // end Poisson spherical harmonic class
//     template <typename FLOAT, template <typename, typename> class
//     sphericalmodel>
//         requires PlanetaryModel::SphericalElasticModel<sphericalmodel<FLOAT,
//         int>>
//     class Spherical_Integrator<FLOAT, sphericalmodel<FLOAT,
//     int>>::Spherical_Integrator(sphericalmodel<FLOAT, int> mymodel,
//     std::vector<int> nnodes, int ord_poly, int l_min, int l_max) :
//     vec_nnodes{nnodes}, polyord{ord_poly}, lmin{l_min}, lmax{l_max},
//     nlayer{mymodel.NumberOfLayers()}, rad_scale{mymodel.OuterRadius()} {
//                                                                                                                                                                                                                                                                                                   using namespace Interpolation;
//     // add first node at r = 0
//     vec_noderadii.push_back(0.0);
//     for (int idx = 0; idx < nlayer; ++idx)
//     {
//         double dx = (mymodel.UpperRadius(idx) - mymodel.LowerRadius(idx)) /
//         (static_cast<double>(vec_numnodes[idx] - 1) * rad_scale); for (int
//         idx2 = 1; idx2 < vec_numnodes[idx]; ++idx2)
//         {

//             vec_noderadii.push_back(mymodel.LowerRadius(idx) / rad_scale + dx
//             * static_cast<double>(idx2));
//         }
//     }
// };

// for (int idxi = 0; idxi < polyord + 1; ++idxi)
// {
//     for (int idxj = 0; idxj < idxi + 1; ++idxj)
//     {
//         for (int idxk = 0; idxk < polyord + 1; ++idxk)
//         {
//             auto mdiff = std::abs(vec_mat_lderiv[idxk](idxi, idxj) -
//             vec_derivatives[idxderiv(idxi, idxj, idxk, polyord)]); if (mdiff
//             > 0)
//             {
//                 std::cout << "Difference: " << mdiff << std::endl;
//             }
//         };
//     };
// };
}   // namespace GravityFunctions

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// still to be done properly with 3d tomography model
// if (vec_alldepths[idxelem * npoly + idxpoly] / 1000.0 >
//         tomo.GetDepths()[0] &&
//     vec_alldepths[idxelem * npoly + idxpoly] / 1000.0 <
//         tomo.GetDepths().back()) {
//    for (int idxt = 0; idxt < grid.NumberOfCoLatitudes(); ++idxt) {
//       for (int idxp = 0; idxp < grid.NumberOfLongitudes(); ++idxp) {
//          // vec_dvs[idxp + idxt * grid.NumberOfLongitudes()] +=
//          // tomo.GetValueAt(vec_alldepths[idxelem * npoly +
//          // idxpoly] / 1000.0,
//          // azimuthtolongitude(grid.LongitudeSpacing() *
//          // static_cast<double>(idxp)),
//          // polartolatitude(grid.CoLatitudes()[idxt])) * 0.3 /
//          // myprem.VS(laynum)(rphys(idxelem, q.X(idxpoly)));
//          // std::cout << tomo.GetValueAt(vec_alldepths[idxelem *
//          // npoly + idxpoly],
//          // azimuthtolongitude(grid.LongitudeSpacing() *
//          // static_cast<double>(idxp)),
//          // polartolatitude(grid.CoLatitudes()[idxt])) <<
//          // std::endl; vec_dvs[idxp + idxt *
//          // grid.NumberOfLongitudes()] += 1.0;
//          // std::cout << "Hello: \n";
//          auto idxuse = idxp + idxt * grid.NumberOfLongitudes();
//          vec_dvs[idxuse] *=
//              backgroundmodel.Density(laynum)(
//                  rphys(idxelem, q.X(idxpoly))) *
//              vec_j[idxelem * (npoly + 1) + idxpoly][idxuse];
//          // std::cout << rphys(idxelem, q.X(idxpoly)) <<
//          // std::endl;
//       }
//    }
//    // std::cout << std::endl;
//    // for (auto &didx : vec_dvs)
//    // {
//    //     std::cout << didx << std::endl;
//    // }
//    // std::cout << "Before FT: " << std::endl;
//    grid.ForwardTransformation(lMax, 0, vec_dvs, dvslm);
//    // std::cout << "After FT: " << std::endl;
// } else {

//    dvslm[0] =
//        backgroundmodel.Density(laynum)(rphys(idxelem, q.X(idxpoly)))
//        * sqrt(4.0 * 3.1415926535);
// }

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
#endif
