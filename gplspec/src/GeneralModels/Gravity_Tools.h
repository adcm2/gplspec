#ifndef GRAVITY_SOLVER_TOOLS_H
#define GRAVITY_SOLVER_TOOLS_H

#include "../Pseudospectral_Matrix_Wrapper.h"
#include "../Pseudospectral_Matrix_Wrapper3D.h"
#include "../Radial_Tools.h"
#include "../Shared_Utilities.h"
#include "Earth_General_Models_1D.h"
#include "../SphericalGeometryPreconditioner.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <GSHTrans/All>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <algorithm>
#include <cmath>
#include <limits>
#include "../Timer_Class.h"

namespace Gravity_Tools {

namespace detail {

struct SolveResult {
   Eigen::VectorXcd solution;
   Eigen::ComputationInfo status = Eigen::InvalidInput;
   Eigen::Index iterations = 0;
   double requested_tolerance = 0.0;
   double recursive_error = std::numeric_limits<double>::quiet_NaN();
   double residual_norm = std::numeric_limits<double>::quiet_NaN();
   double relative_residual = std::numeric_limits<double>::quiet_NaN();
};

template <class Model>
Eigen::SparseMatrix<std::complex<double>>
SphericalPreconditionerMatrix(Model &model, int lmax) {
   Eigen::SparseMatrix<std::complex<double>> matrix =
       -model.SpectralElementInformation().template fullmatrix<std::complex<double>>(
           model.Node_InformationP(), 0, lmax);
   matrix.makeCompressed();
   return matrix;
}

template <class Solver, class Operator>
class PreparedSolver {
 public:
   // The action is borrowed and must outlive this prepared solver.
   PreparedSolver(const Operator &op,
                  const Eigen::SparseMatrix<std::complex<double>> &preconditioner)
       : _operator(&op) {
      _solver.compute(op);
      _solver.preconditioner().addmatrix(preconditioner);
   }
   PreparedSolver(Operator &&,
                  const Eigen::SparseMatrix<std::complex<double>> &) = delete;
   PreparedSolver(const Operator &&,
                  const Eigen::SparseMatrix<std::complex<double>> &) = delete;

   SolveResult Solve(const Eigen::VectorXcd &rhs, double tolerance,
                     bool use_preconditioned_guess,
                     bool calculate_residual = true) {
      _solver.setTolerance(tolerance);
      SolveResult result;
      result.requested_tolerance = tolerance;
      if (use_preconditioned_guess) {
         const Eigen::VectorXcd guess = _solver.preconditioner().solve(rhs);
         result.solution = _solver.solveWithGuess(rhs, guess);
      } else {
         result.solution = _solver.solve(rhs);
      }
      result.status = _solver.info();
      result.iterations = _solver.iterations();
      result.recursive_error = _solver.error();
      if (calculate_residual) {
         const Eigen::VectorXcd residual = rhs - (*_operator * result.solution);
         result.residual_norm = residual.norm();
         const double rhs_norm = rhs.norm();
         result.relative_residual =
             rhs_norm == 0.0
                 ? (result.residual_norm == 0.0
                        ? 0.0
                        : std::numeric_limits<double>::infinity())
                 : result.residual_norm / rhs_norm;
      }
      return result;
   }

 private:
   const Operator *_operator;
   Solver _solver;
};

inline double RadialSourceWeight(const Radial_Tools::RadialMesh &mesh,
                                const GaussQuad::Quadrature1D<double> &q,
                                std::size_t element, std::size_t node) {
   const double radius = mesh.NodeRadius(element, node);
   return q.W(node) * 0.5 * mesh.ElementWidth(element) * radius * radius;
}

}  // namespace detail

// The spherical model has only the monopole density coefficient.
inline Eigen::VectorXcd
FindForce(GeneralEarthModels::spherical_1D &model, const int &lMax) {
   const auto &mesh = model.Node_InformationP();
   const auto q = model.q();
   const std::size_t stride = static_cast<std::size_t>(lMax + 1) * (lMax + 1);
   Eigen::VectorXcd source = Eigen::VectorXcd::Zero(
       mesh.NumberOfGlobalRadialNodes() * stride);

   constexpr double pi = 3.1415926535;
   const double monopole_scale =
       4.0 * pi * model.GravitationalConstant() * std::sqrt(4.0 * pi);
   for (std::size_t element = 0; element < mesh.NumberOfElements(); ++element) {
      for (std::size_t node = 0; node < mesh.NumNodesWithinElement(); ++node) {
         const auto radial = mesh.GlobalNodeIndex(element, node);
         source(radial * stride) +=
             monopole_scale * detail::RadialSourceWeight(mesh, q, element, node) *
             model.Density(element, node);
      }
   }
   return source;
}

// The mapped model transforms each local density sample to scalar harmonics.
inline Eigen::VectorXcd
FindForce(GeneralEarthModels::Density3D &model) {
   using Complex = std::complex<double>;
   const auto &mesh = model.Node_InformationP();
   const auto &grid = model.GSH_Grid();
   const auto &q = model.qP();
   const int lMax = grid.MaxDegree();
   const std::size_t stride = static_cast<std::size_t>(lMax + 1) * (lMax + 1);
   const std::size_t radial_nodes = mesh.NumberOfGlobalRadialNodes();
   Eigen::VectorXcd source = Eigen::VectorXcd::Zero(radial_nodes * stride);

   const auto nonnegative_count =
       GSHTrans::GSHIndices<GSHTrans::NonNegative>(lMax, lMax, 0).Size();
   std::vector<Complex> nonnegative(nonnegative_count);
   std::vector<Complex> full(stride);
   constexpr double pi = 3.1415926535897932;
   const double source_scale = 4.0 * pi * model.GravitationalConstant();

   for (std::size_t element = 0; element < mesh.NumberOfElements(); ++element) {
      for (std::size_t node = 0; node < mesh.NumNodesWithinElement(); ++node) {
         std::fill(nonnegative.begin(), nonnegative.end(), Complex{});
         grid.ForwardTransformation(lMax, 0,
                                    model.DensityAtRadialNode(element, node),
                                    nonnegative);
         GPLSpec::detail::ExpandRealScalarCoefficients(lMax, nonnegative, full);
         const double weight = source_scale *
                               detail::RadialSourceWeight(mesh, q, element, node);
         const Eigen::Index base = mesh.GlobalNodeIndex(element, node) * stride;
         for (std::size_t harmonic = 0; harmonic < stride; ++harmonic) {
            source(base + harmonic) += weight * full[harmonic];
         }
      }
   }

   // Preserve the source domain: non-monopole center terms and both highest-
   // degree endpoint orders are zero, even though the operator has ambient
   // output slots for them.
   for (int degree = 1; degree <= lMax; ++degree) {
      for (int order = -degree; order <= degree; ++order) {
         source(degree * degree + degree + order) = 0.0;
      }
   }
   const std::size_t negative_endpoint = static_cast<std::size_t>(lMax * lMax);
   const std::size_t positive_endpoint = negative_endpoint + 2 * lMax;
   for (std::size_t radial = 0; radial < radial_nodes; ++radial) {
      const Eigen::Index base = radial * stride;
      source(base + negative_endpoint) = 0.0;
      source(base + positive_endpoint) = 0.0;
   }
   return source;
}

// function for finding the force vector:
template <class mapclass>
   requires PlanetaryModel::RadialMappingClass<mapclass>
inline Eigen::VectorXcd
FindBoundaryPerturbationForce(GeneralEarthModels::Density3D &inp_model,
                              mapclass &inp_map) {
   using Real = double;
   using Complex = std::complex<Real>;
   std::size_t nelem = inp_model.Num_Elements();
   int npoly = inp_model.q().N() - 1;
   int matlen = nelem * npoly + 1;   // size of matrix
   int lMax = inp_model.GSH_Grid().MaxDegree();

   auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).Size();
   auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).Size();
   auto rphys = [&inp_model](int idxelem, double x) {
      return Radial_Tools::StandardIntervalMap(
          x, inp_model.Node_Information().ElementLowerRadius(idxelem),
          inp_model.Node_Information().ElementUpperRadius(idxelem));
   };

   //    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _mat_gaussderiv;
   // finding \rho \mathbf{s} in canonical components spatially
   // firstly we need to have the background model density and the perturbation
   // from the perturbing model:
   // go through all layers, find what s is, as we currently just have radial,
   // that should be fairly simple
   int totnum = inp_model.GSH_Grid().NumberOfCoLatitudes() *
                inp_model.GSH_Grid().NumberOfLongitudes();

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
   auto elemwidth = [&inp_model](int idxelem) {
      return inp_model.Node_Information().ElementWidth(idxelem);
   };
   for (int idx = 0; idx < inp_model.q().Points().size(); ++idx) {
      gaussquadpoints.push_back(inp_model.q().X(idx));
      gaussquadweights.push_back(inp_model.q().W(idx));
   }
   for (int idxelem = 0; idxelem < nelem; ++idxelem) {
      double inv2 = 2.0 / elemwidth(idxelem);
      // if (laynum < backgroundmodel.NumberOfLayers()) {
      //    if (!(vec_noderadii[idxelem] < backgroundmodel.UpperRadius(laynum)))
      //    {
      //       laynum += 1;
      //       std::cout << "laynum: " << laynum << ", idxelem: " << idxelem
      //                 << std::endl;
      //    };
      // }

      int laynum = inp_model.Node_Information().LayerNumber(idxelem);
      // looping over quadrature nodes
      int idxpolymax = npoly + 1;
      for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {

         // auto rval = GravityFunctions::StandardIntervalMap(
         //     gaussquadpoints[idxpoly],
         //     inp_model.Node_Information().ElementLowerRadius(idxelem),
         //     inp_model.Node_Information().ElementUpperRadius(idxelem));
         auto rval = inp_model.Node_Information().NodeRadius(idxelem, idxpoly);
         auto invr = 1 / rval;
         auto rval2 = rval * rval;

         /////////////////////////////////////////////////////////////

         // finding spatial q0, qp, qm
         std::vector<Complex> spatial_q0(totnum, 0.0), spatial_qm(totnum, 0.0),
             spatial_qp(totnum, 0.0);

         // auto raddensity = inp_model.Density;
         {
            std::size_t idxspatial = 0;
            std::size_t totnum = inp_model.Node_Information().NumberOfLayers();
            double ballrad = inp_model.Node_Information().OuterRadius();
            double planetrad = inp_model.Node_Information().PlanetRadius();

            // loop through spatial indices
            for (auto it : inp_model.GSH_Grid().CoLatitudes()) {
               for (auto ip : inp_model.GSH_Grid().Longitudes()) {
                  auto raddensity =
                      inp_model.Density_Point(idxelem, idxpoly, idxspatial);
                  if (laynum < totnum - 1) {
                     spatial_q0[idxspatial] =
                         inp_map.RadialMapping(laynum)(rval, it, ip) *
                         raddensity;
                  } else {
                     spatial_q0[idxspatial] =
                         inp_map.RadialMapping(laynum)(planetrad, it, ip) *
                         (ballrad - rval) / (ballrad - planetrad) * raddensity;
                  }
                  ++idxspatial;
               }
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
            inp_model.GSH_Grid().ForwardTransformation(lMax, -1, spatial_qm,
                                                       gsph_qm1);
            inp_model.GSH_Grid().ForwardTransformation(lMax, 0, spatial_q0,
                                                       gsph_q0);
            inp_model.GSH_Grid().ForwardTransformation(lMax, +1, spatial_qp,
                                                       gsph_qp1);
         }

         // step 5: evaluating radial integrals to give Ax

         // 0 component term
         {
            // multiplication factor and index
            std::size_t idx2 = myidxfunc(idxelem, 0, 0, 0);
            // std::size_t idx2 = idxelem * npoly;
            auto mult1 = gaussquadweights[idxpoly] * rval2;

            // loop over radii
            for (int idxn = 0; idxn < npoly + 1; ++idxn) {

               // index and multiplication factor
               std::size_t idx3 = 0;
               auto multfact = inp_model.GaussDerivative(idxn, idxpoly) * mult1;

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
   // const double bigg_db = 6.6743 * std::pow(10.0, -11.0);
   const double pi_db = 3.1415926535;
   double multcomp = 4.0 * pi_db * inp_model.GravitationalConstant();
   // std::size_t idxmaxval = std::pow(lMax + 1, 2) * matlen;
   for (auto idx1 = 0; idx1 < matlen * _size0; ++idx1) {
      vec_output2(idx1) *= multcomp;
   }

   return vec_output2;
};

// function for finding the result of advecting the solution of a problem with
// the mapping defined by inp_map
template <class mapclass>
   requires PlanetaryModel::RadialMappingClass<mapclass>
inline Eigen::VectorXcd
AdvectiveBoundaryPerturbation(GeneralEarthModels::Density3D &inp_model,
                              mapclass &inp_map,
                              const Eigen::VectorXcd &vec_phi) {
   using Real = double;
   using Complex = std::complex<Real>;
   std::size_t nelem = inp_model.Num_Elements();
   int npoly = inp_model.q().N() - 1;
   int matlen = nelem * npoly + 1;   // size of matrix
   int lMax = inp_model.GSH_Grid().MaxDegree();

   auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).Size();
   auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).Size();
   auto rphys = [&inp_model](int idxelem, double x) {
      return Radial_Tools::StandardIntervalMap(
          x, inp_model.Node_Information().ElementLowerRadius(idxelem),
          inp_model.Node_Information().ElementUpperRadius(idxelem));
   };

   //    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _mat_gaussderiv;
   // finding \rho \mathbf{s} in canonical components spatially
   // firstly we need to have the background model density and the perturbation
   // from the perturbing model:
   // go through all layers, find what s is, as we currently just have radial,
   // that should be fairly simple
   int totnum = inp_model.GSH_Grid().NumberOfCoLatitudes() *
                inp_model.GSH_Grid().NumberOfLongitudes();

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
   auto elemwidth = [&inp_model](int idxelem) {
      return inp_model.Node_Information().ElementWidth(idxelem);
   };
   for (int idx = 0; idx < inp_model.q().Points().size(); ++idx) {
      gaussquadpoints.push_back(inp_model.q().X(idx));
      gaussquadweights.push_back(inp_model.q().W(idx));
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
      std::size_t laynum = inp_model.Node_Information().LayerNumber(idxelem);
      for (int idxpoly = idxpolymin; idxpoly < idxpolymax; ++idxpoly) {
         // radius, inverse and squared
         auto rval = inp_model.Node_Information().NodeRadius(idxelem, idxpoly);
         auto invr = 1 / rval;
         auto rval2 = rval * rval;

         ///////////////////////////////////////
         // step 1: finding nabla zeta
         std::vector<std::complex<double>> gsph_nz0(_size0, 0.0),
             gsph_nzp1(_sizepm, 0.0),
             gsph_nzm1(_sizepm, 0.0);   //                  declaration

         // looping over l,m values
         {

            // cache friendly trial
            //  getting first index
            std::size_t idx1 = myidxfunc(idxelem, 0, 0, 0);

            // looping over radii
            for (int idxn = 0; idxn < npoly + 1; ++idxn) {
               std::size_t idx2 = 0;
               auto multfact = inp_model.GaussDerivative(idxn, idxpoly);

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
         //       // std::cout << "invr: " << invr << std::endl;
         //       // }

         ///////////////////////////////////////////////////////////////////////////
         // step 2: transforming (\nabla \zeta) into spatial domain
         FFTWpp::vector<Complex> spatial_nzm1(totnum, 0.0),
             spatial_nzp1(totnum, 0.0),
             spatial_nz0(totnum, 0.0);   // declaration

         // GSPH transforms:
         // lhs.gridused().InverseTransformation(lhs.lMax(), 0, gsph_nz0,
         //                                      spatial_nz0);   // 0th order
         inp_model.GSH_Grid().InverseTransformation(lMax, 0, gsph_nz0,
                                                    spatial_nz0);   // 0th order
         inp_model.GSH_Grid().InverseTransformation(lMax, 1, gsph_nzp1,
                                                    spatial_nzp1);   //+1 order
         inp_model.GSH_Grid().InverseTransformation(lMax, -1, gsph_nzm1,
                                                    spatial_nzm1);   //-1 order

         ///////////////////////////////////////////////////////////////////////////
         // finding s in canonical components spatially
         std::vector<Complex> spatial_q0(totnum, 0.0), spatial_qm(totnum, 0.0),
             spatial_qp(totnum, 0.0);

         {
            std::size_t idxspatial = 0;
            std::size_t totlaynum =
                inp_model.Node_Information().NumberOfLayers();
            double ballrad = inp_model.Node_Information().OuterRadius();
            double planetrad = inp_model.Node_Information().PlanetRadius();

            // loop through spatial indices
            for (auto it : inp_model.GSH_Grid().CoLatitudes()) {
               for (auto ip : inp_model.GSH_Grid().Longitudes()) {
                  // auto idxuse = idxt * grid.NumberOfLongitudes() + idxp;
                  // double radscalval;
                  if (laynum < totlaynum - 1) {
                     spatial_q0[idxspatial] =
                         inp_map.RadialMapping(laynum)(rval, it, ip);
                  } else {
                     spatial_q0[idxspatial] =
                         inp_map.RadialMapping(laynum)(planetrad, it, ip) *
                         (ballrad - rval) / (ballrad - planetrad);
                  }
                  ++idxspatial;
               }
            }
         }

         ///////////////////////////////////////////////////////////////////////////
         // finding the multiple of the two
         std::vector<Complex> spatial_multtot(totnum, 0.0);
         // for (int idxt = 0; idxt < grid.NumberOfCoLatitudes(); ++idxt) {
         //    for (int idxp = 0; idxp < grid.NumberOfLongitudes(); ++idxp) {
         {
            std::size_t idxuse = 0;
            for (auto it : inp_model.GSH_Grid().CoLatitudes()) {
               for (auto ip : inp_model.GSH_Grid().Longitudes()) {
                  // auto idxuse = idxt * grid.NumberOfLongitudes() + idxp;
                  spatial_multtot[idxuse] -=
                      spatial_qm[idxuse] * spatial_nzp1[idxuse];
                  spatial_multtot[idxuse] +=
                      spatial_q0[idxuse] * spatial_nz0[idxuse];
                  spatial_multtot[idxuse] -=
                      spatial_qp[idxuse] * spatial_nzm1[idxuse];
                  ++idxuse;
               }
            }
         }

         ///////////////////////////////////////////////////////////////////////////
         // converting back into spherical harmonics
         std::vector<std::complex<double>> gsph_multtot(_size0, 0.0);
         inp_model.GSH_Grid().ForwardTransformation(lMax, 0, spatial_multtot,
                                                    gsph_multtot);
         std::size_t idx1 = myidxfunc(idxelem, idxpoly, 0, 0);
         std::size_t idx3 = 0;
         for (int idxl = 0; idxl < lMax + 1; ++idxl) {
            // loop over m
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               vec_output2(idx1) = gsph_multtot[idx3];   // global

               // increment indices
               ++idx1;
               ++idx3;
            }
         }
      }
   }

   return vec_output2;
};

template <class Grid>
inline detail::SolveResult
SolveGravitationalPotential(GeneralEarthModels::spherical_1D &inp_model,
                            Grid &grid, double relerr = std::pow(10.0, -6.0)) {
   using Complex = std::complex<double>;
   using BICGSTAB =
       Eigen::BiCGSTAB<MatrixReplacement<std::complex<double>>,
                       Eigen::SphericalGeometryPreconditioner<Complex>>;
   MatrixReplacement<std::complex<double>> mymatrix(inp_model, grid);
   Eigen::VectorXcd vec_fullforce =
       Gravity_Tools::FindForce(inp_model, grid.MaxDegree());
   const auto preconditioner =
       detail::SphericalPreconditionerMatrix(inp_model, grid.MaxDegree());
   detail::PreparedSolver<BICGSTAB, decltype(mymatrix)> solver(mymatrix,
                                                               preconditioner);
   return solver.Solve(vec_fullforce, relerr, false);
}

template <class Grid>
inline Eigen::VectorXcd
FindGravitationalPotential(GeneralEarthModels::spherical_1D &inp_model,
                           Grid &grid, double relerr = std::pow(10.0, -6.0)) {
   return SolveGravitationalPotential(inp_model, grid, relerr).solution;
};

// 3D solver
inline detail::SolveResult
SolveGravitationalPotential(GeneralEarthModels::Density3D &inp_model,
                            double relerr = std::pow(10.0, -6.0)) {
   using Complex = std::complex<double>;
   using CONJG = Eigen::ConjugateGradient<
       MatrixReplacement3D<Complex>, Eigen::Lower | Eigen::Upper,
       Eigen::SphericalGeometryPreconditioner<Complex>>;
   MatrixReplacement3D<std::complex<double>> mymatrix(inp_model);
   Eigen::VectorXcd vec_fullforce = Gravity_Tools::FindForce(inp_model);
   const auto preconditioner = detail::SphericalPreconditionerMatrix(
       inp_model, inp_model.GSH_Grid().MaxDegree());
   detail::PreparedSolver<CONJG, decltype(mymatrix)> solver(mymatrix,
                                                            preconditioner);
   return solver.Solve(vec_fullforce, relerr, true);
}

inline auto FindGravitationalPotential(GeneralEarthModels::Density3D &inp_model,
                                      double relerr = std::pow(10.0, -6.0)) {
   const auto result = SolveGravitationalPotential(inp_model, relerr);
   return inp_model.SingleEigenVectorToGeneralFormat(result.solution);
};

// 3D solver
inline auto
SphericalHarmonicSensitivityKernel(GeneralEarthModels::Density3D &inp_model,
                                   int l, int m,
                                   double relerr = std::pow(10.0, -6.0)) {
   using Complex = std::complex<double>;
   using CONJG = Eigen::ConjugateGradient<
       MatrixReplacement3D<Complex>, Eigen::Lower | Eigen::Upper,
       Eigen::SphericalGeometryPreconditioner<Complex>>;

   MatrixReplacement3D<std::complex<double>> mymatrix(inp_model);

   // Eigen::VectorXcd vec_fullforce = Gravity_Tools::FindForce(inp_model);
   // force vector is simply a minus one in the lm component on the boundary. We
   // simply need to find its index.
   int npoly = inp_model.q().N() - 1;
   int lMax = inp_model.GSH_Grid().MaxDegree();
   int nelem = inp_model.Node_InformationP().NumberOfElements();
   int matlen = nelem * npoly + 1;   // size of matrix
   Eigen::VectorXcd vec_fullforce =
       Eigen::VectorXcd::Zero(std::pow(lMax + 1, 2) * matlen);
   std::size_t lmidx = (npoly + npoly * (nelem - 1)) * std::pow(lMax + 1, 2);
   lmidx += l * l + l - m;
   double pi_db = 3.1415926535;
   double ballrad = inp_model.Node_InformationP().OuterRadius();
   vec_fullforce(lmidx) =
       4.0 * pi_db * inp_model.GravitationalConstant() * std::pow(-1.0, m);

   const auto preconditioner = detail::SphericalPreconditionerMatrix(
       inp_model, inp_model.GSH_Grid().MaxDegree());
   detail::PreparedSolver<CONJG, decltype(mymatrix)> solver(mymatrix,
                                                            preconditioner);
   const auto result = solver.Solve(vec_fullforce, relerr, true, false);
   return inp_model.SingleEigenVectorToGeneralFormat(result.solution);
};

// 3D solver
inline auto
SphericalHarmonicSensitivityKernel(GeneralEarthModels::Density3D &inp_model,
                                   std::vector<int> l, std::vector<int> m,
                                   std::vector<double> multval,
                                   double relerr = std::pow(10.0, -6.0)) {
   assert(l.size() == m.size() && "Need same size for l and m");
   assert(multval.size() == m.size() && "Need same size for l and m");
   using Complex = std::complex<double>;
   using BICGSTAB =
       Eigen::BiCGSTAB<MatrixReplacement3D<std::complex<double>>,
                       Eigen::SphericalGeometryPreconditioner<Complex>>;
   MatrixReplacement3D<std::complex<double>> mymatrix(inp_model);

   // Eigen::VectorXcd vec_fullforce = Gravity_Tools::FindForce(inp_model);
   // force vector is simply a minus one in the lm component on the boundary. We
   // simply need to find its index.
   int npoly = inp_model.q().N() - 1;
   int lMax = inp_model.GSH_Grid().MaxDegree();
   int nelem = inp_model.Node_Information().NumberOfElements();
   int matlen = nelem * npoly + 1;   // size of matrix
   Eigen::VectorXcd vec_fullforce =
       Eigen::VectorXcd::Zero(std::pow(lMax + 1, 2) * matlen);
   std::size_t lmidx = (npoly + npoly * (nelem - 1)) * std::pow(lMax + 1, 2);
   for (int i = 0; i < l.size(); ++i) {
      auto lidx = lmidx + l[i] * l[i] + l[i] - m[i];
      double pi_db = 3.1415926535;
      double ballrad = inp_model.Node_Information().OuterRadius();
      vec_fullforce(lidx) = 4.0 * pi_db * inp_model.GravitationalConstant() *
                            std::pow(-1.0, m[i]) * multval[i];
   }

   const auto preconditioner = detail::SphericalPreconditionerMatrix(
       inp_model, inp_model.GSH_Grid().MaxDegree());
   detail::PreparedSolver<BICGSTAB, decltype(mymatrix)> solver(mymatrix,
                                                               preconditioner);
   const auto result = solver.Solve(vec_fullforce, relerr, true, false);
   return inp_model.SingleEigenVectorToGeneralFormat(result.solution);
};

// 3D solver
inline auto
FindGravitationalPotentialPerturbation(
    GeneralEarthModels::Density3D &inp_model,
    GeneralEarthModels::MappingPerturbation &inp_map,
    double relerr1 = std::pow(10.0, -6.0),
    double relerr2 = std::pow(10.0, -2.0)) {
   using Complex = std::complex<double>;
   using BICGSTAB =
       Eigen::BiCGSTAB<MatrixReplacement3D<std::complex<double>>,
                       Eigen::SphericalGeometryPreconditioner<Complex>>;
   using CONJG = Eigen::ConjugateGradient<
       MatrixReplacement3D<Complex>, Eigen::Lower | Eigen::Upper,
       Eigen::SphericalGeometryPreconditioner<Complex>>;
   Timer _timer;
   // pseudospectral matrix replacement
   MatrixReplacement3D<std::complex<double>> mymatrix(inp_model);
   MatrixReplacement3D<std::complex<double>> matrixpert(
       inp_model, inp_map, false);

   // get force vector
   Eigen::VectorXcd vec_fullforce = Gravity_Tools::FindForce(inp_model);
   // Eigen::VectorXcd vec_ff2 = vec_fullforce * 0.0;
   // std::cout << vec_ff2.size() << "\n";
   // Eigen::VectorXcd vec_ff2 =
   //     Eigen::VectorXcd::Constant(vec_fullforce.size(), 1.0);

   const auto preconditioner = detail::SphericalPreconditionerMatrix(
       inp_model, inp_model.GSH_Grid().MaxDegree());
   detail::PreparedSolver<CONJG, decltype(mymatrix)> solver(mymatrix,
                                                            preconditioner);
   _timer.start();
   const auto base_result = solver.Solve(vec_fullforce, relerr1, true, false);
   _timer.stop("Time for base solution");
   Eigen::VectorXcd vec_pertforce = -(matrixpert * base_result.solution);

   _timer.start();
   const auto perturbation_result =
       solver.Solve(vec_pertforce, relerr2, true, false);
   _timer.stop("Time for perturbation solution");
   Eigen::VectorXcd fullsol =
       base_result.solution + perturbation_result.solution;
   auto retsol = inp_model.SingleEigenVectorToGeneralFormat(fullsol);
   return retsol;
};

// 3D solver
template <class mapclass>
   requires PlanetaryModel::RadialMappingClass<mapclass>
inline auto
FindGravitationalPotentialClassicalPerturbation(
    GeneralEarthModels::Density3D &inp_model, mapclass &inp_map,
    double relerr = std::pow(10.0, -6.0)) {
   using Complex = std::complex<double>;
   using BICGSTAB =
       Eigen::BiCGSTAB<MatrixReplacement3D<std::complex<double>>,
                       Eigen::SphericalGeometryPreconditioner<Complex>>;

   MatrixReplacement3D<std::complex<double>> mymatrix(inp_model);

   // get force vector
   Eigen::VectorXcd vec_fullforce = Gravity_Tools::FindForce(inp_model);
   // Eigen::VectorXcd vec_ff2 = vec_fullforce * 0.0;
   // std::cout << vec_ff2.size() << "\n";
   // Eigen::VectorXcd vec_ff2 =
   //     Eigen::VectorXcd::Constant(vec_fullforce.size(), 1.0);

   const auto preconditioner = detail::SphericalPreconditionerMatrix(
       inp_model, inp_model.GSH_Grid().MaxDegree());
   detail::PreparedSolver<BICGSTAB, decltype(mymatrix)> solver(mymatrix,
                                                               preconditioner);
   const auto base_result = solver.Solve(vec_fullforce, relerr, true, false);

   Eigen::VectorXcd testforce =
       FindBoundaryPerturbationForce(inp_model, inp_map);
   const auto perturbation_result = solver.Solve(testforce, relerr, true, false);
   Eigen::VectorXcd advectpert =
       AdvectiveBoundaryPerturbation(inp_model, inp_map, base_result.solution);

   Eigen::VectorXcd fullsol = base_result.solution + advectpert -
                              perturbation_result.solution;
   auto retsol = inp_model.SingleEigenVectorToGeneralFormat(fullsol);
   return retsol;
};

// The normalized radial Green function is evaluated as two cumulative moments.
// Scaling each moment at its current boundary keeps the powers bounded by one.
inline auto
GravitationalSphericalIntegral(const GeneralEarthModels::Density3D &model) {
   using Complex = std::complex<double>;
   constexpr double pi = 3.141592653589793238462643383279502884;

   const auto &mesh = model.Node_Information();
   const auto &quadrature = model.qP();
   const auto &grid = model.GSH_Grid();
   const std::size_t elements = mesh.NumberOfElements();
   const std::size_t nodes_per_element = mesh.NumNodesWithinElement();
   const int lmax = grid.MaxDegree();
   const std::size_t stride = static_cast<std::size_t>(lmax + 1) * (lmax + 1);
   const std::size_t spatial_size = grid.NumberOfLongitudes() *
                                    grid.NumberOfCoLatitudes();

   // Transform each element-local density sample once. Shared interface nodes
   // intentionally occur in both adjacent elements with their own material.
   std::vector<Complex> density_coefficients(elements * nodes_per_element * stride);
   std::vector<Complex> density_on_grid(spatial_size);
   std::vector<Complex> coefficients(stride);
   for (std::size_t element = 0; element < elements; ++element) {
      for (std::size_t node = 0; node < nodes_per_element; ++node) {
         for (std::size_t spatial = 0; spatial < spatial_size; ++spatial) {
            density_on_grid[spatial] =
                Complex{model.Density_Point(element, node, spatial), 0.0};
         }
         std::fill(coefficients.begin(), coefficients.end(), Complex{});
         grid.ForwardTransformation(lmax, 0, density_on_grid, coefficients);
         const std::size_t sample = (element * nodes_per_element + node) * stride;
         std::copy(coefficients.begin(), coefficients.end(),
                   density_coefficients.begin() + sample);
      }
   }

   std::vector<std::vector<Complex>> potential(
       elements + 1, std::vector<Complex>(stride, Complex{}));
   std::vector<Complex> inward(elements + 1, Complex{});
   std::vector<Complex> outward(elements + 1, Complex{});

   for (int degree = 0; degree <= lmax; ++degree) {
      for (int order = -degree; order <= degree; ++order) {
         const std::size_t harmonic =
             static_cast<std::size_t>(degree * degree + degree + order);
         std::fill(inward.begin(), inward.end(), Complex{});
         std::fill(outward.begin(), outward.end(), Complex{});

         // inward[e] = r_e^(-l-1) * integral_0^r_e rho_lm(s) s^(l+2) ds.
         for (std::size_t element = 0; element < elements; ++element) {
            const double lower = mesh.ElementLowerRadius(element);
            const double upper = mesh.ElementUpperRadius(element);
            const double ratio = lower / upper;
            Complex segment{};
            for (std::size_t node = 0; node < nodes_per_element; ++node) {
               const std::size_t sample =
                   (element * nodes_per_element + node) * stride + harmonic;
               const double radius = mesh.NodeRadius(element, node);
               segment += quadrature.W(node) * density_coefficients[sample] *
                          std::pow(radius / upper, degree + 2);
            }
            inward[element + 1] =
                inward[element] * std::pow(ratio, degree + 1) +
                (0.5 * upper * (upper - lower)) * segment;
         }

         // outward[e] = r_e^l * integral_r_e^R rho_lm(s) s^(1-l) ds.
         // Only the monopole has a nonzero contribution at the origin.
         for (std::size_t element = elements; element-- > 1;) {
            const double lower = mesh.ElementLowerRadius(element);
            const double upper = mesh.ElementUpperRadius(element);
            Complex segment{};
            for (std::size_t node = 0; node < nodes_per_element; ++node) {
               const std::size_t sample =
                   (element * nodes_per_element + node) * stride + harmonic;
               const double scaled_radius =
                   mesh.NodeRadius(element, node) / lower;
               segment += quadrature.W(node) * density_coefficients[sample] *
                          std::pow(scaled_radius, 1 - degree);
            }
            segment *= 0.5 * lower *
                       (upper - lower);
            const double ratio = lower / upper;
            outward[element] =
                segment + outward[element + 1] * std::pow(ratio, degree);
         }

         if (degree == 0) {
            const double lower = mesh.ElementLowerRadius(0);
            const double upper = mesh.ElementUpperRadius(0);
            Complex center_segment{};
            for (std::size_t node = 0; node < nodes_per_element; ++node) {
               const std::size_t sample = harmonic + node * stride;
               center_segment += quadrature.W(node) *
                                 density_coefficients[sample] *
                                 mesh.NodeRadius(0, node);
            }
            outward[0] = outward[1] +
                         0.5 * (upper - lower) * center_segment;
         }

         const double green_scale = -4.0 * pi * model.GravitationalConstant() /
                                    (2.0 * degree + 1.0);
         for (std::size_t radial = 0; radial <= elements; ++radial) {
            potential[radial][harmonic] =
                green_scale * (inward[radial] + outward[radial]);
            if (lmax > 0 && degree == lmax && std::abs(order) == degree) {
               // Preserve both reserved highest-order coefficients.
               potential[radial][harmonic] = Complex{};
            }
         }
      }
   }

   return potential;
}

namespace detail {
inline double HomogeneousSphereScalarPotential(
    const GeneralEarthModels::Density3D &model, double radius) {
   constexpr double pi = 3.141592653589793238462643383279502884;
   const double body_radius = model.Node_Information().PlanetRadius();
   const double density = model.Density_Point(0, 0, 0);
   const double scale = 2.0 * pi * model.GravitationalConstant() * density / 3.0;
   if (radius <= body_radius) {
      return scale * (radius * radius - 3.0 * body_radius * body_radius);
   }
   return -2.0 * scale * std::pow(body_radius, 3.0) / radius;
}
}  // namespace detail

// Mesh-boundary values are the spherical-harmonic Y00 coefficient.
inline auto
HomogeneousSphereIntegral(const GeneralEarthModels::Density3D &model) {
   constexpr double sqrt_four_pi = 3.544907701811032054596334966682290366;
   const auto &mesh = model.Node_Information();
   std::vector<double> coefficient(mesh.NumberOfElements() + 1);
   for (std::size_t element = 0; element < mesh.NumberOfElements(); ++element) {
      coefficient[element] =
          sqrt_four_pi * detail::HomogeneousSphereScalarPotential(
                             model, mesh.ElementLowerRadius(element));
   }
   coefficient.back() =
       sqrt_four_pi * detail::HomogeneousSphereScalarPotential(
                          model, mesh.OuterRadius());
   return coefficient;
}

// Arbitrary-radius values are scalar potentials, as used by spatial outputs.
inline auto
HomogeneousSphereIntegral(const GeneralEarthModels::Density3D &model,
                          const std::vector<double> &radii) {
   std::vector<double> potential(radii.size());
   for (std::size_t index = 0; index < radii.size(); ++index) {
      potential[index] =
          detail::HomogeneousSphereScalarPotential(model, radii[index]);
   }
   return potential;
}

}   // namespace Gravity_Tools

#endif
