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
#include <stdexcept>
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

inline void ValidateSensitivityTolerance(double tolerance) {
   if (!std::isfinite(tolerance) || tolerance <= 0.0) {
      throw std::invalid_argument(
          "sensitivity solver tolerance must be finite and positive");
   }
}

inline Eigen::VectorXcd
SensitivitySelectorRhs(const GeneralEarthModels::Density3D &model,
                       const std::vector<int> &degrees,
                       const std::vector<int> &orders,
                       const std::vector<double> &multipliers) {
   if (degrees.size() != orders.size() || degrees.size() != multipliers.size()) {
      throw std::invalid_argument(
          "sensitivity selector degrees, orders, and weights must have equal lengths");
   }

   const int lmax = model.GSH_Grid().MaxDegree();
   const std::size_t stride =
       static_cast<std::size_t>(lmax + 1) * (lmax + 1);
   const auto &mesh = model.Node_InformationP();
   const std::size_t radial_count = mesh.NumberOfGlobalRadialNodes();
   if (radial_count == 0 ||
       radial_count > static_cast<std::size_t>(
                          std::numeric_limits<Eigen::Index>::max()) /
                          stride) {
      throw std::length_error("sensitivity RHS dimensions are not representable");
   }

   Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(
       static_cast<Eigen::Index>(radial_count * stride));
   const std::size_t outer_index = (radial_count - 1) * stride;
   constexpr double legacy_pi = 3.1415926535;
   const double prefactor = 4.0 * legacy_pi * model.GravitationalConstant();
   for (std::size_t i = 0; i < degrees.size(); ++i) {
      const int degree = degrees[i];
      const int order = orders[i];
      if (degree < 0 || degree > lmax || order < -degree || order > degree) {
         throw std::out_of_range("sensitivity harmonic selector is out of range");
      }
      if (!std::isfinite(multipliers[i])) {
         throw std::invalid_argument("sensitivity selector weight must be finite");
      }
      const std::size_t harmonic = static_cast<std::size_t>(
          degree * degree + degree - order);
      const double phase = order % 2 == 0 ? 1.0 : -1.0;
      // Assignment preserves the established last-value-wins behavior for duplicates.
      rhs[static_cast<Eigen::Index>(outer_index + harmonic)] =
          prefactor * phase * multipliers[i];
   }
   return rhs;
}

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
FindBoundaryPerturbationForce(GeneralEarthModels::Density3D &model,
                              mapclass &mapping) {
   using Complex = std::complex<double>;
   const auto &mesh = model.Node_Information();
   const auto &grid = model.GSH_Grid();
   const std::size_t element_count = model.Num_Elements();
   const int order = model.q().N() - 1;
   const int degree = grid.MaxDegree();
   const std::size_t stride =
       GSHTrans::GSHIndices<GSHTrans::All>(degree, degree, 0).Size();
   const std::size_t node_count = element_count * order + 1;
   constexpr double pi = 3.141592653589793238462643383279502884;
   const double scale = -4.0 * pi * model.GravitationalConstant();
   Eigen::VectorXcd force = Eigen::VectorXcd::Zero(node_count * stride);
   std::vector<Complex> spatial(grid.NumberOfCoLatitudes() *
                                    grid.NumberOfLongitudes(), Complex{});
   std::vector<Complex> coefficients(stride, Complex{});

   for (std::size_t element = 0; element < element_count; ++element) {
      const std::size_t layer = mesh.LayerNumber(element);
      for (int j = 0; j <= order; ++j) {
         const double radius = mesh.NodeRadius(element, j);
         const double mapped_radius =
             layer + 1 < mesh.NumberOfLayers() ? radius : mesh.PlanetRadius();
         std::size_t point = 0;
         for (const double theta : grid.CoLatitudes()) {
            for (const double phi : grid.Longitudes()) {
               const double taper =
                   layer + 1 < mesh.NumberOfLayers()
                       ? 1.0
                       : (mesh.OuterRadius() - radius) /
                             (mesh.OuterRadius() - mesh.PlanetRadius());
               const auto displacement =
                   mapping.RadialMapping(layer)(mapped_radius, theta, phi) *
                   taper;
               spatial[point] =
                   model.Density_Point(element, j, point) * displacement;
               ++point;
            }
         }
         std::fill(coefficients.begin(), coefficients.end(), Complex{});
         grid.ForwardTransformation(degree, 0, spatial, coefficients);

         // GaussDerivative is with respect to the reference coordinate. Its
         // 2/Delta-r factor cancels the element map's Delta-r/2 quadrature.
         const double radial_weight = model.q().W(j) * radius * radius;
         for (int i = 0; i <= order; ++i) {
            const double weak_weight =
                scale * model.GaussDerivative(i, j) * radial_weight;
            const std::size_t global_node = element * order + i;
            force.segment(global_node * stride, stride).noalias() +=
                weak_weight * Eigen::Map<const Eigen::VectorXcd>(
                                  coefficients.data(),
                                  static_cast<Eigen::Index>(stride));
         }
      }
   }
   return force;
}

// Pull the radial displacement back against the radial derivative of a field.
template <class mapclass>
   requires PlanetaryModel::RadialMappingClass<mapclass>
inline Eigen::VectorXcd
AdvectiveBoundaryPerturbation(GeneralEarthModels::Density3D &model,
                              mapclass &mapping,
                              const Eigen::VectorXcd &potential) {
   using Complex = std::complex<double>;
   const auto &mesh = model.Node_Information();
   const auto &grid = model.GSH_Grid();
   const std::size_t element_count = model.Num_Elements();
   const int order = model.q().N() - 1;
   const int degree = grid.MaxDegree();
   const std::size_t stride =
       GSHTrans::GSHIndices<GSHTrans::All>(degree, degree, 0).Size();
   const std::size_t node_count = element_count * order + 1;
   Eigen::VectorXcd advected = Eigen::VectorXcd::Zero(node_count * stride);
   const std::size_t spatial_size =
       grid.NumberOfCoLatitudes() * grid.NumberOfLongitudes();
   std::vector<Complex> derivative_coefficients(stride, Complex{});
   std::vector<Complex> spatial_derivative(spatial_size, Complex{});
   std::vector<Complex> spatial_displacement(spatial_size, Complex{});
   std::vector<Complex> spatial_product(spatial_size, Complex{});
   std::vector<Complex> product_coefficients(stride, Complex{});

   for (std::size_t element = 0; element < element_count; ++element) {
      const std::size_t layer = mesh.LayerNumber(element);
      const double inverse_width = 2.0 / mesh.ElementWidth(element);
      const int first_node = element == 0 ? 1 : 0;
      const int end_node = element + 1 == element_count ? order + 1 : order;
      for (int j = first_node; j < end_node; ++j) {
         std::fill(derivative_coefficients.begin(),
                   derivative_coefficients.end(), Complex{});
         for (int k = 0; k <= order; ++k) {
            const std::size_t global_node = element * order + k;
            for (std::size_t coefficient = 0; coefficient < stride;
                 ++coefficient) {
               derivative_coefficients[coefficient] +=
                   model.GaussDerivative(k, j) *
                   potential[global_node * stride + coefficient];
            }
         }
         for (auto &value : derivative_coefficients) value *= inverse_width;
         std::fill(spatial_derivative.begin(), spatial_derivative.end(), Complex{});
         grid.InverseTransformation(degree, 0, derivative_coefficients,
                                    spatial_derivative);

         const double radius = mesh.NodeRadius(element, j);
         const double mapped_radius =
             layer + 1 < mesh.NumberOfLayers() ? radius : mesh.PlanetRadius();
         std::size_t point = 0;
         for (const double theta : grid.CoLatitudes()) {
            for (const double phi : grid.Longitudes()) {
               const double taper =
                   layer + 1 < mesh.NumberOfLayers()
                       ? 1.0
                       : (mesh.OuterRadius() - radius) /
                             (mesh.OuterRadius() - mesh.PlanetRadius());
               spatial_displacement[point] =
                   mapping.RadialMapping(layer)(mapped_radius, theta, phi) *
                   taper;
               spatial_product[point] =
                   spatial_displacement[point] * spatial_derivative[point];
               ++point;
            }
         }
         std::fill(product_coefficients.begin(), product_coefficients.end(),
                   Complex{});
         grid.ForwardTransformation(degree, 0, spatial_product,
                                    product_coefficients);
         const std::size_t global_node = element * order + j;
         for (std::size_t coefficient = 0; coefficient < stride; ++coefficient) {
            advected[global_node * stride + coefficient] =
                product_coefficients[coefficient];
         }
      }
   }
   return advected;
}

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
   detail::ValidateSensitivityTolerance(relerr);
   using Complex = std::complex<double>;
   using CONJG = Eigen::ConjugateGradient<
       MatrixReplacement3D<Complex>, Eigen::Lower | Eigen::Upper,
       Eigen::SphericalGeometryPreconditioner<Complex>>;

   MatrixReplacement3D<std::complex<double>> mymatrix(inp_model);

   const Eigen::VectorXcd vec_fullforce = detail::SensitivitySelectorRhs(
       inp_model, std::vector<int>{l}, std::vector<int>{m},
       std::vector<double>{1.0});

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
   detail::ValidateSensitivityTolerance(relerr);
   using Complex = std::complex<double>;
   using BICGSTAB =
       Eigen::BiCGSTAB<MatrixReplacement3D<std::complex<double>>,
                       Eigen::SphericalGeometryPreconditioner<Complex>>;
   MatrixReplacement3D<std::complex<double>> mymatrix(inp_model);

   const Eigen::VectorXcd vec_fullforce = detail::SensitivitySelectorRhs(
       inp_model, l, m, multval);

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
