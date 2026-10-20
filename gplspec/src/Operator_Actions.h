#ifndef GPLSPEC_OPERATOR_ACTIONS_H
#define GPLSPEC_OPERATOR_ACTIONS_H

#include <Eigen/Dense>
#include <FFTWpp/Ranges>
#include <GSHTrans/All>

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "GeneralModels/Earth_Density_Models_3D.h"
#include "GeneralModels/Earth_General_Models_1D.h"
#include "Shared_Utilities.h"

namespace GPLSpec::OperatorActions {

using Complex = std::complex<double>;

inline std::size_t ScalarIndex(int degree, int order) {
   return static_cast<std::size_t>(degree * degree + degree + order);
}

inline std::size_t SpinIndex(int degree, int order) {
   return static_cast<std::size_t>(degree * degree + degree + order - 1);
}

inline double AngularFactor(int degree) {
   return std::sqrt(static_cast<double>(degree) * (degree + 1) / 2.0);
}

template <class Model, class Grid, class Derived>
Eigen::VectorXcd ApplySpherical1D(const Model &model, const Grid &grid,
                                  const Eigen::MatrixBase<Derived> &input,
                                  bool include_boundary) {
   const auto &mesh = model.Node_Information();
   const auto quadrature = model.q();
   const auto derivative = model.GaussDerivative();
   const int degree_max = grid.MaxDegree();
   const std::size_t harmonic_count =
       static_cast<std::size_t>(degree_max + 1) * (degree_max + 1);
   const std::size_t radial_count = mesh.NumberOfGlobalRadialNodes();
   const Eigen::Index dimension = radial_count * harmonic_count;
   if (input.size() != dimension)
      throw std::invalid_argument("spherical operator input has wrong size");

   const Eigen::VectorXcd x = input.template cast<Complex>();
   Eigen::VectorXcd result = Eigen::VectorXcd::Zero(dimension);
   const std::size_t spin_count =
       GSHTrans::GSHIndices<GSHTrans::All>(degree_max, degree_max, 1).Size();
   const std::size_t spatial_count =
       grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();
   Eigen::Matrix3cd tensor;
   tensor << 0.0, 0.0, -1.0,
             0.0, 1.0, 0.0,
            -1.0, 0.0, 0.0;

   std::vector<Complex> radial_gradient(harmonic_count),
       plus_gradient(spin_count), minus_gradient(spin_count),
       harmonic_minus(spin_count), harmonic_zero(harmonic_count),
       harmonic_plus(spin_count);
   FFTWpp::vector<Complex> spatial_minus(spatial_count),
       spatial_zero(spatial_count), spatial_plus(spatial_count),
       flux_minus(spatial_count), flux_zero(spatial_count),
       flux_plus(spatial_count);

   for (std::size_t element = 0; element < mesh.NumberOfElements(); ++element) {
      const double width = mesh.ElementWidth(element);
      for (std::size_t node = 0; node < mesh.NumNodesWithinElement(); ++node) {
         const double radius = detail::StandardIntervalMap(
             quadrature.X(node), mesh.ElementLowerRadius(element),
             mesh.ElementUpperRadius(element));
         std::fill(radial_gradient.begin(), radial_gradient.end(), Complex{});
         std::fill(plus_gradient.begin(), plus_gradient.end(), Complex{});
         std::fill(minus_gradient.begin(), minus_gradient.end(), Complex{});

         for (std::size_t local = 0; local < mesh.NumNodesWithinElement(); ++local) {
            const std::size_t radial = mesh.GlobalNodeIndex(element, local);
            const double scale = derivative(local, node) * 2.0 / width;
            for (int degree = 0; degree <= degree_max; ++degree) {
               for (int order = -degree; order <= degree; ++order) {
                  radial_gradient[ScalarIndex(degree, order)] +=
                      scale * x[radial * harmonic_count + ScalarIndex(degree, order)];
               }
            }
         }
         for (int degree = 1; degree <= degree_max; ++degree) {
            // The spherical 1D tensor is block diagonal: its spin flux has no
            // radial coupling. Cancel the spin derivative's 1/r against the
            // weak angular measure r before evaluating the origin node.
            const double scale = AngularFactor(degree);
            for (int order = -degree; order <= degree; ++order) {
               const auto spin = SpinIndex(degree, order);
               const auto radial = mesh.GlobalNodeIndex(element, node);
               const Complex value = scale * x[radial * harmonic_count + ScalarIndex(degree, order)];
               plus_gradient[spin] = value;
               minus_gradient[spin] = value;
            }
         }

         std::fill(spatial_minus.begin(), spatial_minus.end(), Complex{});
         std::fill(spatial_zero.begin(), spatial_zero.end(), Complex{});
         std::fill(spatial_plus.begin(), spatial_plus.end(), Complex{});
         grid.InverseTransformation(degree_max, -1, minus_gradient, spatial_minus);
         grid.InverseTransformation(degree_max, 0, radial_gradient, spatial_zero);
         grid.InverseTransformation(degree_max, 1, plus_gradient, spatial_plus);

         std::fill(flux_minus.begin(), flux_minus.end(), Complex{});
         std::fill(flux_zero.begin(), flux_zero.end(), Complex{});
         std::fill(flux_plus.begin(), flux_plus.end(), Complex{});
         for (std::size_t point = 0; point < spatial_count; ++point) {
            detail::ContractCanonicalTensorVector(
                tensor, spatial_minus[point], spatial_zero[point],
                spatial_plus[point], flux_minus[point], flux_zero[point],
                flux_plus[point]);
         }

         std::fill(harmonic_minus.begin(), harmonic_minus.end(), Complex{});
         std::fill(harmonic_zero.begin(), harmonic_zero.end(), Complex{});
         std::fill(harmonic_plus.begin(), harmonic_plus.end(), Complex{});
         grid.ForwardTransformation(degree_max, -1, flux_minus, harmonic_minus);
         grid.ForwardTransformation(degree_max, 0, flux_zero, harmonic_zero);
         grid.ForwardTransformation(degree_max, 1, flux_plus, harmonic_plus);

         const double weight = quadrature.W(node);
         for (std::size_t local = 0; local < mesh.NumNodesWithinElement(); ++local) {
            const std::size_t radial = mesh.GlobalNodeIndex(element, local);
            const double scale = derivative(local, node) *
                                 weight * radius * radius;
            for (std::size_t harmonic = 0; harmonic < harmonic_count; ++harmonic)
               result[radial * harmonic_count + harmonic] -=
                   scale * harmonic_zero[harmonic];
         }
         const std::size_t radial = mesh.GlobalNodeIndex(element, node);
         for (int degree = 1; degree <= degree_max; ++degree) {
            const double scale = width * 0.5 * weight * AngularFactor(degree);
            for (int order = -degree; order <= degree; ++order) {
               const auto spin = SpinIndex(degree, order);
               result[radial * harmonic_count + ScalarIndex(degree, order)] -=
                   scale * (harmonic_minus[spin] + harmonic_plus[spin]);
            }
         }
      }
   }

   if (include_boundary) {
      const std::size_t radial = radial_count - 1;
      for (int degree = 0; degree <= degree_max; ++degree) {
         for (int order = -degree; order <= degree; ++order) {
            const std::size_t index = radial * harmonic_count + ScalarIndex(degree, order);
            result[index] -= (degree + 1.0) * mesh.OuterRadius() * x[index];
         }
      }
   }
   return result;
}

using TensorField3D =
    std::vector<std::vector<std::vector<Eigen::Matrix3cd>>>;

template <class Model, class Derived>
Eigen::VectorXcd ApplyMapped3D(const Model &model,
                               const TensorField3D &tensor_field,
                               const Eigen::MatrixBase<Derived> &input,
                               bool include_boundary) {
   const auto &mesh = model.Node_Information();
   const auto &grid = model.GSH_Grid();
   const auto &quadrature = model.q();
   const auto derivative = model.GaussDerivative();
   const int degree_max = grid.MaxDegree();
   const std::size_t harmonic_count =
       static_cast<std::size_t>(degree_max + 1) * (degree_max + 1);
   const std::size_t radial_count = mesh.NumberOfGlobalRadialNodes();
   const Eigen::Index dimension = radial_count * harmonic_count;
   if (input.size() != dimension)
      throw std::invalid_argument("mapped operator input has wrong size");

   const Eigen::VectorXcd x = input.template cast<Complex>();
   Eigen::VectorXcd result = Eigen::VectorXcd::Zero(dimension);
   const std::size_t spin_count =
       GSHTrans::GSHIndices<GSHTrans::All>(degree_max, degree_max, 1).Size();
   const std::size_t spatial_count =
       grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();

   std::vector<Complex> radial_gradient(harmonic_count),
       plus_gradient(spin_count), minus_gradient(spin_count),
       harmonic_minus(spin_count), harmonic_zero(harmonic_count),
       harmonic_plus(spin_count);
   FFTWpp::vector<Complex> spatial_minus(spatial_count),
       spatial_zero(spatial_count), spatial_plus(spatial_count),
       flux_minus(spatial_count), flux_zero(spatial_count),
       flux_plus(spatial_count);

   for (std::size_t element = 0; element < mesh.NumberOfElements(); ++element) {
      const double width = mesh.ElementWidth(element);
      for (std::size_t node = 0; node < mesh.NumNodesWithinElement(); ++node) {
         const double radius = mesh.NodeRadius(element, node);
         std::fill(radial_gradient.begin(), radial_gradient.end(), Complex{});
         std::fill(plus_gradient.begin(), plus_gradient.end(), Complex{});
         std::fill(minus_gradient.begin(), minus_gradient.end(), Complex{});

         for (std::size_t local = 0; local < mesh.NumNodesWithinElement(); ++local) {
            const std::size_t radial = mesh.GlobalNodeIndex(element, local);
            const double scale = derivative(local, node) *
                                 radius * 2.0 / width;
            for (int degree = 0; degree <= degree_max; ++degree) {
               for (int order = -degree; order <= degree; ++order) {
                  radial_gradient[ScalarIndex(degree, order)] +=
                      scale * x[radial * harmonic_count + ScalarIndex(degree, order)];
               }
            }
         }
         const std::size_t radial = mesh.GlobalNodeIndex(element, node);
         for (int degree = 1; degree <= degree_max; ++degree) {
            const double scale = AngularFactor(degree);
            for (int order = -degree; order <= degree; ++order) {
               const auto spin = SpinIndex(degree, order);
               const Complex value = scale * x[radial * harmonic_count + ScalarIndex(degree, order)];
               plus_gradient[spin] = value;
               minus_gradient[spin] = value;
            }
         }

         std::fill(spatial_minus.begin(), spatial_minus.end(), Complex{});
         std::fill(spatial_zero.begin(), spatial_zero.end(), Complex{});
         std::fill(spatial_plus.begin(), spatial_plus.end(), Complex{});
         grid.InverseTransformation(degree_max, -1, minus_gradient, spatial_minus);
         grid.InverseTransformation(degree_max, 0, radial_gradient, spatial_zero);
         grid.InverseTransformation(degree_max, 1, plus_gradient, spatial_plus);

         std::fill(flux_minus.begin(), flux_minus.end(), Complex{});
         std::fill(flux_zero.begin(), flux_zero.end(), Complex{});
         std::fill(flux_plus.begin(), flux_plus.end(), Complex{});
         for (std::size_t point = 0; point < spatial_count; ++point) {
            detail::ContractCanonicalTensorVector(
                tensor_field[element][node][point], spatial_minus[point],
                spatial_zero[point], spatial_plus[point], flux_minus[point],
                flux_zero[point], flux_plus[point]);
         }

         std::fill(harmonic_minus.begin(), harmonic_minus.end(), Complex{});
         std::fill(harmonic_zero.begin(), harmonic_zero.end(), Complex{});
         std::fill(harmonic_plus.begin(), harmonic_plus.end(), Complex{});
         grid.ForwardTransformation(degree_max, -1, flux_minus, harmonic_minus);
         grid.ForwardTransformation(degree_max, 0, flux_zero, harmonic_zero);
         grid.ForwardTransformation(degree_max, 1, flux_plus, harmonic_plus);

         const double weight = quadrature.W(node);
         for (std::size_t local = 0; local < mesh.NumNodesWithinElement(); ++local) {
            const std::size_t output_radial = mesh.GlobalNodeIndex(element, local);
            const double scale = derivative(local, node) *
                                 weight * radius;
            for (std::size_t harmonic = 0; harmonic < harmonic_count; ++harmonic)
               result[output_radial * harmonic_count + harmonic] -=
                   scale * harmonic_zero[harmonic];
         }
         const double angular_scale = width * 0.5 * weight;
         for (int degree = 1; degree <= degree_max; ++degree) {
            const double scale = angular_scale * AngularFactor(degree);
            for (int order = -degree; order <= degree; ++order) {
               const auto spin = SpinIndex(degree, order);
               result[radial * harmonic_count + ScalarIndex(degree, order)] -=
                   scale * (harmonic_minus[spin] + harmonic_plus[spin]);
            }
         }
      }
   }

   if (include_boundary) {
      const std::size_t radial = radial_count - 1;
      for (int degree = 0; degree <= degree_max; ++degree) {
         for (int order = -degree; order <= degree; ++order) {
            const std::size_t index = radial * harmonic_count + ScalarIndex(degree, order);
            result[index] -= (degree + 1.0) * mesh.OuterRadius() * x[index];
         }
      }
   }
   return result;
}

}  // namespace GPLSpec::OperatorActions

#endif
