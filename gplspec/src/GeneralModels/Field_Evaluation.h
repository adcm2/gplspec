#ifndef GPLSPEC_FIELD_EVALUATION_H
#define GPLSPEC_FIELD_EVALUATION_H

#include "Density_Model_Return.h"
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace GPLSpec::FieldEvaluation {

enum class Quantity {
   Potential,
   SensitivityKernel,
   ReferentialDensity,
   PhysicalDensity,
   DegreePower,
   SurfaceRadius
};
enum class Representation {
   AngularGrid,
   HarmonicCoefficients,
   RotatedSlice,
   DegreePower,
   RadialProfile
};
enum class TraceSide { None, Inner, Outer };
enum class BoundaryKind { None, Numerical, Material, BodyVacuum, SurfaceWithoutExterior };

inline constexpr std::size_t NoIndex = std::numeric_limits<std::size_t>::max();

struct Rotation {
   double alpha = 0.0;
   double beta = 0.0;
   double gamma = 0.0;
};

struct Sample {
   Quantity quantity = Quantity::Potential;
   Representation representation = Representation::AngularGrid;
   BoundaryKind boundary = BoundaryKind::None;
   TraceSide side = TraceSide::None;
   std::size_t element = NoIndex;
   std::size_t numerical_interface = NoIndex;
   std::size_t material_interface = NoIndex;
   std::size_t inner_element = NoIndex;
   std::size_t outer_element = NoIndex;
   int material = -1;
   int inner_material = -1;
   int outer_material = -1;
   int degree = -1;
   int order = 0;
   int radial_node = -1;
   double reference_radius_m = std::numeric_limits<double>::quiet_NaN();
   double physical_radius_m = std::numeric_limits<double>::quiet_NaN();
   double theta_rad = std::numeric_limits<double>::quiet_NaN();
   double phi_rad = std::numeric_limits<double>::quiet_NaN();
   double x_m = std::numeric_limits<double>::quiet_NaN();
   double y_m = std::numeric_limits<double>::quiet_NaN();
   double z_m = std::numeric_limits<double>::quiet_NaN();
   double reference_x_m = std::numeric_limits<double>::quiet_NaN();
   double reference_y_m = std::numeric_limits<double>::quiet_NaN();
   double reference_z_m = std::numeric_limits<double>::quiet_NaN();
   std::complex<double> value = 0.0;
};

namespace detail {

inline void ValidateQuantity(Quantity quantity) {
   switch (quantity) {
   case Quantity::Potential:
   case Quantity::SensitivityKernel:
   case Quantity::ReferentialDensity:
   case Quantity::PhysicalDensity:
   case Quantity::DegreePower:
   case Quantity::SurfaceRadius:
      return;
   }
   throw std::invalid_argument("unknown field quantity");
}

inline bool RequiresSolution(Quantity quantity) {
   return quantity == Quantity::Potential ||
          quantity == Quantity::SensitivityKernel ||
          quantity == Quantity::DegreePower;
}

inline bool IsDensity(Quantity quantity) {
   return quantity == Quantity::ReferentialDensity ||
          quantity == Quantity::PhysicalDensity;
}

inline bool IsPotential(Quantity quantity) {
   return quantity == Quantity::Potential ||
          quantity == Quantity::SensitivityKernel;
}

inline bool IsSurfaceTrace(const Sample &sample) {
   return sample.boundary == BoundaryKind::BodyVacuum ||
          sample.boundary == BoundaryKind::SurfaceWithoutExterior;
}

inline void ValidateRotation(const Rotation &rotation) {
   if (!std::isfinite(rotation.alpha) || !std::isfinite(rotation.beta) ||
       !std::isfinite(rotation.gamma)) {
      throw std::invalid_argument("rotation angles must be finite");
   }
}

struct InterfaceInfo {
   BoundaryKind kind = BoundaryKind::Numerical;
   std::size_t id = NoIndex;
   std::size_t material_id = NoIndex;
   std::size_t inner_element = NoIndex;
   std::size_t outer_element = NoIndex;
   int inner_material = -1;
   int outer_material = -1;
};

inline Sample SurfaceOnlyTrace(const Radial_Tools::RadialMesh &mesh,
                               std::size_t element, double length_norm) {
   Sample sample;
   sample.element = element;
   sample.material = mesh.LayerNumber(static_cast<int>(element));
   sample.radial_node = mesh.PolynomialOrder();
   sample.reference_radius_m =
       mesh.NodeRadius(static_cast<int>(element), mesh.PolynomialOrder()) *
       length_norm;
   sample.boundary = BoundaryKind::SurfaceWithoutExterior;
   sample.side = TraceSide::Inner;
   sample.inner_element = element;
   sample.inner_material = sample.material;
   return sample;
}

inline std::vector<InterfaceInfo>
Interfaces(const GeneralEarthModels::Density3D &model) {
   const auto &mesh = model.Node_Information();
   std::vector<InterfaceInfo> result(mesh.NumberOfElements() > 0
                                         ? mesh.NumberOfElements() - 1
                                         : 0);
   std::size_t next_material_id = 0;
   for (std::size_t i = 0; i < result.size(); ++i) {
      const int inner = mesh.LayerNumber(static_cast<int>(i));
      const int outer = mesh.LayerNumber(static_cast<int>(i + 1));
      auto &info = result[i];
      info.id = i;
      info.inner_element = i;
      info.outer_element = i + 1;
      info.inner_material = inner;
      info.outer_material = outer;
      const double radius = mesh.NodeRadius(static_cast<int>(i),
                                            mesh.PolynomialOrder());
      const bool at_body_surface =
          std::abs(radius - mesh.PlanetRadius()) <=
          16.0 * std::numeric_limits<double>::epsilon() *
              std::max(1.0, std::abs(mesh.PlanetRadius()));
      const bool material_change = inner != outer;
      if (material_change) {
         info.material_id = next_material_id++;
      }
      if (at_body_surface && material_change &&
          outer == mesh.NumberOfLayers() - 1) {
         info.kind = BoundaryKind::BodyVacuum;
      } else if (material_change) {
         info.kind = BoundaryKind::Material;
      }
   }
   return result;
}

inline Sample AtLowerTrace(const GeneralEarthModels::Density3D &model,
                           std::size_t element,
                           const std::vector<InterfaceInfo> &interfaces) {
   const auto &mesh = model.Node_Information();
   Sample sample;
   sample.element = element;
   sample.material = mesh.LayerNumber(static_cast<int>(element));
   sample.radial_node = 0;
   sample.reference_radius_m = mesh.NodeRadius(static_cast<int>(element), 0) *
                               model.LengthNorm();
   if (element > 0) {
      const auto &info = interfaces[element - 1];
      sample.boundary = info.kind;
      sample.side = TraceSide::Outer;
      sample.numerical_interface = info.id;
      sample.material_interface = info.material_id;
      sample.inner_element = info.inner_element;
      sample.outer_element = info.outer_element;
      sample.inner_material = info.inner_material;
      sample.outer_material = info.outer_material;
   }
   return sample;
}

inline Sample AtUpperTrace(const GeneralEarthModels::Density3D &model,
                           std::size_t element,
                           const std::vector<InterfaceInfo> &interfaces) {
   const auto &mesh = model.Node_Information();
   Sample sample;
   sample.element = element;
   sample.material = mesh.LayerNumber(static_cast<int>(element));
   sample.radial_node = mesh.PolynomialOrder();
   sample.reference_radius_m =
       mesh.NodeRadius(static_cast<int>(element), mesh.PolynomialOrder()) *
       model.LengthNorm();
   if (element + 1 < static_cast<std::size_t>(mesh.NumberOfElements())) {
      const auto &info = interfaces[element];
      sample.boundary = info.kind;
      sample.side = TraceSide::Inner;
      sample.numerical_interface = info.id;
      sample.material_interface = info.material_id;
      sample.inner_element = info.inner_element;
      sample.outer_element = info.outer_element;
      sample.inner_material = info.inner_material;
      sample.outer_material = info.outer_material;
   } else {
      const double radius = mesh.NodeRadius(static_cast<int>(element),
                                            mesh.PolynomialOrder());
      const bool at_surface =
          std::abs(radius - mesh.PlanetRadius()) <=
          16.0 * std::numeric_limits<double>::epsilon() *
              std::max(1.0, std::abs(mesh.PlanetRadius()));
      if (at_surface) {
         sample = SurfaceOnlyTrace(mesh, element, model.LengthNorm());
      }
   }
   return sample;
}

inline void SetCoordinates(Sample &sample, double physical_radius_m,
                           double theta, double phi) {
   sample.physical_radius_m = physical_radius_m;
   sample.theta_rad = theta;
   sample.phi_rad = phi;
   sample.x_m = physical_radius_m * std::sin(theta) * std::cos(phi);
   sample.y_m = physical_radius_m * std::sin(theta) * std::sin(phi);
   sample.z_m = physical_radius_m * std::cos(theta);
   sample.reference_x_m = sample.reference_radius_m * std::sin(theta) * std::cos(phi);
   sample.reference_y_m = sample.reference_radius_m * std::sin(theta) * std::sin(phi);
   sample.reference_z_m = sample.reference_radius_m * std::cos(theta);
}

inline void RotateCoefficients(
    int lmax,
    const std::vector<std::complex<double>> &input,
    std::vector<std::complex<double>> &output,
    const std::vector<Eigen::MatrixXcd> &matrices) {
   const std::size_t size = static_cast<std::size_t>((lmax + 1) * (lmax + 1));
   if (input.size() != size) {
      throw std::invalid_argument("harmonic vector has the wrong size");
   }
   output.resize(size);
   for (int l = 0; l <= lmax; ++l) {
      const std::size_t offset = static_cast<std::size_t>(l * l);
      for (int m = -l; m <= l; ++m) {
         std::complex<double> value = 0.0;
         for (int mp = -l; mp <= l; ++mp) {
            value += matrices[l](m + l, mp + l) * input[offset + mp + l];
         }
         output[offset + m + l] = value;
      }
   }
}

inline std::vector<std::complex<double>>
HarmonicsAt(double theta, int lmax) {
   const GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All,
                         GSHTrans::Single, GSHTrans::Single,
                         GSHTrans::ColumnMajor>
       wigner(lmax, lmax, 0, theta);
   std::vector<std::complex<double>> values(
       static_cast<std::size_t>((lmax + 1) * (lmax + 1)));
   const auto d = wigner[0, 0];
   for (int l = 0; l <= lmax; ++l) {
      for (int m = -l; m <= l; ++m) {
         values[static_cast<std::size_t>(l * l + l + m)] =
             d[l][m];
      }
   }
   return values;
}

inline std::complex<double>
EvaluateHarmonicSeries(int lmax,
                       const std::vector<std::complex<double>> &coefficients,
                       const std::vector<std::complex<double>> &at_theta,
                       double phi) {
   std::complex<double> value = 0.0;
   for (int l = 0; l <= lmax; ++l) {
      for (int m = -l; m <= l; ++m) {
         value += coefficients[static_cast<std::size_t>(l * l + l + m)] *
                  at_theta[static_cast<std::size_t>(l * l + l + m)] *
                  std::polar(1.0, static_cast<double>(m) * phi);
      }
   }
   return value;
}

inline std::size_t GlobalOffset(const GeneralEarthModels::Density3D &model,
                                std::size_t element, int node,
                                std::size_t stride) {
   return static_cast<std::size_t>(model.Node_Information().GlobalNodeIndex(
              static_cast<int>(element), node)) *
          stride;
}

}   // namespace detail

// Emits both endpoint traces of every element. Both sides of a shared
// interface are emitted; potential values come from the
// same global degree of freedom. The callback can serialize or reduce samples
// without materializing the full field.
template <class Callback>
void ForEachAngularGridSample(
    const GeneralEarthModels::Density3D &model,
    std::span<const std::complex<double>> flat_solution, Quantity quantity,
    Callback &&callback, const Rotation *rotation = nullptr) {
   detail::ValidateQuantity(quantity);
   if (rotation != nullptr) detail::ValidateRotation(*rotation);
   if (quantity == Quantity::DegreePower)
      throw std::invalid_argument("degree power has no angular spatial samples");
   const bool needs_solution = detail::RequiresSolution(quantity);
   const int lmax = model.GSH_Grid().MaxDegree();
   const std::size_t stride = static_cast<std::size_t>((lmax + 1) * (lmax + 1));
   const auto &mesh = model.Node_Information();
   if (needs_solution && flat_solution.size() !=
                             static_cast<std::size_t>(mesh.NumberOfGlobalRadialNodes()) *
                                 stride) {
      throw std::invalid_argument("flat solution does not match model dimensions");
   }
   if (!needs_solution && !flat_solution.empty()) {
      throw std::invalid_argument("density evaluation does not take a solution");
   }
   const auto interfaces = detail::Interfaces(model);
   const auto spatial_size = model.GSH_Grid().Longitudes().size() *
                             model.GSH_Grid().CoLatitudes().size();
   const auto longitudes = model.GSH_Grid().Longitudes();
   const auto colatitudes = model.GSH_Grid().CoLatitudes();
   std::vector<std::complex<double>> coefficients(stride), rotated(stride);
   std::vector<std::complex<double>> nonnegative(
       static_cast<std::size_t>((lmax + 1) * (lmax + 2) / 2));
   std::vector<double> spatial(spatial_size), mapping(spatial_size);
   std::vector<std::complex<double>> complex_spatial(spatial_size);
   std::vector<Eigen::MatrixXcd> matrices;
   if (rotation != nullptr) {
      matrices = GPLSpec::detail::RotationOutputMatrices(
          lmax, rotation->alpha, rotation->beta, rotation->gamma);
   }

   auto emit = [&](std::size_t element, int node, Sample base) {
      if (quantity == Quantity::SurfaceRadius &&
          !detail::IsSurfaceTrace(base))
         return;
      const std::size_t offset = detail::GlobalOffset(model, element, node, stride);
      const auto mapping_node =
          model.MappingAtRadialNode(static_cast<int>(element), node);
      std::copy(mapping_node.begin(), mapping_node.end(), mapping.begin());
      if (rotation != nullptr) {
         std::fill(nonnegative.begin(), nonnegative.end(),
                   std::complex<double>{});
         model.GSH_Grid().ForwardTransformation(lmax, 0, mapping, nonnegative);
         GPLSpec::detail::ExpandRealScalarCoefficients(lmax, nonnegative,
                                                       coefficients);
         detail::RotateCoefficients(lmax, coefficients, rotated, matrices);
         model.GSH_Grid().InverseTransformation(lmax, 0, rotated,
                                                complex_spatial);
         for (std::size_t k = 0; k < spatial_size; ++k)
            mapping[k] = complex_spatial[k].real();
      }
      if (detail::IsDensity(quantity)) {
         for (std::size_t k = 0; k < spatial_size; ++k) {
            const double rho = model.Density_Point(static_cast<int>(element), node,
                                                    static_cast<int>(k));
            spatial[k] = quantity == Quantity::PhysicalDensity
                             ? rho / model.Jacobian_Point(static_cast<int>(element),
                                                         node, static_cast<int>(k))
                             : rho;
         }
         if (rotation != nullptr) {
            std::fill(nonnegative.begin(), nonnegative.end(),
                      std::complex<double>{});
            model.GSH_Grid().ForwardTransformation(lmax, 0, spatial, nonnegative);
            GPLSpec::detail::ExpandRealScalarCoefficients(lmax, nonnegative,
                                                          coefficients);
            detail::RotateCoefficients(lmax, coefficients, rotated, matrices);
            model.GSH_Grid().InverseTransformation(lmax, 0, rotated,
                                                   complex_spatial);
            for (std::size_t k = 0; k < spatial_size; ++k)
               spatial[k] = complex_spatial[k].real();
         }
      } else if (detail::IsPotential(quantity)) {
         std::copy_n(flat_solution.begin() + static_cast<std::ptrdiff_t>(offset),
                     stride, coefficients.begin());
         if (rotation != nullptr) {
            detail::RotateCoefficients(lmax, coefficients, rotated, matrices);
            model.GSH_Grid().InverseTransformation(lmax, 0, rotated,
                                                   complex_spatial);
         } else {
            model.GSH_Grid().InverseTransformation(lmax, 0, coefficients,
                                                   complex_spatial);
         }
      }
      for (std::size_t i = 0; i < colatitudes.size(); ++i) {
         for (std::size_t j = 0; j < longitudes.size(); ++j) {
            const std::size_t k = i * longitudes.size() + j;
            Sample sample = base;
            sample.quantity = quantity;
            sample.representation = Representation::AngularGrid;
            const double theta = colatitudes[i];
            const double phi = longitudes[j];
            const double physical_radius =
                (mesh.NodeRadius(static_cast<int>(element), node) + mapping[k]) *
                model.LengthNorm();
            detail::SetCoordinates(sample, physical_radius, theta, phi);
            if (quantity == Quantity::Potential) {
               sample.value = model.PotentialNorm() * complex_spatial[k];
            } else if (quantity == Quantity::SensitivityKernel) {
               sample.value = model.PotentialNorm() / model.MassNorm() * complex_spatial[k];
            } else if (detail::IsDensity(quantity)) {
               sample.value = model.DensityNorm() * spatial[k];
            } else {
               sample.value = physical_radius;
            }
            callback(sample);
         }
      }
   };

   for (std::size_t e = 0; e < static_cast<std::size_t>(mesh.NumberOfElements()); ++e) {
      emit(e, 0, detail::AtLowerTrace(model, e, interfaces));
      emit(e, mesh.PolynomialOrder(), detail::AtUpperTrace(model, e, interfaces));
   }
}

template <class Callback>
void ForEachAngularGridSample(
    const GeneralEarthModels::Density3D &model,
    const Eigen::VectorXcd &flat_solution, Quantity quantity,
    Callback &&callback, const Rotation *rotation = nullptr) {
   ForEachAngularGridSample(
       model, std::span<const std::complex<double>>(flat_solution.data(),
                                                   static_cast<std::size_t>(flat_solution.size())),
       quantity, std::forward<Callback>(callback), rotation);
}

template <class Callback>
void ForEachRadialHarmonicSample(
    const GeneralEarthModels::Density3D &model,
    std::span<const std::complex<double>> flat_solution, Quantity quantity,
    Callback &&callback) {
   detail::ValidateQuantity(quantity);
   if (quantity == Quantity::SurfaceRadius)
      throw std::invalid_argument("surface radius has no radial harmonic samples");
   const bool needs_solution = detail::RequiresSolution(quantity);
   const int lmax = model.GSH_Grid().MaxDegree();
   const std::size_t stride = static_cast<std::size_t>((lmax + 1) * (lmax + 1));
   const auto &mesh = model.Node_Information();
   if (needs_solution && flat_solution.size() !=
                             static_cast<std::size_t>(mesh.NumberOfGlobalRadialNodes()) *
                                 stride) {
      throw std::invalid_argument("flat solution does not match model dimensions");
   }
   if (!needs_solution && !flat_solution.empty()) {
      throw std::invalid_argument("density evaluation does not take a solution");
   }
   const auto interfaces = detail::Interfaces(model);
   std::vector<std::complex<double>> coefficients(stride);
   std::vector<std::complex<double>> nonnegative(
       static_cast<std::size_t>((lmax + 1) * (lmax + 2) / 2));
   std::vector<double> spatial(model.GSH_Grid().Longitudes().size() *
                               model.GSH_Grid().CoLatitudes().size());

   auto emit = [&](std::size_t element, int node, Sample base) {
      const std::size_t offset = detail::GlobalOffset(model, element, node, stride);
      if (needs_solution) {
         std::copy_n(flat_solution.begin() + static_cast<std::ptrdiff_t>(offset),
                     stride, coefficients.begin());
      } else if (detail::IsDensity(quantity)) {
         for (std::size_t k = 0; k < spatial.size(); ++k) {
            double value = model.Density_Point(static_cast<int>(element), node,
                                               static_cast<int>(k));
            if (quantity == Quantity::PhysicalDensity) {
               value /= model.Jacobian_Point(static_cast<int>(element), node,
                                             static_cast<int>(k));
            }
            spatial[k] = value;
         }
         std::fill(nonnegative.begin(), nonnegative.end(),
                   std::complex<double>{});
         model.GSH_Grid().ForwardTransformation(lmax, 0, spatial, nonnegative);
         GPLSpec::detail::ExpandRealScalarCoefficients(lmax, nonnegative,
                                                       coefficients);
      }
      const double scale = quantity == Quantity::Potential
                               ? model.PotentialNorm()
                               : quantity == Quantity::SensitivityKernel
                                     ? model.PotentialNorm() / model.MassNorm()
                                     : model.DensityNorm();
      base.quantity = quantity;
      base.representation = quantity == Quantity::DegreePower
                                ? Representation::DegreePower
                                : Representation::HarmonicCoefficients;
      for (int l = 0; l <= lmax; ++l) {
         double degree_power = 0.0;
         if (quantity == Quantity::DegreePower) {
            const double monopole_power = std::norm(coefficients[0]);
            if (!std::isfinite(monopole_power) || monopole_power == 0.0)
               throw std::domain_error(
                   "degree power is undefined for a zero or non-finite monopole");
            for (int m = -l; m <= l; ++m)
               degree_power +=
                   std::norm(coefficients[static_cast<std::size_t>(l * l + l + m)]) /
                   monopole_power;
         }
         if (quantity == Quantity::DegreePower) {
            Sample sample = base;
            sample.degree = l;
            sample.order = 0;
            sample.value = degree_power;
            callback(sample);
            continue;
         }
         for (int m = -l; m <= l; ++m) {
            Sample sample = base;
            sample.degree = l;
            sample.order = m;
            sample.value = scale * coefficients[static_cast<std::size_t>(l * l + l + m)];
            callback(sample);
         }
      }
   };

   for (std::size_t e = 0; e < static_cast<std::size_t>(mesh.NumberOfElements()); ++e) {
      emit(e, 0, detail::AtLowerTrace(model, e, interfaces));
      emit(e, mesh.PolynomialOrder(), detail::AtUpperTrace(model, e, interfaces));
   }
}

template <class Callback>
void ForEachElementBoundaryHarmonics(
    const GeneralEarthModels::Density3D &model,
    const std::vector<std::vector<std::complex<double>>> &boundary_values,
    Quantity quantity, Callback &&callback) {
   detail::ValidateQuantity(quantity);
   if (quantity != Quantity::Potential && quantity != Quantity::SensitivityKernel)
      throw std::invalid_argument("boundary harmonics require potential or sensitivity");
   const auto &mesh = model.Node_Information();
   const int lmax = model.GSH_Grid().MaxDegree();
   const std::size_t stride = static_cast<std::size_t>((lmax + 1) * (lmax + 1));
   const auto element_count = static_cast<std::size_t>(mesh.NumberOfElements());
   if (boundary_values.size() != element_count + 1)
      throw std::invalid_argument("boundary harmonics require one value per mesh boundary");
   for (const auto &row : boundary_values)
      if (row.size() != stride)
         throw std::invalid_argument("boundary harmonics have the wrong degree count");
   const auto interfaces = detail::Interfaces(model);
   const double scale = quantity == Quantity::Potential
                            ? model.PotentialNorm()
                            : model.PotentialNorm() / model.MassNorm();
   auto emit = [&](std::size_t boundary, Sample base) {
      base.quantity = quantity;
      base.representation = Representation::HarmonicCoefficients;
      for (int l = 0; l <= lmax; ++l)
         for (int m = -l; m <= l; ++m) {
            Sample sample = base;
            sample.degree = l;
            sample.order = m;
            sample.value =
                scale * boundary_values[boundary][static_cast<std::size_t>(l*l+l+m)];
            callback(sample);
         }
   };
   for (std::size_t element = 0; element < element_count; ++element) {
      emit(element,
           detail::AtLowerTrace(model, element, interfaces));
      emit(element + 1,
           detail::AtUpperTrace(model, element, interfaces));
   }
}

template <class Callback>
void ForEachScalarRadialProfile(
    const GeneralEarthModels::Density3D &model,
    const std::vector<double> &global_values, Quantity quantity,
    Callback &&callback) {
   if (quantity != Quantity::Potential ||
       global_values.size() != static_cast<std::size_t>(model.Node_Information().NumberOfGlobalRadialNodes()))
      throw std::invalid_argument("scalar profile requires global potential node values");
   const auto &mesh = model.Node_Information();
   const auto interfaces = detail::Interfaces(model);
   auto emit = [&](std::size_t element, int node, Sample sample) {
      const auto global = static_cast<std::size_t>(mesh.GlobalNodeIndex(
          static_cast<int>(element), node));
      sample.quantity = quantity;
      sample.representation = Representation::RadialProfile;
      sample.degree = -1;
      sample.order = 0;
      sample.value = model.PotentialNorm() * global_values[global];
      callback(sample);
   };
   for (std::size_t element = 0;
        element < static_cast<std::size_t>(mesh.NumberOfElements()); ++element) {
      emit(element, 0, detail::AtLowerTrace(model, element, interfaces));
      emit(element, mesh.PolynomialOrder(),
           detail::AtUpperTrace(model, element, interfaces));
   }
}

template <class Callback>
void ForEachElementBoundaryMonopole(
    const GeneralEarthModels::Density3D &model,
    const std::vector<double> &boundary_values, Quantity quantity,
    Callback &&callback) {
   const auto &mesh = model.Node_Information();
   const auto element_count = static_cast<std::size_t>(mesh.NumberOfElements());
   if (quantity != Quantity::Potential || boundary_values.size() != element_count + 1)
      throw std::invalid_argument("boundary monopole requires one potential coefficient per boundary");
   const auto interfaces = detail::Interfaces(model);
   auto emit = [&](std::size_t boundary, Sample sample) {
      sample.quantity = quantity;
      sample.representation = Representation::HarmonicCoefficients;
      sample.degree = 0;
      sample.order = 0;
      sample.value = model.PotentialNorm() * boundary_values[boundary];
      callback(sample);
   };
   for (std::size_t element = 0; element < element_count; ++element) {
      emit(element, detail::AtLowerTrace(model, element, interfaces));
      emit(element + 1,
           detail::AtUpperTrace(model, element, interfaces));
   }
}

// Exports an analytic scalar profile sampled at one modeled angular ray. The
// supplied radii are normalized physical radii, as returned by
// Density3D::PhysicalRadius_Line; shared interfaces are emitted from both
// adjacent element traces with their common analytic value.
template <class Callback>
void ForEachMappedRadialProfile(
    const GeneralEarthModels::Density3D &model,
    const std::vector<double> &physical_radii_normalized,
    const std::vector<double> &normalized_potential,
    std::size_t spatial_index, Callback &&callback) {
   const auto &mesh = model.Node_Information();
   const auto element_count = static_cast<std::size_t>(mesh.NumberOfElements());
   if (physical_radii_normalized.size() != element_count + 1 ||
       normalized_potential.size() != element_count + 1)
      throw std::invalid_argument("mapped radial profile must contain all boundary nodes");
   const auto longitude_count = model.GSH_Grid().Longitudes().size();
   const auto colatitudes = model.GSH_Grid().CoLatitudes();
   const auto longitudes = model.GSH_Grid().Longitudes();
   if (spatial_index >= colatitudes.size() * longitude_count)
      throw std::invalid_argument("mapped radial profile angular index is out of range");
   const double theta = colatitudes[spatial_index / longitude_count];
   const double phi = longitudes[spatial_index % longitude_count];
   const auto interfaces = detail::Interfaces(model);
   auto emit = [&](std::size_t element, int node, std::size_t global,
                   Sample base) {
      base.quantity = Quantity::Potential;
      base.representation = Representation::RadialProfile;
      base.radial_node = node;
      base.theta_rad = theta;
      base.phi_rad = phi;
      detail::SetCoordinates(base,
          physical_radii_normalized[global] * model.LengthNorm(), theta, phi);
      base.reference_radius_m = mesh.NodeRadius(static_cast<int>(element), node) *
                                model.LengthNorm();
      base.value = model.PotentialNorm() * normalized_potential[global];
      callback(base);
   };
   for (std::size_t element = 0; element < element_count; ++element) {
      emit(element, 0, element, detail::AtLowerTrace(model, element, interfaces));
      emit(element, mesh.PolynomialOrder(), element + 1,
           detail::AtUpperTrace(model, element, interfaces));
   }
}

template <class Callback>
void ForEachRadialHarmonicSample(
    const GeneralEarthModels::Density3D &model,
    const Eigen::VectorXcd &flat_solution, Quantity quantity,
    Callback &&callback) {
   ForEachRadialHarmonicSample(
       model, std::span<const std::complex<double>>(flat_solution.data(),
                                                   static_cast<std::size_t>(flat_solution.size())),
       quantity, std::forward<Callback>(callback));
}

// Arbitrary-angle route for bounded slices such as the equator. Unlike the
// Gauss-grid path, this evaluates O(lmax^2) harmonics per requested azimuth;
// callers should keep the requested angular sample set to a slice/profile.
template <class Callback>
void ForEachRotatedSliceSample(
    const GeneralEarthModels::Density3D &model,
    std::span<const std::complex<double>> flat_solution, Quantity quantity,
    const Rotation &rotation, std::span<const double> azimuths,
    Callback &&callback) {
   detail::ValidateQuantity(quantity);
   detail::ValidateRotation(rotation);
   if (quantity == Quantity::DegreePower)
      throw std::invalid_argument("degree power has no rotated spatial samples");
   const bool needs_solution = detail::RequiresSolution(quantity);
   const int lmax = model.GSH_Grid().MaxDegree();
   const std::size_t stride = static_cast<std::size_t>((lmax + 1) * (lmax + 1));
   const auto &mesh = model.Node_Information();
   if (needs_solution && flat_solution.size() !=
                             static_cast<std::size_t>(mesh.NumberOfGlobalRadialNodes()) *
                                 stride) {
      throw std::invalid_argument("flat solution does not match model dimensions");
   }
   if (!needs_solution && !flat_solution.empty()) {
      throw std::invalid_argument("density evaluation does not take a solution");
   }
   for (double phi : azimuths) {
      if (!std::isfinite(phi))
         throw std::invalid_argument("slice azimuths must be finite");
   }
   const auto interfaces = detail::Interfaces(model);
   const auto spatial_size = model.GSH_Grid().Longitudes().size() *
                             model.GSH_Grid().CoLatitudes().size();
   std::vector<std::complex<double>> coefficients(stride), rotated(stride),
       rotated_mapping(stride);
   std::vector<std::complex<double>> nonnegative(
       static_cast<std::size_t>((lmax + 1) * (lmax + 2) / 2));
   std::vector<double> spatial(spatial_size), mapping(spatial_size);
   const auto matrices = GPLSpec::detail::RotationOutputMatrices(
       lmax, rotation.alpha, rotation.beta, rotation.gamma);
   const auto harmonics = detail::HarmonicsAt(std::numbers::pi / 2.0, lmax);

   auto emit = [&](std::size_t element, int node, Sample base) {
      if (quantity == Quantity::SurfaceRadius &&
          !detail::IsSurfaceTrace(base))
         return;
      const std::size_t offset = detail::GlobalOffset(model, element, node, stride);
      const auto mapping_node =
          model.MappingAtRadialNode(static_cast<int>(element), node);
      std::copy(mapping_node.begin(), mapping_node.end(), mapping.begin());
      std::fill(nonnegative.begin(), nonnegative.end(),
                std::complex<double>{});
      model.GSH_Grid().ForwardTransformation(lmax, 0, mapping, nonnegative);
      GPLSpec::detail::ExpandRealScalarCoefficients(lmax, nonnegative,
                                                    coefficients);
      detail::RotateCoefficients(lmax, coefficients, rotated, matrices);
      std::copy(rotated.begin(), rotated.end(), rotated_mapping.begin());
      if (needs_solution) {
         std::copy_n(flat_solution.begin() + static_cast<std::ptrdiff_t>(offset),
                     stride, coefficients.begin());
      } else if (detail::IsDensity(quantity)) {
         for (std::size_t k = 0; k < spatial_size; ++k) {
            double value = model.Density_Point(static_cast<int>(element), node,
                                               static_cast<int>(k));
            if (quantity == Quantity::PhysicalDensity)
               value /= model.Jacobian_Point(static_cast<int>(element), node,
                                             static_cast<int>(k));
            spatial[k] = value;
         }
         std::fill(nonnegative.begin(), nonnegative.end(),
                   std::complex<double>{});
         model.GSH_Grid().ForwardTransformation(lmax, 0, spatial, nonnegative);
         GPLSpec::detail::ExpandRealScalarCoefficients(lmax, nonnegative,
                                                       coefficients);
      }
      if (quantity != Quantity::SurfaceRadius)
         detail::RotateCoefficients(lmax, coefficients, rotated, matrices);
      const double scale = quantity == Quantity::Potential
                               ? model.PotentialNorm()
                               : quantity == Quantity::SensitivityKernel
                                     ? model.PotentialNorm() / model.MassNorm()
                                     : model.DensityNorm();
      base.quantity = quantity;
      base.representation = Representation::RotatedSlice;
      for (double phi : azimuths) {
         Sample sample = base;
         sample.theta_rad = std::numbers::pi / 2.0;
         sample.phi_rad = phi;
         sample.reference_x_m = 0.0;
         sample.reference_y_m = 0.0;
         sample.reference_z_m = 0.0;
         const auto ref_r = base.reference_radius_m;
         const double map_value = detail::EvaluateHarmonicSeries(
             lmax, rotated_mapping, harmonics, phi).real();
         detail::SetCoordinates(sample,
             (mesh.NodeRadius(static_cast<int>(element), node) + map_value) *
                 model.LengthNorm(),
             sample.theta_rad, phi);
         sample.reference_radius_m = ref_r;
         sample.value = quantity == Quantity::SurfaceRadius
                            ? sample.physical_radius_m
                            : scale * detail::EvaluateHarmonicSeries(
                                          lmax, rotated, harmonics, phi);
         callback(sample);
      }
   };

   for (std::size_t e = 0; e < static_cast<std::size_t>(mesh.NumberOfElements()); ++e) {
      emit(e, 0, detail::AtLowerTrace(model, e, interfaces));
      emit(e, mesh.PolynomialOrder(), detail::AtUpperTrace(model, e, interfaces));
   }
}

template <class Callback>
void ForEachRotatedSliceSample(
    const GeneralEarthModels::Density3D &model,
    const Eigen::VectorXcd &flat_solution, Quantity quantity,
    const Rotation &rotation, std::span<const double> azimuths,
    Callback &&callback) {
   ForEachRotatedSliceSample(
       model, std::span<const std::complex<double>>(flat_solution.data(),
           static_cast<std::size_t>(flat_solution.size())), quantity, rotation,
       azimuths, std::forward<Callback>(callback));
}

}   // namespace GPLSpec::FieldEvaluation

#endif
