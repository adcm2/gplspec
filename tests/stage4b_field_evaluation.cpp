#include "../gplspec/All"
#include "../gplspec/src/GeneralModels/Field_Bundle.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char *message) {
   if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const char *message) {
   if (std::abs(actual - expected) > tolerance * std::max(1.0, std::abs(expected)))
      throw std::runtime_error(message);
}

struct TwoLayerSource {
   int NumberOfLayers() const { return 2; }
   double LowerRadius(int layer) const { return layer == 0 ? 0.0 : 0.5; }
   double UpperRadius(int layer) const { return layer == 0 ? 0.5 : 1.0; }
   double OuterRadius() const { return 1.0; }
   double LengthNorm() const { return 2.0; }
   double MassNorm() const { return 64.0; }
   double TimeNorm() const { return 1.0; }
   auto Density(int layer) const {
      const double value = layer == 0 ? 4.0 : 10.0;
      return [value](double) { return value; };
   }
};

struct RadialExpansion {
   auto RadialMapping(int) const {
      return [](double radius, double, double) { return 0.1 * radius; };
   }
};

struct AngularRadialExpansion {
   auto RadialMapping(int) const {
      return [](double radius, double theta, double phi) {
         return 0.1 * radius * (1.0 + 0.2 * std::sin(theta) * std::cos(phi));
      };
   }
};
}   // namespace

int main(int argc, char **argv) {
   const auto high_degree = GPLSpec::FieldEvaluation::detail::HarmonicsAt(
       std::numbers::pi / 2.0, 101);
   auto highest_order_value = [](int l) {
      const double log_amplitude = 0.5 *
          (std::log((2.0 * l + 1.0) / (4.0 * std::numbers::pi)) +
           std::lgamma(2.0 * l + 1.0) - 2.0 * std::lgamma(l + 1.0) -
           2.0 * l * std::log(2.0));
      const double sign = l % 2 == 0 ? 1.0 : -1.0;
      return sign * std::exp(log_amplitude);
   };
   near(high_degree[100 * 100 + 100 + 100].real(),
        highest_order_value(100), 2.0e-12,
        "degree-100 highest-order equatorial harmonic underflowed");
   near(high_degree[101 * 101 + 101 + 101].real(),
        highest_order_value(101), 2.0e-12,
        "positive highest order lost its Condon-Shortley sign");
   near(high_degree[101 * 101 + 101 - 101].real(),
        -highest_order_value(101), 2.0e-12,
        "negative highest order lost its harmonic sign relation");

   using GeneralEarthModels::Density3D;
   using namespace GPLSpec::FieldEvaluation;
   const PlanetaryModel::TomographyZeroModel tomography;
   const TwoLayerSource source;
   const RadialExpansion mapping;
   const Density3D model(source, tomography, mapping, 2, 3, 0.25, 1.2);
   const auto &mesh = model.Node_Information();
   const int lmax = model.GSH_Grid().MaxDegree();
   const std::size_t stride = static_cast<std::size_t>((lmax + 1) * (lmax + 1));
   std::vector<std::complex<double>> solution(
       static_cast<std::size_t>(mesh.NumberOfGlobalRadialNodes()) * stride);
   for (int node = 0; node < mesh.NumberOfGlobalRadialNodes(); ++node) {
      solution[static_cast<std::size_t>(node) * stride] = 2.0;
      // A real l=1, |m|=1 field, used to check a nonconstant rotated slice.
      solution[static_cast<std::size_t>(node) * stride + 1] = -1.0;
      solution[static_cast<std::size_t>(node) * stride + 3] = 1.0;
   }

   std::vector<Sample> potential;
   ForEachAngularGridSample(model, solution, Quantity::Potential,
                            [&](const Sample &sample) { potential.push_back(sample); });
   const double y00 = 1.0 / std::sqrt(4.0 * std::numbers::pi);
   const double y11scale = std::sqrt(3.0 / (8.0 * std::numbers::pi));
   for (const auto &sample : potential) {
      if (sample.representation != Representation::AngularGrid) continue;
      const double expected = model.PotentialNorm() *
          (2.0 * y00 - 2.0 * y11scale * std::sin(sample.theta_rad) *
                              std::cos(sample.phi_rad));
      near(sample.value.real(), expected, 2.0e-12,
           "angular evaluator changed the orthonormal harmonic field");
      near(sample.value.imag(), 0.0, 2.0e-12,
           "real harmonic fixture acquired an imaginary field");
   }

   std::vector<Sample> body_vacuum;
   for (const auto &sample : potential) {
      if (sample.boundary == BoundaryKind::BodyVacuum &&
          std::abs(sample.reference_radius_m - 2.0) < 1.0e-12)
         body_vacuum.push_back(sample);
   }
   const std::size_t angular_count = model.GSH_Grid().Longitudes().size() *
                                     model.GSH_Grid().CoLatitudes().size();
   std::vector<Sample> outer_cohort, inner_cohort;
   ForEachAngularGridSample(model, solution, Quantity::Potential,
                            [&](const Sample &sample) { outer_cohort.push_back(sample); },
                            nullptr, TraceSelection::Outer);
   ForEachAngularGridSample(model, solution, Quantity::Potential,
                            [&](const Sample &sample) { inner_cohort.push_back(sample); },
                            nullptr, TraceSelection::Inner);
   require(outer_cohort.size() == (mesh.NumberOfElements() + 1) * angular_count &&
               inner_cohort.size() == (mesh.NumberOfElements() + 1) * angular_count,
           "single-side angular cohorts must retain one endpoint per boundary");
   require(potential.size() == 2 * mesh.NumberOfElements() * angular_count,
           "default angular trace selection must retain both endpoint sides");
   require(outer_cohort.front().side == TraceSide::None &&
               inner_cohort.back().side == TraceSide::None,
           "matched single-side cohorts must retain their external endpoint");
   require(body_vacuum.size() == 2 * angular_count,
           "body-vacuum surface must emit both material traces");
   for (std::size_t i = 0; i < angular_count; ++i) {
      const auto &inner = body_vacuum[i];
      const auto &outer = body_vacuum[i + angular_count];
      require(inner.side == TraceSide::Inner && outer.side == TraceSide::Outer,
              "body-vacuum trace sides are not explicit");
      require(inner.inner_element == outer.inner_element &&
                  inner.outer_element == outer.outer_element,
              "body-vacuum adjacent element IDs disagree");
      require(inner.material != outer.material,
              "body-vacuum material IDs must distinguish density traces");
      near(inner.value.real(), outer.value.real(), 0.0,
           "shared potential DOF must remain continuous across density traces");
      near(inner.physical_radius_m, 1.1 * inner.reference_radius_m, 2.0e-12,
           "mapped physical radius is incorrect");
      near(inner.reference_x_m,
           inner.reference_radius_m * std::sin(inner.theta_rad) *
               std::cos(inner.phi_rad),
           2.0e-12, "reference coordinates are incorrect");
   }

   std::vector<Sample> material_traces;
   for (const auto &sample : potential) {
      if (sample.boundary == BoundaryKind::Material &&
          std::abs(sample.reference_radius_m - 1.0) < 1.0e-12)
         material_traces.push_back(sample);
   }
   require(material_traces.size() == 2 * angular_count,
           "material interface must emit both element-local traces");
   for (std::size_t i = 0; i < angular_count; ++i) {
      const auto &inner = material_traces[i];
      const auto &outer = material_traces[i + angular_count];
      require(inner.side == TraceSide::Inner && outer.side == TraceSide::Outer &&
                  inner.material != outer.material,
              "material interface IDs or sides are incorrect");
      near(inner.value.real(), outer.value.real(), 0.0,
           "potential must use one shared degree of freedom at material jumps");
   }

   std::vector<Sample> density;
   ForEachAngularGridSample(model, std::span<const std::complex<double>>{},
                            Quantity::PhysicalDensity,
                            [&](const Sample &sample) { density.push_back(sample); });
   std::vector<Sample> referential_density;
   ForEachAngularGridSample(model, std::span<const std::complex<double>>{},
                            Quantity::ReferentialDensity,
                            [&](const Sample &sample) {
                               referential_density.push_back(sample);
                            });
   require(density.size() == referential_density.size(),
           "physical and referential density sample layouts differ");
   for (std::size_t i = 0; i < density.size(); ++i) {
      const int element = static_cast<int>(density[i].element);
      const int node = density[i].radial_node;
      const double jacobian = model.Jacobian_Point(element, node, 0);
      near(density[i].value.real() * jacobian,
           referential_density[i].value.real(), 2.0e-12,
           "density views must preserve rho_physical*J=rho_referential");
   }
   bool saw_nonunit_jacobian = false;
   bool saw_inner_material = false;
   bool saw_outer_material = false;
   for (const auto &sample : density) {
      if (sample.boundary != BoundaryKind::BodyVacuum) continue;
      const std::size_t element = sample.element;
      const int node = sample.radial_node;
      const double rho_ref = model.Density_Point(static_cast<int>(element), node, 0);
      const double jacobian = model.Jacobian_Point(static_cast<int>(element), node, 0);
      if (std::abs(jacobian - 1.0) > 1.0e-6) saw_nonunit_jacobian = true;
      near(sample.value.real() * jacobian, rho_ref * model.DensityNorm(), 2.0e-12,
           "physical density must equal referential density divided by J");
      if (rho_ref > 0.0) saw_inner_material = true;
      if (rho_ref == 0.0) saw_outer_material = true;
   }
   require(saw_nonunit_jacobian && saw_inner_material && saw_outer_material,
           "fixture must exercise nonunit J and both one-sided densities");

   const Rotation z_quarter_turn{0.0, 0.0, std::numbers::pi / 2.0};
   std::vector<Sample> rotated;
   ForEachAngularGridSample(model, solution, Quantity::Potential,
                            [&](const Sample &sample) { rotated.push_back(sample); },
                            &z_quarter_turn);
   for (const auto &sample : rotated) {
      const double expected = model.PotentialNorm() *
          (2.0 * y00 + 2.0 * y11scale * std::sin(sample.theta_rad) *
                              std::sin(sample.phi_rad));
      near(sample.value.real(), expected, 2.0e-11,
           "rotated nonconstant harmonic field has the wrong orientation");
   }

   // Exercise nonconstant real mapping and physical density through the
   // complex-output inverse expected for full rotated coefficients.
   const AngularRadialExpansion angular_mapping;
   const Density3D angular_model(source, tomography, angular_mapping, 2, 3,
                                 0.25, 1.2);
   const auto &angular_grid = angular_model.GSH_Grid();
   const int angular_lmax = angular_grid.MaxDegree();
   const std::size_t angular_stride =
       static_cast<std::size_t>((angular_lmax + 1) * (angular_lmax + 1));
   const std::size_t spatial_count = angular_grid.Longitudes().size() *
                                     angular_grid.CoLatitudes().size();
   const auto mapping_rotation = GPLSpec::detail::RotationOutputMatrices(
       angular_lmax, z_quarter_turn.alpha, z_quarter_turn.beta,
       z_quarter_turn.gamma);
   auto complex_rotated_real_coefficients = [&](const std::vector<double> &field) {
      std::vector<std::complex<double>> nonnegative(
          static_cast<std::size_t>((angular_lmax + 1) * (angular_lmax + 2) / 2));
      std::vector<std::complex<double>> coefficients(angular_stride),
          rotated_coefficients(angular_stride);
      angular_grid.ForwardTransformation(angular_lmax, 0, field, nonnegative);
      GPLSpec::detail::ExpandRealScalarCoefficients(
          angular_lmax, nonnegative, coefficients);
      GPLSpec::FieldEvaluation::detail::RotateCoefficients(
          angular_lmax, coefficients, rotated_coefficients, mapping_rotation);
      return rotated_coefficients;
   };
   auto complex_rotated_real_field = [&](const std::vector<double> &field) {
      const auto coefficients = complex_rotated_real_coefficients(field);
      std::vector<std::complex<double>> rotated_field(spatial_count);
      angular_grid.InverseTransformation(angular_lmax, 0, coefficients,
                                         rotated_field);
      return rotated_field;
   };
   auto spatial_index = [&](double theta, double phi) {
      const auto theta_it = std::find_if(
          angular_grid.CoLatitudes().begin(), angular_grid.CoLatitudes().end(),
          [&](double value) { return std::abs(value - theta) < 1.0e-13; });
      const auto phi_it = std::find_if(
          angular_grid.Longitudes().begin(), angular_grid.Longitudes().end(),
          [&](double value) { return std::abs(value - phi) < 1.0e-13; });
      require(theta_it != angular_grid.CoLatitudes().end() &&
                  phi_it != angular_grid.Longitudes().end(),
              "rotated sample is not located on its declared angular grid");
      return static_cast<std::size_t>(theta_it - angular_grid.CoLatitudes().begin()) *
                 angular_grid.Longitudes().size() +
             static_cast<std::size_t>(phi_it - angular_grid.Longitudes().begin());
   };

   std::vector<Sample> surface_unrotated, surface_rotated;
   ForEachAngularGridSample(
       angular_model, std::span<const std::complex<double>>{},
       Quantity::SurfaceRadius,
       [&](const Sample &sample) { surface_unrotated.push_back(sample); });
   ForEachAngularGridSample(
       angular_model, std::span<const std::complex<double>>{},
       Quantity::SurfaceRadius,
       [&](const Sample &sample) { surface_rotated.push_back(sample); },
       &z_quarter_turn);
   require(surface_rotated.size() == surface_unrotated.size(),
           "rotated and unrotated mapping samples differ in size");
   const auto [surface_min, surface_max] = std::minmax_element(
       surface_unrotated.begin(), surface_unrotated.end(),
       [](const Sample &a, const Sample &b) {
          return a.value.real() < b.value.real();
       });
   require(surface_max->value.real() - surface_min->value.real() > 1.0e-6,
           "mapping regression must include angular variation");
   std::map<std::pair<int, int>, std::vector<std::complex<double>>> rotated_mapping;
   for (const auto &sample : surface_rotated) {
      const auto key = std::pair{static_cast<int>(sample.element), sample.radial_node};
      if (!rotated_mapping.contains(key)) {
         const auto mapping = angular_model.MappingAtRadialNode(key.first, key.second);
         rotated_mapping.emplace(key, complex_rotated_real_field(mapping));
      }
      const std::size_t k = spatial_index(sample.theta_rad, sample.phi_rad);
      const double expected_radius =
          (angular_model.Node_Information().NodeRadius(key.first, key.second) +
           rotated_mapping.at(key)[k].real()) * angular_model.LengthNorm();
      const double coordinate_tolerance =
          1.0e-12 * std::abs(expected_radius) +
          1.0e-14 * std::max(1.0, std::abs(sample.reference_radius_m));
      require(std::abs(sample.value.real() - expected_radius) <=
                  coordinate_tolerance,
              "rotated mapping differs from complex-output inverse expectation");
   }

   std::vector<Sample> density_rotated;
   ForEachAngularGridSample(
       angular_model, std::span<const std::complex<double>>{},
       Quantity::PhysicalDensity,
       [&](const Sample &sample) { density_rotated.push_back(sample); },
       &z_quarter_turn);
   std::map<std::pair<int, int>, std::vector<std::complex<double>>> rotated_density;
   for (const auto &sample : density_rotated) {
      const auto key = std::pair{static_cast<int>(sample.element), sample.radial_node};
      if (!rotated_density.contains(key)) {
         std::vector<double> physical_density(spatial_count);
         for (std::size_t k = 0; k < spatial_count; ++k) {
            physical_density[k] =
                angular_model.Density_Point(key.first, key.second,
                                            static_cast<int>(k)) /
                angular_model.Jacobian_Point(key.first, key.second,
                                             static_cast<int>(k));
         }
         rotated_density.emplace(key,
                                 complex_rotated_real_field(physical_density));
      }
      const std::size_t k = spatial_index(sample.theta_rad, sample.phi_rad);
      const double expected_density =
          angular_model.DensityNorm() * rotated_density.at(key)[k].real();
      const double density_tolerance =
          1.0e-12 * std::abs(expected_density) +
          1.0e-14 * angular_model.DensityNorm();
      require(std::abs(sample.value.real() - expected_density) <=
                  density_tolerance,
              "rotated physical density differs from complex-output inverse expectation");
   }

   std::vector<Sample> radial_density;
   ForEachRadialHarmonicSample(
       angular_model, std::span<const std::complex<double>>{},
       Quantity::PhysicalDensity,
       [&](const Sample &sample) { radial_density.push_back(sample); });
   std::map<std::pair<int, int>, std::vector<std::complex<double>>>
       radial_density_coefficients;
   for (const auto &sample : radial_density) {
      const auto key = std::pair{static_cast<int>(sample.element), sample.radial_node};
      if (!radial_density_coefficients.contains(key)) {
         std::vector<double> physical_density(spatial_count);
         for (std::size_t k = 0; k < spatial_count; ++k)
            physical_density[k] = angular_model.Density_Point(
                key.first, key.second, static_cast<int>(k)) /
                angular_model.Jacobian_Point(key.first, key.second,
                                             static_cast<int>(k));
         std::vector<std::complex<double>> nonnegative(
             static_cast<std::size_t>((angular_lmax + 1) * (angular_lmax + 2) / 2));
         std::vector<std::complex<double>> coefficients(angular_stride);
         angular_grid.ForwardTransformation(angular_lmax, 0, physical_density,
                                            nonnegative);
         GPLSpec::detail::ExpandRealScalarCoefficients(
             angular_lmax, nonnegative, coefficients);
         radial_density_coefficients.emplace(key, std::move(coefficients));
      }
      const std::size_t coefficient = static_cast<std::size_t>(
          sample.degree * sample.degree + sample.degree + sample.order);
      const double expected_density = angular_model.DensityNorm() *
          radial_density_coefficients.at(key)[coefficient].real();
      const double tolerance = 1.0e-12 * std::abs(expected_density) +
                               1.0e-14 * angular_model.DensityNorm();
      require(std::abs(sample.value.real() - expected_density) <= tolerance,
              "radial physical-density harmonics disagree with fresh transform");
   }

   const std::vector<double> slice_azimuths{-1.0, 0.0, 1.0};
   std::vector<Sample> rotated_density_slice;
   ForEachRotatedSliceSample(
       angular_model, std::span<const std::complex<double>>{},
       Quantity::PhysicalDensity, z_quarter_turn, slice_azimuths,
       [&](const Sample &sample) { rotated_density_slice.push_back(sample); });
   std::map<std::pair<int, int>, std::vector<std::complex<double>>>
       slice_density_coefficients, slice_mapping_coefficients;
   std::vector<std::complex<double>> slice_harmonics =
       GPLSpec::FieldEvaluation::detail::HarmonicsAt(
           std::numbers::pi / 2.0, angular_lmax);
   for (const auto &sample : rotated_density_slice) {
      const auto key = std::pair{static_cast<int>(sample.element), sample.radial_node};
      if (!slice_density_coefficients.contains(key)) {
         std::vector<double> physical_density(spatial_count), mapping(spatial_count);
         const auto mapping_node = angular_model.MappingAtRadialNode(key.first, key.second);
         std::copy(mapping_node.begin(), mapping_node.end(), mapping.begin());
         for (std::size_t k = 0; k < spatial_count; ++k)
            physical_density[k] = angular_model.Density_Point(
                key.first, key.second, static_cast<int>(k)) /
                angular_model.Jacobian_Point(key.first, key.second,
                                             static_cast<int>(k));
         slice_mapping_coefficients.emplace(
             key, complex_rotated_real_coefficients(mapping));
         slice_density_coefficients.emplace(
             key, complex_rotated_real_coefficients(physical_density));
      }
      const double expected_map =
          GPLSpec::FieldEvaluation::detail::EvaluateHarmonicSeries(
              angular_lmax, slice_mapping_coefficients.at(key),
              slice_harmonics, sample.phi_rad).real();
      const double expected_radius =
          (angular_model.Node_Information().NodeRadius(key.first, key.second) +
           expected_map) * angular_model.LengthNorm();
      const double coordinate_tolerance =
          1.0e-12 * std::abs(expected_radius) +
          1.0e-14 * std::max(1.0, std::abs(sample.reference_radius_m));
      require(std::abs(sample.physical_radius_m - expected_radius) <=
                  coordinate_tolerance,
              "rotated-slice mapping differs from complex inverse expectation");
      const double expected_density = angular_model.DensityNorm() *
          GPLSpec::FieldEvaluation::detail::EvaluateHarmonicSeries(
              angular_lmax, slice_density_coefficients.at(key),
              slice_harmonics, sample.phi_rad).real();
      const double density_tolerance =
          1.0e-12 * std::abs(expected_density) +
          1.0e-14 * angular_model.DensityNorm();
      require(std::abs(sample.value.real() - expected_density) <=
                  density_tolerance,
              "rotated-slice density differs from complex inverse expectation");
   }

   const auto colatitudes = model.GSH_Grid().CoLatitudes();
   require(std::none_of(colatitudes.begin(), colatitudes.end(), [](double theta) {
                           return std::abs(theta - std::numbers::pi / 2.0) < 1.0e-14;
                        }),
           "slice fixture must require an off-grid equatorial latitude");
   std::vector<std::complex<double>> axial_solution(solution.size());
   for (int node = 0; node < mesh.NumberOfGlobalRadialNodes(); ++node)
      axial_solution[static_cast<std::size_t>(node) * stride + 2] = 1.0;
   const Rotation tilt_x{0.0, std::numbers::pi / 2.0, 0.0};
   const std::vector<double> azimuths{-1.0, 0.0, 1.0};
   std::vector<Sample> slice;
   ForEachRotatedSliceSample(model, axial_solution, Quantity::Potential, tilt_x,
                             azimuths,
                             [&](const Sample &sample) { slice.push_back(sample); });
   require(slice.size() == static_cast<std::size_t>(mesh.NumberOfElements()) * 6,
           "equatorial slice must emit both endpoint traces at each azimuth");
   const double y10scale = std::sqrt(3.0 / (4.0 * std::numbers::pi));
   for (const auto &sample : slice) {
      require(sample.representation == Representation::RotatedSlice,
              "off-grid angular samples must be labeled as a rotated slice");
      near(sample.theta_rad, std::numbers::pi / 2.0, 0.0,
           "rotated slice substituted a grid latitude for the requested equator");
      near(sample.value.real(), -model.PotentialNorm() * y10scale *
                                    std::cos(sample.phi_rad),
           2.0e-11,
           "tilted axial harmonic did not rotate to the analytic negative x-axis field");
   }

   std::vector<Sample> unrotated_slice, explicit_identity_slice;
   ForEachEquatorialSliceSample(
       model, axial_solution, Quantity::Potential, azimuths,
       [&](const Sample &sample) { unrotated_slice.push_back(sample); });
   ForEachRotatedSliceSample(
       model, axial_solution, Quantity::Potential, Rotation{}, azimuths,
       [&](const Sample &sample) { explicit_identity_slice.push_back(sample); });
   require(unrotated_slice.size() == explicit_identity_slice.size(),
           "unrotated and explicit Euler slices must share the sample cohort");
   for (std::size_t i = 0; i < unrotated_slice.size(); ++i) {
      const auto &plain = unrotated_slice[i];
      const auto &legacy = explicit_identity_slice[i];
      std::vector<double> mapping_values(
          model.GSH_Grid().Longitudes().size() *
              model.GSH_Grid().CoLatitudes().size());
      const auto mapping_node = model.MappingAtRadialNode(
          static_cast<int>(plain.element), plain.radial_node);
      std::copy(mapping_node.begin(), mapping_node.end(), mapping_values.begin());
      const int mapping_lmax = model.GSH_Grid().MaxDegree();
      std::vector<std::complex<double>> mapping_nonnegative(
          static_cast<std::size_t>((mapping_lmax + 1) * (mapping_lmax + 2) / 2));
      std::vector<std::complex<double>> mapping_coefficients(
          static_cast<std::size_t>((mapping_lmax + 1) * (mapping_lmax + 1)));
      model.GSH_Grid().ForwardTransformation(mapping_lmax, 0, mapping_values,
                                             mapping_nonnegative);
      GPLSpec::detail::ExpandRealScalarCoefficients(
          mapping_lmax, mapping_nonnegative, mapping_coefficients);
      const double mapping_value =
          GPLSpec::FieldEvaluation::detail::EvaluateHarmonicSeries(
              mapping_lmax, mapping_coefficients,
              GPLSpec::FieldEvaluation::detail::HarmonicsAt(
                  std::numbers::pi / 2.0, model.GSH_Grid().MaxDegree()),
              plain.phi_rad).real();
      const double expected_plain_radius =
          (mesh.NodeRadius(static_cast<int>(plain.element), plain.radial_node) +
           mapping_value) * model.LengthNorm();
      near(plain.physical_radius_m, expected_plain_radius, 1.0e-12,
           "unrotated slice mapping differs from its unrotated coefficient field");
      const double legacy_mapping_factor =
          std::sqrt(3.1415926535 / std::numbers::pi);
      const double expected_legacy_radius =
          (mesh.NodeRadius(static_cast<int>(legacy.element), legacy.radial_node) +
           mapping_value * legacy_mapping_factor) * model.LengthNorm();
      near(legacy.physical_radius_m, expected_legacy_radius, 1.0e-12,
           "explicit Euler slice changed its historical mapping normalization");
   }

   std::vector<Sample> radial;
   ForEachRadialHarmonicSample(model, solution, Quantity::Potential,
                               [&](const Sample &sample) { radial.push_back(sample); });
   require(!radial.empty(), "radial harmonic profile emitted no samples");
   for (const auto &sample : radial) {
      require(sample.representation == Representation::HarmonicCoefficients &&
                  sample.degree >= 0 && std::abs(sample.order) <= sample.degree,
              "harmonic sample metadata is invalid");
      const std::size_t global = static_cast<std::size_t>(
          mesh.GlobalNodeIndex(static_cast<int>(sample.element), sample.radial_node));
      const std::size_t coefficient = static_cast<std::size_t>(
          sample.degree * sample.degree + sample.degree + sample.order);
      near(sample.value.real(),
           model.PotentialNorm() *
               solution[global * stride + coefficient].real(),
           0.0, "radial harmonic evaluator selected the wrong flat coefficient");
   }

   std::vector<Sample> sensitivity;
   ForEachAngularGridSample(model, solution, Quantity::SensitivityKernel,
                            [&](const Sample &sample) {
                               sensitivity.push_back(sample);
                            });
   require(sensitivity.size() == potential.size(),
           "sensitivity and potential spatial sample layouts differ");
   for (std::size_t i = 0; i < potential.size(); ++i) {
      near(sensitivity[i].value.real(),
           potential[i].value.real() / model.MassNorm(), 2.0e-12,
           "sensitivity samples must retain PotentialNorm/MassNorm units");
   }

   std::vector<Sample> degree_power;
   ForEachRadialHarmonicSample(model, solution, Quantity::DegreePower,
                               [&](const Sample &sample) {
                                  degree_power.push_back(sample);
                               });
   require(!degree_power.empty() &&
               degree_power.front().representation == Representation::DegreePower,
           "degree power must use its own dimensionless representation");
   require(degree_power.size() == static_cast<std::size_t>(mesh.NumberOfElements()) *
                                      (model.Poly_Order() + 1) * (lmax + 1),
           "degree power must include every element-local radial node");
   bool saw_interior_power_node = false;
   for (const auto &sample : degree_power) {
      if (sample.radial_node > 0 && sample.radial_node < model.Poly_Order()) {
         saw_interior_power_node = true;
         require(sample.side == TraceSide::None &&
                     sample.boundary == BoundaryKind::None,
                 "interior power samples must not claim an interface side");
      }
      const double expected = sample.degree == 0 ? 1.0
                              : sample.degree == 1 ? 0.5
                                                   : 0.0;
      near(sample.value.real(), expected, 1.0e-14,
           "degree power normalization differs from sum_m |a_lm|^2/|a_00|^2");
   }
   require(saw_interior_power_node,
           "degree power must retain interior local-node metadata");

   std::vector<Sample> surface_radius;
   ForEachAngularGridSample(model, std::span<const std::complex<double>>{},
                            Quantity::SurfaceRadius,
                            [&](const Sample &sample) {
                               surface_radius.push_back(sample);
                            });
   require(surface_radius.size() == 2 * angular_count,
           "surface radius must emit both modeled surface traces only");
   for (const auto &sample : surface_radius) {
      require(detail::IsSurfaceTrace(sample),
              "surface-radius evaluator emitted a non-surface sample");
      near(sample.value.real(), sample.physical_radius_m, 2.0e-12,
           "surface radius value is not the mapped physical radius");
   }

   const Radial_Tools::RadialMesh surface_mesh(source, model.q(), 0.25, 1.0,
                                               false);
   const std::size_t last =
       static_cast<std::size_t>(surface_mesh.NumberOfElements() - 1);
   const auto surface_trace = detail::SurfaceOnlyTrace(surface_mesh, last,
                                                        model.LengthNorm());
   require(surface_trace.boundary == BoundaryKind::SurfaceWithoutExterior &&
               surface_trace.side == TraceSide::Inner &&
               surface_trace.inner_element == last &&
               surface_trace.outer_element == NoIndex &&
               surface_trace.outer_material == -1,
           "surface-only mesh must expose only its modeled inner trace");

   if (argc == 2) {
      const std::filesystem::path output = argv[1];
      GPLSpec::FieldBundle::RunMetadata metadata;
      metadata.workflow = "stage4b_two_layer_visual_fixture";
      metadata.dataset = "two_layer_field";
      std::vector<double> slice_azimuths;
      for (double phi : angular_grid.Longitudes())
         slice_azimuths.push_back(phi);
      GPLSpec::FieldBundle::WriteAngular(
          output / "potential", model, solution, Quantity::Potential, metadata);
      GPLSpec::FieldBundle::WriteAngular(
          output / "density_physical", model,
          std::span<const std::complex<double>>{}, Quantity::PhysicalDensity,
          metadata);
      GPLSpec::FieldBundle::WriteAngular(
          output / "density_referential", model,
          std::span<const std::complex<double>>{}, Quantity::ReferentialDensity,
          metadata);
      GPLSpec::FieldBundle::WriteRadialHarmonics(
          output / "degree_power", model, solution, Quantity::DegreePower,
          metadata);
      GPLSpec::FieldBundle::WriteSlice(
          output / "potential_slice", model, solution, Quantity::Potential,
          z_quarter_turn, slice_azimuths, metadata);
      GPLSpec::FieldBundle::WriteSlice(
          output / "density_slice", model,
          std::span<const std::complex<double>>{}, Quantity::PhysicalDensity,
          z_quarter_turn, slice_azimuths, metadata);
   } else if (argc != 1) {
      throw std::runtime_error("usage: stage4b_field_evaluation [fixture-output-directory]");
   }
   return 0;
}
