#include "../gplspec/All"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
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
}   // namespace

int main() {
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
   for (const auto &sample : degree_power) {
      const double expected = sample.degree == 0 ? 1.0
                              : sample.degree == 1 ? 0.5
                                                   : 0.0;
      near(sample.value.real(), expected, 1.0e-14,
           "degree power normalization differs from sum_m |a_lm|^2/|a_00|^2");
   }

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

   return 0;
}
