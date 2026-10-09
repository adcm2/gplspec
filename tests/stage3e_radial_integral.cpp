#include <gplspec/All>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Complex = std::complex<double>;
constexpr double pi = std::numbers::pi_v<double>;
constexpr double sqrt_four_pi = 3.544907701811032054596334966682290366;

void require(bool condition, const std::string &message) {
   if (!condition) throw std::runtime_error(message);
}

void require_analytic(double actual, double expected, double scale,
                      const std::string &message) {
   const double tolerance = 1e-11 * std::abs(expected) + 1e-13 * scale;
   require(std::abs(actual - expected) <= tolerance, message);
}

double homogeneous_scalar_potential(double G, double density, double body_radius,
                                    double radius) {
   if (radius <= body_radius) {
      return 2.0 * pi * G * density / 3.0 *
             (radius * radius - 3.0 * body_radius * body_radius);
   }
   return -4.0 * pi * G * density * std::pow(body_radius, 3.0) /
          (3.0 * radius);
}

double layered_scalar_potential(const std::vector<double> &radii,
                                const std::vector<double> &density,
                                double G, double radius) {
   double inner_moment = 0.0;
   double outer_moment = 0.0;
   double lower = 0.0;
   for (std::size_t layer = 0; layer < radii.size(); ++layer) {
      const double upper = radii[layer];
      const double interior_upper = std::min(radius, upper);
      if (interior_upper > lower) {
         inner_moment += density[layer] / 3.0 *
                         (std::pow(interior_upper, 3.0) - std::pow(lower, 3.0));
      }
      const double exterior_lower = std::max(radius, lower);
      if (upper > exterior_lower) {
         outer_moment += density[layer] / 2.0 *
                         (upper * upper - exterior_lower * exterior_lower);
      }
      lower = upper;
   }
   const double interior_term = radius == 0.0 ? 0.0 : inner_moment / radius;
   return -4.0 * pi * G * (interior_term + outer_moment);
}

struct R2YlmTomography {
   double body_radius;
   double outer_radius;
   double length_norm;
   int order;

   double GetValueAt(double depth, double longitude, double latitude,
                     bool = false) const {
      const double radius = body_radius - depth / (length_norm / 1000.0);
      const double cos_theta = std::sin(latitude * pi / 180.0);
      const double sin_theta_squared = 1.0 - cos_theta * cos_theta;
      const double angular =
          order == 0
              ? std::sqrt(5.0 / (16.0 * pi)) *
                    (3.0 * cos_theta * cos_theta - 1.0)
              : std::sqrt(15.0 / (32.0 * pi)) * sin_theta_squared *
                    std::cos(2.0 * longitude * pi / 180.0);
      return (radius / body_radius) * (radius / body_radius) * angular;
   }
   std::vector<double> GetDepths() const { return {-1.0, 10000.0}; }
   std::vector<double> GetLatitudes() const { return {-90.0, 90.0}; }
   std::vector<double> GetLongitudes() const { return {0.0, 360.0}; }
   std::vector<double> GetValues() const { return {0.0}; }
};

double r2_y20_potential(double G, double rho0, double body_radius,
                        double radius) {
   constexpr int degree = 2;
   if (radius <= body_radius) {
      const double inner_outer =
          std::pow(radius, degree + 2) / (2.0 * degree + 3.0) +
          std::pow(radius, degree) *
              (body_radius * body_radius - radius * radius) / 2.0;
      return -4.0 * pi * G * rho0 /
             ((2.0 * degree + 1.0) * std::pow(body_radius, degree)) *
             inner_outer;
   }
   return -4.0 * pi * G * rho0 *
          std::pow(body_radius, degree + 3) /
          ((2.0 * degree + 1.0) * (2.0 * degree + 3.0) *
           std::pow(radius, degree + 1));
}

double volume_weighted_potential_rms(const GeneralEarthModels::Density3D &model,
                                     const Eigen::VectorXcd &solution,
                                     double body_radius) {
   const auto &mesh = model.Node_Information();
   const auto &q = model.qP();
   const std::size_t stride =
       static_cast<std::size_t>(model.GSH_Grid().MaxDegree() + 1) *
       (model.GSH_Grid().MaxDegree() + 1);
   double error_integral = 0.0;
   double volume_integral = 0.0;
   const double physical_scale = model.PotentialNorm();
   for (std::size_t element = 0; element < mesh.NumberOfElements(); ++element) {
      for (std::size_t node = 0; node < mesh.NumNodesWithinElement(); ++node) {
         const double radius = mesh.NodeRadius(element, node);
         const double weight = q.W(node) * mesh.ElementWidth(element) * 0.5 *
                               radius * radius;
         const std::size_t radial = mesh.GlobalNodeIndex(element, node);
         const double exact_y00 = sqrt_four_pi * homogeneous_scalar_potential(
                                                       model.GravitationalConstant(),
                                                       model.Density_Point(0, 0, 0),
                                                       body_radius, radius);
         double coefficient_error_norm2 = 0.0;
         for (std::size_t harmonic = 0; harmonic < stride; ++harmonic) {
            const Complex expected = harmonic == 0 ? Complex{exact_y00, 0.0}
                                                   : Complex{};
            coefficient_error_norm2 +=
                std::norm(solution[radial * stride + harmonic] - expected);
         }
         error_integral += weight * coefficient_error_norm2 *
                           physical_scale * physical_scale;
         volume_integral += weight;
      }
   }
   return std::sqrt(error_integral / (4.0 * pi * volume_integral));
}

}  // namespace

int main() {
   using GeneralEarthModels::Density3D;
   using GeneralEarthModels::SimpleModels::spherical_model;
   using Gravity_Tools::GravitationalSphericalIntegral;
   using Gravity_Tools::HomogeneousSphereIntegral;
   using Gravity_Tools::SolveGravitationalPotential;

   constexpr double body_radius = 0.63;
   constexpr double outer_radius = 1.2;
   const Density3D sphere = Density3D::SphericalHomogeneousPlanet(
       body_radius, 2.0, 6371000.0, 3600.0, 5.972e24, 0.2, outer_radius);
   const double G = sphere.GravitationalConstant();
   const double rho0 = sphere.Density_Point(0, 0, 0);
   const double potential_scale =
       std::abs(homogeneous_scalar_potential(G, rho0, body_radius, 0.0));

   // The arbitrary-radius overload is a scalar field; the mesh overload is
   // the Y00 coefficient. Include both sides of the nonunit sphere boundary.
   const std::vector<double> radii{0.0, 0.31, body_radius, 0.91, outer_radius,
                                   2.5};
   const auto scalar_values = HomogeneousSphereIntegral(sphere, radii);
   for (std::size_t i = 0; i < radii.size(); ++i) {
      require_analytic(scalar_values[i],
                       homogeneous_scalar_potential(G, rho0, body_radius,
                                                    radii[i]),
                       potential_scale,
                       "arbitrary-radius homogeneous sphere formula differs");
   }
   const auto coefficient_values = HomogeneousSphereIntegral(sphere);
   const auto &sphere_mesh = sphere.Node_Information();
   for (std::size_t radial = 0; radial < coefficient_values.size(); ++radial) {
      const double radius = radial == sphere_mesh.NumberOfElements()
                                ? sphere_mesh.OuterRadius()
                                : sphere_mesh.ElementLowerRadius(radial);
      require_analytic(coefficient_values[radial],
                       sqrt_four_pi * homogeneous_scalar_potential(
                                          G, rho0, body_radius, radius),
                       sqrt_four_pi * potential_scale,
                       "mesh-boundary Y00 sphere formula differs");
   }

   const std::vector<double> layer_radii{0.35, body_radius};
   const std::vector<double> layer_density{2.5, 1.2};
   const Density3D layered = Density3D::SphericalLayeredPlanet(
       layer_radii, layer_density, 6371000.0, 3600.0, 5.972e24, 0.2,
       outer_radius);
   const auto layered_integral = GravitationalSphericalIntegral(layered);
   const auto &layer_mesh = layered.Node_Information();
   double layered_scale = 0.0;
   std::vector<double> expected_layered(layered_integral.size());
   for (std::size_t radial = 0; radial < layered_integral.size(); ++radial) {
      const double radius = radial == layer_mesh.NumberOfElements()
                                ? layer_mesh.OuterRadius()
                                : layer_mesh.ElementLowerRadius(radial);
      expected_layered[radial] =
          sqrt_four_pi * layered_scalar_potential(layer_radii, layer_density,
                                                  layered.GravitationalConstant(),
                                                  radius);
      layered_scale = std::max(layered_scale, std::abs(expected_layered[radial]));
   }
   for (std::size_t radial = 0; radial < layered_integral.size(); ++radial) {
      require_analytic(layered_integral[radial][0].real(),
                       expected_layered[radial], layered_scale,
                       "layered monopole Green integral differs from shell formula");
      for (std::size_t harmonic = 1; harmonic < layered_integral[radial].size();
           ++harmonic) {
         require(std::abs(layered_integral[radial][harmonic]) <=
                     1e-13 * layered_scale,
                 "spherical layered model produced a nonmonopole integral");
      }
   }
   double expected_mass = 0.0;
   double lower = 0.0;
   for (std::size_t layer = 0; layer < layer_radii.size(); ++layer) {
      expected_mass += 4.0 * pi * layered.MassNorm() * layer_density[layer] /
                       3.0 *
                       (std::pow(layer_radii[layer], 3.0) -
                        std::pow(lower, 3.0));
      lower = layer_radii[layer];
   }
   require_analytic(layered.Mass(), expected_mass, std::abs(expected_mass),
                    "layered density mass differs from closed-form shells: actual=" +
                        std::to_string(layered.Mass()) + " expected=" +
                        std::to_string(expected_mass));

   // An axisymmetric positive density with rho_20(r) proportional to r^2
   // checks both Green-function powers independently of the monopole case.
   constexpr int harmonic_order = 4;
   constexpr int lmax = 2;
   const double length_norm = 6371000.0;
   const R2YlmTomography r2_tomography{body_radius, outer_radius, length_norm,
                                       0};
   const auto base = spherical_model::HomogeneousSphere(
       body_radius, 2.0, length_norm, 3600.0, 5.972e24);
   const Density3D angular = Density3D::SphericalThreeDimensionalPlanetFromModel(
       base, r2_tomography, harmonic_order, lmax, 0.2, outer_radius);
   const auto &angular_grid = angular.GSH_Grid();
   std::vector<double> spatial(angular_grid.NumberOfLongitudes() *
                               angular_grid.NumberOfCoLatitudes());
   std::vector<Complex> spatial_complex(spatial.size());
   std::vector<Complex> harmonic_coefficients((lmax + 1) * (lmax + 1));
   const std::size_t l2_m0 = 2 * 2 + 2;
   for (std::size_t element = 0; element < angular.Num_Elements(); ++element) {
      for (std::size_t node = 0; node < angular.Poly_Order() + 1; ++node) {
         const double radius = angular.Node_Information().NodeRadius(element, node);
         for (std::size_t sample = 0; sample < spatial.size(); ++sample) {
            spatial[sample] = angular.Density_Point(element, node, sample);
            spatial_complex[sample] = Complex{spatial[sample], 0.0};
         }
         std::fill(harmonic_coefficients.begin(), harmonic_coefficients.end(),
                   Complex{});
         angular_grid.ForwardTransformation(lmax, 0, spatial_complex,
                                            harmonic_coefficients);
         const double expected =
             0.005 * rho0 * (radius / body_radius) * (radius / body_radius);
         const bool interior_trace =
             angular.Node_Information().LayerNumber(element) <
             static_cast<int>(base.NumberOfLayers());
         if (interior_trace) {
            require_analytic(harmonic_coefficients[l2_m0].real(), expected,
                             std::max(1.0, std::abs(expected)),
                             "analytic r^2 Y20 density coefficient differs: actual_micro=" +
                                 std::to_string(1e6 * harmonic_coefficients[l2_m0].real()) +
                                 " expected_micro=" + std::to_string(1e6 * expected) +
                                 " radius=" + std::to_string(radius) +
                                 " sample_min=" +
                                 std::to_string(*std::min_element(spatial.begin(),
                                                                  spatial.end())) +
                                 " sample_max=" +
                                 std::to_string(*std::max_element(spatial.begin(),
                                                                  spatial.end())) +
                                 " c4_micro=" + std::to_string(1e6 * harmonic_coefficients[4].real()) +
                                 " c5_micro=" + std::to_string(1e6 * harmonic_coefficients[5].real()) +
                                 " c6_micro=" + std::to_string(1e6 * harmonic_coefficients[6].real()) +
                                 " c7_micro=" + std::to_string(1e6 * harmonic_coefficients[7].real()) +
                                 " c8_micro=" + std::to_string(1e6 * harmonic_coefficients[8].real()));
         } else {
            require(std::abs(harmonic_coefficients[l2_m0]) <= 1e-13,
                    "vacuum exterior trace acquired a density harmonic");
         }
      }
   }
   const auto angular_integral = GravitationalSphericalIntegral(angular);
   const auto &angular_mesh = angular.Node_Information();
   double harmonic_scale = 0.0;
   std::vector<double> expected_harmonic(angular_integral.size());
   const double rho20_scale = 0.005 * rho0;
   for (std::size_t radial = 0; radial < angular_integral.size(); ++radial) {
      const double radius = radial == angular_mesh.NumberOfElements()
                                ? angular_mesh.OuterRadius()
                                : angular_mesh.ElementLowerRadius(radial);
      expected_harmonic[radial] =
          r2_y20_potential(angular.GravitationalConstant(), rho20_scale,
                           body_radius, radius);
      harmonic_scale = std::max(harmonic_scale,
                                std::abs(expected_harmonic[radial]));
   }
   for (std::size_t radial = 0; radial < angular_integral.size(); ++radial) {
      require_analytic(angular_integral[radial][l2_m0].real(),
                       expected_harmonic[radial], harmonic_scale,
                       "l=2 Green integral differs from analytic radial moments");
   }

   const R2YlmTomography nyquist_tomography{body_radius, outer_radius,
                                            length_norm, 2};
   const Density3D nyquist_model =
       Density3D::SphericalThreeDimensionalPlanetFromModel(
           base, nyquist_tomography, harmonic_order, lmax, 0.2, outer_radius);
   const auto &nyquist_grid = nyquist_model.GSH_Grid();
   std::vector<Complex> nyquist_density_coefficients(
       static_cast<std::size_t>((lmax + 1) * (lmax + 1)), Complex{});
   const auto nyquist_density = nyquist_model.DensityAtRadialNode(1, 2);
   std::vector<Complex> nyquist_density_complex(nyquist_density.size());
   std::transform(nyquist_density.begin(), nyquist_density.end(),
                  nyquist_density_complex.begin(),
                  [](double value) { return Complex{value, 0.0}; });
   nyquist_grid.ForwardTransformation(lmax, 0, nyquist_density_complex,
                                      nyquist_density_coefficients);
   const std::size_t l2_mminus2 = 2 * 2 + 2 - 2;
   const std::size_t l2_mplus2 = 2 * 2 + 2 + 2;
   require(std::abs(nyquist_density_coefficients[l2_mminus2]) >
               1e-13 * rho0,
           "test field did not populate the negative highest-order coefficient");
   const auto nyquist_integral = GravitationalSphericalIntegral(nyquist_model);
   for (const auto &radial : nyquist_integral) {
      require(radial[l2_mminus2] == Complex{} &&
                  radial[l2_mplus2] == Complex{},
              "Green integral did not preserve both reserved Nyquist slots");
   }

   // The continuum comparison is distinct from integral algebra: use the
   // declared p=2, lmax=2 family and require volume-RMS refinement convergence.
   std::vector<double> errors;
   for (double max_step : {0.24, 0.12, 0.06}) {
      Density3D refined = Density3D::OneDimensionalPlanetFromModel(
          base, 2, lmax, max_step, outer_radius);
      const auto solved = SolveGravitationalPotential(refined, 1e-12);
      require(solved.status == Eigen::Success,
              "homogeneous sphere solver did not converge");
      errors.push_back(volume_weighted_potential_rms(
          refined, solved.solution, body_radius));
   }
   require(errors[0] > 1e-12 && errors[1] > 1e-12 && errors[2] > 1e-12,
           "mesh-convergence fixture is at the roundoff floor");
   require(errors[0] / errors[1] >= 1.5 && errors[1] / errors[2] >= 1.5,
           "volume-weighted physical potential RMS did not converge under refinement");
   std::cout << std::setprecision(8)
             << "Volume-weighted physical potential RMS (m^2/s^2), h=0.24, "
                "0.12, 0.06: "
             << errors[0] << ", " << errors[1] << ", " << errors[2]
             << " (ratios " << errors[0] / errors[1] << ", "
             << errors[1] / errors[2] << ")\n";
}
