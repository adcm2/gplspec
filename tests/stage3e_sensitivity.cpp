#include <gplspec/All>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Complex = std::complex<double>;
using SensitivityField =
    std::vector<std::vector<std::vector<Complex>>>;

void require(bool condition, const std::string &message) {
   if (!condition) throw std::runtime_error(message);
}

struct DensityY20Tomography {
   double epsilon;
   double body_radius;
   double length_norm;

   double GetValueAt(double depth, double, double latitude,
                     bool = false) const {
      const double radius = body_radius - depth / (length_norm / 1000.0);
      const double sin_latitude = std::sin(latitude * std::numbers::pi_v<double> / 180.0);
      const double y20 = std::sqrt(5.0 / (16.0 * std::numbers::pi_v<double>)) *
                         (3.0 * sin_latitude * sin_latitude - 1.0);
      const double radial = radius / body_radius;
      // Density3D applies (1 + 0.005*GetValueAt(...)).
      return (epsilon / 0.005) * radial * radial * y20;
   }
   std::vector<double> GetDepths() const { return {-1.0, 10000.0}; }
   std::vector<double> GetLatitudes() const { return {-90.0, 90.0}; }
   std::vector<double> GetLongitudes() const { return {0.0, 360.0}; }
   std::vector<double> GetValues() const { return {0.0}; }
};

std::size_t harmonic_index(int degree, int order) {
   return static_cast<std::size_t>(degree * degree + degree + order);
}

template <class Model>
Complex outer_flat_coefficient(const Model &model,
                               const Eigen::VectorXcd &solution,
                               int degree, int order) {
   const auto &mesh = model.Node_InformationP();
   const std::size_t stride = static_cast<std::size_t>(
       model.GSH_Grid().MaxDegree() + 1) *
                              (model.GSH_Grid().MaxDegree() + 1);
   const std::size_t index =
       (mesh.NumberOfGlobalRadialNodes() - 1) * stride +
       harmonic_index(degree, order);
   return solution[static_cast<Eigen::Index>(index)];
}

Complex outer_nested_coefficient(const SensitivityField &field,
                                 int degree, int order) {
   return field.back().back()[harmonic_index(degree, order)];
}

double field_norm(const SensitivityField &field) {
   double sum = 0.0;
   for (const auto &element : field)
      for (const auto &node : element)
         for (const auto &value : node) sum += std::norm(value);
   return std::sqrt(sum);
}

void require_near(double actual, double expected, double relative,
                  const std::string &label) {
   const double tolerance = relative * std::abs(expected);
   if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
      throw std::runtime_error(label + ": actual=" + std::to_string(actual) +
                               " expected=" + std::to_string(expected) +
                               " tolerance=" + std::to_string(tolerance));
   }
}

void check_selector_rhs(GeneralEarthModels::Density3D &model) {
   constexpr double legacy_pi = 3.1415926535;
   const int lmax = model.GSH_Grid().MaxDegree();
   const std::size_t stride = static_cast<std::size_t>(lmax + 1) * (lmax + 1);
   const std::size_t outer =
       (model.Node_InformationP().NumberOfGlobalRadialNodes() - 1) * stride;
   const double base = 4.0 * legacy_pi * model.GravitationalConstant();

   const Eigen::VectorXcd rhs = Gravity_Tools::detail::SensitivitySelectorRhs(
       model, {2, 2, 2}, {1, -1, 0}, {3.0, 4.0, 2.0});
   require(rhs.size() == static_cast<Eigen::Index>(
                             model.Node_InformationP().NumberOfGlobalRadialNodes() *
                             stride),
           "selector RHS has the wrong dimension");
   require(rhs[static_cast<Eigen::Index>(outer + harmonic_index(2, -1))] ==
               Complex{-3.0 * base, 0.0},
           "m=+1 index or phase differs");
   require(rhs[static_cast<Eigen::Index>(outer + harmonic_index(2, 1))] ==
               Complex{-4.0 * base, 0.0},
           "m=-1 index or phase differs");
   require(rhs[static_cast<Eigen::Index>(outer + harmonic_index(2, 0))] ==
               Complex{2.0 * base, 0.0},
           "m=0 selector differs");

   const Eigen::VectorXcd duplicates =
       Gravity_Tools::detail::SensitivitySelectorRhs(
           model, {2, 2}, {0, 0}, {1.0, 2.0});
   require(duplicates[static_cast<Eigen::Index>(outer + harmonic_index(2, 0))] ==
               Complex{2.0 * base, 0.0},
           "duplicate selectors no longer use last-value-wins");
   const Eigen::VectorXcd empty =
       Gravity_Tools::detail::SensitivitySelectorRhs(model, {}, {}, {});
   require(empty.isZero(0.0), "empty selector list must produce zero RHS");

   auto must_reject = [&](const std::vector<int> &degrees,
                          const std::vector<int> &orders,
                          const std::vector<double> &weights,
                          const std::string &label) {
      bool rejected = false;
      try {
         (void)Gravity_Tools::detail::SensitivitySelectorRhs(
             model, degrees, orders, weights);
      } catch (const std::invalid_argument &) {
         rejected = true;
      } catch (const std::out_of_range &) {
         rejected = true;
      }
      require(rejected, label);
   };
   must_reject({2}, {}, {1.0}, "selector length mismatch was accepted");
   must_reject({-1}, {0}, {1.0}, "negative degree was accepted");
   must_reject({3}, {0}, {1.0}, "degree above lmax was accepted");
   must_reject({2}, {3}, {1.0}, "order outside degree was accepted");
   must_reject({2}, {0}, {std::numeric_limits<double>::infinity()},
               "non-finite selector weight was accepted");
   bool rejected_tolerance = false;
   try {
      (void)Gravity_Tools::SphericalHarmonicSensitivityKernel(
          model, 2, 0, 0.0);
   } catch (const std::invalid_argument &) {
      rejected_tolerance = true;
   }
   require(rejected_tolerance, "non-positive sensitivity tolerance was accepted");
   rejected_tolerance = false;
   try {
      (void)Gravity_Tools::SphericalHarmonicSensitivityKernel(
          model, std::vector<int>{2}, std::vector<int>{0},
          std::vector<double>{1.0}, std::numeric_limits<double>::quiet_NaN());
   } catch (const std::invalid_argument &) {
      rejected_tolerance = true;
   }
   require(rejected_tolerance, "non-finite sensitivity tolerance was accepted");
}

}  // namespace

int main() {
   constexpr int order = 5;
   constexpr int lmax = 2;
   constexpr double body_radius = 0.63;
   constexpr double outer_radius = 1.2;
   constexpr double max_element_step = 0.1;
   constexpr double solve_tolerance = 1e-13;
   constexpr double amplitude_coarse = 1e-4;
   constexpr double amplitude_fine = 5e-5;
   constexpr double analytic_relative_tolerance = 1e-4;
   constexpr double kernel_relative_tolerance = 1e-11;
   constexpr double vector_route_relative_tolerance = 1e-8;
   constexpr double legacy_pi = 3.1415926535;
   constexpr double full_pi = std::numbers::pi_v<double>;
   constexpr double G_physical = 6.6743015e-11;
   constexpr double length_norm = 6371000.0;
   constexpr double time_norm = 3600.0;
   constexpr double mass_norm = 5.972e24;

   using SimpleModel = GeneralEarthModels::SimpleModels::spherical_model;
   const SimpleModel material = SimpleModel::HomogeneousSphere(
       body_radius, 2.0, length_norm, time_norm, mass_norm);
   auto base = GeneralEarthModels::Density3D::OneDimensionalPlanetFromModel(
       material, order, lmax, max_element_step, outer_radius);
   check_selector_rhs(base);

   auto zero_model =
       GeneralEarthModels::Density3D::SphericalThreeDimensionalPlanetFromModel(
           material, DensityY20Tomography{0.0, body_radius, length_norm}, order,
           lmax, max_element_step, outer_radius);
   const Eigen::VectorXcd base_source = Gravity_Tools::FindForce(base);
   const Eigen::VectorXcd zero_source = Gravity_Tools::FindForce(zero_model);
   require((zero_source - base_source).norm() == 0.0,
           "zero-amplitude tomography changed the source");

   auto plus_coarse =
       GeneralEarthModels::Density3D::SphericalThreeDimensionalPlanetFromModel(
           material,
           DensityY20Tomography{amplitude_coarse, body_radius, length_norm},
           order, lmax, max_element_step, outer_radius);
   auto minus_coarse =
       GeneralEarthModels::Density3D::SphericalThreeDimensionalPlanetFromModel(
           material,
           DensityY20Tomography{-amplitude_coarse, body_radius, length_norm},
           order, lmax, max_element_step, outer_radius);
   auto plus_fine =
       GeneralEarthModels::Density3D::SphericalThreeDimensionalPlanetFromModel(
           material, DensityY20Tomography{amplitude_fine, body_radius, length_norm},
           order, lmax, max_element_step, outer_radius);
   auto minus_fine =
       GeneralEarthModels::Density3D::SphericalThreeDimensionalPlanetFromModel(
           material, DensityY20Tomography{-amplitude_fine, body_radius, length_norm},
           order, lmax, max_element_step, outer_radius);

   const auto base_result =
       Gravity_Tools::SolveGravitationalPotential(base, solve_tolerance);
   const auto plus_coarse_result =
       Gravity_Tools::SolveGravitationalPotential(plus_coarse, solve_tolerance);
   const auto minus_coarse_result =
       Gravity_Tools::SolveGravitationalPotential(minus_coarse, solve_tolerance);
   const auto plus_fine_result =
       Gravity_Tools::SolveGravitationalPotential(plus_fine, solve_tolerance);
   const auto minus_fine_result =
       Gravity_Tools::SolveGravitationalPotential(minus_fine, solve_tolerance);
   for (const auto *result : {&base_result, &plus_coarse_result,
                              &minus_coarse_result, &plus_fine_result,
                              &minus_fine_result}) {
      require(result->status == Eigen::Success,
              "full density solve did not converge");
      require(result->solution.allFinite(), "full density solve is non-finite");
   }

   const Complex base_outer = outer_flat_coefficient(base, base_result.solution, 2, 0);
   const Complex plus_outer_coarse =
       outer_flat_coefficient(plus_coarse, plus_coarse_result.solution, 2, 0);
   const Complex minus_outer_coarse =
       outer_flat_coefficient(minus_coarse, minus_coarse_result.solution, 2, 0);
   const Complex plus_outer_fine =
       outer_flat_coefficient(plus_fine, plus_fine_result.solution, 2, 0);
   const Complex minus_outer_fine =
       outer_flat_coefficient(minus_fine, minus_fine_result.solution, 2, 0);
   const double potential_norm = base.PotentialNorm();
   const double derivative_coarse =
       (plus_outer_coarse.real() - minus_outer_coarse.real()) /
       (2.0 * amplitude_coarse) * potential_norm;
   const double derivative_fine =
       (plus_outer_fine.real() - minus_outer_fine.real()) /
       (2.0 * amplitude_fine) * potential_norm;

   const double density_norm = base.DensityNorm();
   const double rho_physical = base.Density_Point(0, 0, 0) * density_norm;
   const double body_physical = body_radius * length_norm;
   const double outer_physical = outer_radius * length_norm;
   const double analytic_physical =
       -4.0 * full_pi * G_physical * rho_physical *
       std::pow(body_physical, 5) /
       (35.0 * std::pow(outer_physical, 3));
   require(std::abs(analytic_physical) > 0.0,
           "analytic density response unexpectedly vanishes");
   require_near(derivative_coarse, analytic_physical,
                analytic_relative_tolerance,
                "coarse full-model density derivative differs from analytic value");
   require_near(derivative_fine, analytic_physical,
                analytic_relative_tolerance,
                "fine full-model density derivative differs from analytic value");
   require_near(derivative_coarse, derivative_fine,
                analytic_relative_tolerance,
                "successive full-model density derivatives disagree");

   const double symmetry_limit =
       1e-8 * std::abs(base_outer.real() * potential_norm) +
       1e-12 * G_physical * base.Mass() / outer_physical;
   const double coarse_mean_error =
       std::abs(0.5 * (plus_outer_coarse.real() + minus_outer_coarse.real()) -
                base_outer.real()) * potential_norm;
   const double fine_mean_error =
       std::abs(0.5 * (plus_outer_fine.real() + minus_outer_fine.real()) -
                base_outer.real()) * potential_norm;
   require(coarse_mean_error <= symmetry_limit && fine_mean_error <= symmetry_limit,
           "plus/minus density solutions violate central symmetry");

   const auto single_kernel = Gravity_Tools::SphericalHarmonicSensitivityKernel(
       base, 2, 0, solve_tolerance);
   const auto weighted_kernel =
       Gravity_Tools::SphericalHarmonicSensitivityKernel(
           base, {2, 2}, {0, 0},
           {1.0, 2.0}, solve_tolerance);
   const double single_norm = field_norm(single_kernel);
   double weighted_difference2 = 0.0;
   for (std::size_t e = 0; e < single_kernel.size(); ++e)
      for (std::size_t j = 0; j < single_kernel[e].size(); ++j)
         for (std::size_t k = 0; k < single_kernel[e][j].size(); ++k)
            weighted_difference2 +=
                std::norm(weighted_kernel[e][j][k] - 2.0 * single_kernel[e][j][k]);
   require(std::sqrt(weighted_difference2) <=
               vector_route_relative_tolerance * 2.0 * single_norm,
           "weighted vector sensitivity differs from its duplicate last value");

   const auto &mesh = base.Node_InformationP();
   const auto &quadrature = base.qP();
   const auto &kernel = single_kernel;
   Complex kernel_contraction{};
   for (std::size_t element = 0; element < mesh.NumberOfElements(); ++element) {
      for (std::size_t node = 0; node < mesh.NumNodesWithinElement(); ++node) {
         const double radius = mesh.NodeRadius(element, node);
         const double base_density = base.Density_Point(element, node, 0);
         if (base_density == 0.0) continue;
         const double radial_weight =
             0.5 * mesh.ElementWidth(element) * quadrature.W(node) *
             radius * radius;
         const double delta_rho_lm =
             base_density * std::pow(radius / body_radius, 2);
         // The solved kernel already contains 4*pi_legacy*G in its boundary
         // RHS. The primal 4*pi*G source factor cancels against that adjoint
         // normalization, leaving this weighted density contraction.
         kernel_contraction +=
             kernel[element][node][harmonic_index(2, 0)] * delta_rho_lm *
             radial_weight;
      }
   }
   const double kernel_physical = kernel_contraction.real() * potential_norm;
   const double legacy_adjusted_analytic = analytic_physical * legacy_pi / full_pi;
   require_near(kernel_physical, legacy_adjusted_analytic,
                kernel_relative_tolerance,
                "independent kernel volume contraction differs from analytic value");

   std::cout << "analytic physical derivative=" << analytic_physical << " m^2/s^2\n"
             << "kernel quadrature=" << kernel_physical << " m^2/s^2\n"
             << "finite differences coarse/fine=" << derivative_coarse << "/"
             << derivative_fine << " m^2/s^2\n"
             << "FD relative errors="
             << std::abs(derivative_coarse / analytic_physical - 1.0) << "/"
             << std::abs(derivative_fine / analytic_physical - 1.0) << "\n";
}
