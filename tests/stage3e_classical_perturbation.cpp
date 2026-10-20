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
constexpr double pi = std::numbers::pi_v<double>;

void require(bool ok, const std::string &message) {
   if (!ok) throw std::runtime_error(message);
}

std::size_t harmonic(int degree, int order) {
   return static_cast<std::size_t>(degree * degree + degree + order);
}

struct RadialAngularMap {
   double amplitude;
   auto RadialMapping(int) const {
      const double value = amplitude;
      return [value](double radius, double theta, double) {
         return value * radius * std::cos(theta);
      };
   }
};

struct LayerRadialMap {
   std::vector<double> scale;
   auto RadialMapping(int layer) const {
      const double value = scale.at(static_cast<std::size_t>(layer));
      return [value](double radius, double, double) { return value * radius; };
   }
};

std::vector<Complex> flatten(const GeneralEarthModels::Density3D &model,
                             const std::vector<std::vector<std::vector<Complex>>> &values) {
   const auto &mesh = model.Node_Information();
   const int order = static_cast<int>(mesh.PolynomialOrder());
   const std::size_t stride = static_cast<std::size_t>(
       (model.GSH_Grid().MaxDegree() + 1) * (model.GSH_Grid().MaxDegree() + 1));
   std::vector<Complex> flat(mesh.NumberOfGlobalRadialNodes() * stride, Complex{});
   for (std::size_t element = 0; element < mesh.NumberOfElements(); ++element) {
      for (int node = 0; node <= order; ++node) {
         const std::size_t global = mesh.GlobalNodeIndex(element, node);
         for (std::size_t lm = 0; lm < stride; ++lm)
            flat[global * stride + lm] = values[element][node][lm];
      }
   }
   return flat;
}

Eigen::VectorXcd as_eigen(const std::vector<Complex> &values) {
   Eigen::VectorXcd out(static_cast<Eigen::Index>(values.size()));
   for (std::size_t i = 0; i < values.size(); ++i)
      out[static_cast<Eigen::Index>(i)] = values[i];
   return out;
}

double physical_rms(const GeneralEarthModels::Density3D &model,
                    const Eigen::VectorXcd &field) {
   const auto &mesh = model.Node_Information();
   const auto &q = model.q();
   const std::size_t stride = static_cast<std::size_t>(
       (model.GSH_Grid().MaxDegree() + 1) * (model.GSH_Grid().MaxDegree() + 1));
   double norm2 = 0.0, volume = 0.0;
   for (std::size_t e = 0; e < mesh.NumberOfElements(); ++e) {
      for (std::size_t j = 0; j < mesh.NumNodesWithinElement(); ++j) {
         const double r = mesh.NodeRadius(e, j);
         const double weight = q.W(j) * 0.5 * mesh.ElementWidth(e) * r * r;
         const std::size_t radial = mesh.GlobalNodeIndex(e, j);
         for (std::size_t lm = 0; lm < stride; ++lm)
            norm2 += weight * std::norm(field[radial * stride + lm]);
         volume += weight;
      }
   }
   return model.PotentialNorm() * std::sqrt(norm2 / (4.0 * pi * volume));
}

Eigen::VectorXcd central(const Eigen::VectorXcd &plus,
                         const Eigen::VectorXcd &minus, double step) {
   return (plus - minus) / (2.0 * step);
}

Eigen::VectorXcd classical(const GeneralEarthModels::Density3D &model,
                           const RadialAngularMap &mapping, double tolerance) {
   const auto result = Gravity_Tools::SolveClassicalMappedPerturbation(
       const_cast<GeneralEarthModels::Density3D &>(model), mapping, tolerance);
   if (result.solves.size() != 2 ||
       result.solves[0].status != Eigen::Success ||
       result.solves[1].status != Eigen::Success ||
       !result.solution.allFinite())
      throw std::runtime_error(
          "classical perturbation did not retain both successful solve diagnostics");
   return as_eigen(flatten(model,
                           Gravity_Tools::ToHarmonicField(model, result)));
}

void check_source_and_advection() {
   constexpr int lmax = 2;
   auto model = GeneralEarthModels::Density3D::SphericalLayeredPlanet(
       {0.5, 1.0}, {2.0, 1.0}, 6371000.0, 3600.0, 5.972e24, 0.25, 1.2);
   const auto &mesh = model.Node_Information();
   RadialAngularMap source_map{0.2 * std::sqrt(3.0 / (4.0 * pi))};
   const Eigen::VectorXcd source =
       Gravity_Tools::FindBoundaryPerturbationForce(model, source_map);
   const int p = static_cast<int>(model.q().N()) - 1;
   const std::size_t stride = static_cast<std::size_t>((lmax + 1) * (lmax + 1));
   const std::size_t y10 = harmonic(1, 0);
   Eigen::VectorXcd analytic_source = Eigen::VectorXcd::Zero(source.size());
   for (std::size_t e = 0; e < mesh.NumberOfElements(); ++e) {
      for (int j = 0; j <= p; ++j) {
         const double r = mesh.NodeRadius(e, j);
         const double rho = model.Density_Point(e, j, 0);
         const double coefficient = 0.2 * r * rho;
         for (int i = 0; i <= p; ++i) {
            // The analytic Y10 coefficient is integrated separately on each
            // one-sided element trace; D_ref*w_ref contains the map factors.
            const double factor = -4.0 * pi * model.GravitationalConstant() *
                                  model.GaussDerivative(i, j) * model.q().W(j) *
                                  r * r * coefficient;
            const std::size_t global = mesh.GlobalNodeIndex(e, i);
            analytic_source[global * stride + y10] += factor;
         }
      }
   }
   const double source_error = (source - analytic_source).norm();
   require(source_error <= 1e-11 * analytic_source.norm() +
                                1e-13 * std::max(1.0, analytic_source.norm()),
           "boundary source differs from analytic Y10 element quadrature");

   std::size_t interface = mesh.NumberOfElements();
   for (std::size_t e = 1; e < mesh.NumberOfElements(); ++e) {
      if (mesh.LayerNumber(e - 1) != mesh.LayerNumber(e)) {
         interface = e;
         break;
      }
   }
   require(interface < mesh.NumberOfElements(), "layer interface missing");
   const auto left_density = model.Density_Point(interface - 1, mesh.PolynomialOrder(), 0);
   const auto right_density = model.Density_Point(interface, 0, 0);
   require(left_density != right_density, "fixture lost its one-sided density jump");
   require(source.allFinite(), "boundary source contains nonfinite entries");
   const RadialAngularMap zero_map{0.0};
   require(Gravity_Tools::FindBoundaryPerturbationForce(model, zero_map).norm() == 0.0,
           "zero displacement produced a boundary source");

   Eigen::VectorXcd phi = Eigen::VectorXcd::Zero(source.size());
   const std::size_t y20 = harmonic(2, 0);
   for (std::size_t e = 0; e < mesh.NumberOfElements(); ++e) {
      for (std::size_t j = 0; j < mesh.NumNodesWithinElement(); ++j) {
         const std::size_t global = mesh.GlobalNodeIndex(e, j);
         const double r = mesh.NodeRadius(e, j);
         phi[global * stride + y20] = r * r;
      }
   }
   std::vector<double> scales(mesh.NumberOfLayers());
   for (std::size_t layer = 0; layer < scales.size(); ++layer)
      scales[layer] = 0.1 + 0.05 * layer;
   LayerRadialMap radial_map{scales};
   const Eigen::VectorXcd advected =
       Gravity_Tools::AdvectiveBoundaryPerturbation(model, radial_map, phi);
   Eigen::VectorXcd expected = Eigen::VectorXcd::Zero(advected.size());
   for (std::size_t e = 0; e < mesh.NumberOfElements(); ++e) {
      const std::size_t layer = mesh.LayerNumber(e);
      const bool vacuum = layer + 1 == mesh.NumberOfLayers();
      const int first = e == 0 ? 1 : 0;
      const int end = e + 1 == mesh.NumberOfElements()
                          ? static_cast<int>(mesh.NumNodesWithinElement())
                          : static_cast<int>(mesh.PolynomialOrder());
      for (int j = first; j < end; ++j) {
         const double r = mesh.NodeRadius(e, j);
         const double s = vacuum
                              ? scales[layer] * mesh.PlanetRadius() *
                                    (mesh.OuterRadius() - r) /
                                    (mesh.OuterRadius() - mesh.PlanetRadius())
                              : scales[layer] * r;
         const std::size_t global = mesh.GlobalNodeIndex(e, j);
         expected[global * stride + y20] = 2.0 * r * s;
      }
   }
   const double adv_error = (advected - expected).norm();
   require(adv_error <= 1e-11 * expected.norm() +
                              1e-13 * std::max(1.0, expected.norm()),
           "manufactured radial advection differs from its analytic value");
   require(advected.allFinite() && advected.head(stride).norm() == 0.0,
           "advection origin or finite-value contract failed");
   const LayerRadialMap zero_radial_map{std::vector<double>(scales.size(), 0.0)};
   require(Gravity_Tools::AdvectiveBoundaryPerturbation(model, zero_radial_map, phi)
                   .norm() == 0.0,
           "zero displacement produced advection");
}

struct FullMappedTangent {
   Eigen::VectorXcd full;
   Eigen::VectorXcd approximate;
   double full_rms;
   double approximate_rms;
   double error_rms;
   double half_step_error_rms;
   double step_error_ratio;
   double noise_floor;
   double half_step_noise_floor;
};

FullMappedTangent mapped_tangent(double max_step, double step) {
   constexpr int p = 5;
   constexpr int lmax = 2;
   constexpr double outer = 1.5;
   constexpr double solver_tolerance = 1e-12;
   const double length = 6371000.0;
   const double time = 3600.0;
   const double mass = 5.972e24;
   const auto material = GeneralEarthModels::SimpleModels::spherical_model::HomogeneousSphere(
       1.0, 1.0, length, time, mass);
   const RadialAngularMap zero{0.0};
   auto base = GeneralEarthModels::Density3D(
       material, PlanetaryModel::TomographyZeroModel(), zero, p, lmax,
       max_step, outer);
   const auto base_source = Gravity_Tools::FindForce(base);
   const auto base_result =
       Gravity_Tools::SolvePotential(base, solver_tolerance);
   require(base_result.final_solve().status == Eigen::Success,
           "base solve failed");
   const double base_rms = physical_rms(base, base_result.solution);
   Eigen::VectorXcd full_estimates[2], classical_estimates[2];
   double error_rms[2] = {};
   double noise_floors[2] = {};
   const double amplitudes[2] = {step, step / 2.0};
   for (int k = 0; k < 2; ++k) {
      const double epsilon = amplitudes[k];
      const RadialAngularMap plus_map{epsilon}, minus_map{-epsilon};
      auto plus = GeneralEarthModels::Density3D(
          material, PlanetaryModel::TomographyZeroModel(), plus_map, p, lmax,
          max_step, outer);
      auto minus = GeneralEarthModels::Density3D(
          material, PlanetaryModel::TomographyZeroModel(), minus_map, p, lmax,
          max_step, outer);
      const auto plus_source = Gravity_Tools::FindForce(plus);
      const auto minus_source = Gravity_Tools::FindForce(minus);
      require((plus_source - base_source).norm() <=
                  1e-11 * base_source.norm() + 1e-13 * std::max(1.0, base_source.norm()) &&
                  (minus_source - base_source).norm() <=
                  1e-11 * base_source.norm() + 1e-13 * std::max(1.0, base_source.norm()),
              "mapped referential source is not invariant");
      const auto plus_result =
          Gravity_Tools::SolvePotential(plus, solver_tolerance);
      const auto minus_result =
          Gravity_Tools::SolvePotential(minus, solver_tolerance);
      require(plus_result.final_solve().status == Eigen::Success &&
                  minus_result.final_solve().status == Eigen::Success,
              "full mapped solve failed");
      full_estimates[k] = central(plus_result.solution, minus_result.solution, epsilon);
      const auto boundary_source =
          Gravity_Tools::FindBoundaryPerturbationForce(base, plus_map);
      require(boundary_source.norm() >
                  1e-8 * std::max(1.0, base_source.norm()),
              "boundary-source signal is too small in source-vector units");
      const auto plus_classical = classical(base, plus_map, solver_tolerance);
      const auto minus_classical = classical(base, minus_map, solver_tolerance);
      classical_estimates[k] = central(plus_classical, minus_classical, epsilon);
      error_rms[k] = physical_rms(base, full_estimates[k] - classical_estimates[k]);
      noise_floors[k] =
          std::max(1e-10 * std::max(1.0, base_rms),
                   100.0 * solver_tolerance * base_rms / epsilon);
      const double classical_rms = physical_rms(base, classical_estimates[k]);
      require(classical_rms > 1e-8 * std::max(1.0, base_rms),
              "classical tangent signal is too small");
      require(error_rms[k] <= 2e-4 * classical_rms + noise_floors[k],
              "individual classical tangent error exceeds the predeclared bound");
   }
   const double full_rms = physical_rms(base, full_estimates[0]);
   const double classical_rms = physical_rms(base, classical_estimates[0]);
   require(full_rms > 1e-8 * std::max(1.0, base_rms),
           "full mapped tangent signal is too small");
   const double estimate_difference = physical_rms(base, full_estimates[0] - full_estimates[1]);
   require(estimate_difference <= 1e-4 * classical_rms + noise_floors[0] +
                                       noise_floors[1],
           "successive full-mapped tangent estimates disagree");
   return {full_estimates[0], classical_estimates[0], full_rms,
           classical_rms, error_rms[0], error_rms[1],
           error_rms[0] / std::max(error_rms[1], std::numeric_limits<double>::min()),
           noise_floors[0], noise_floors[1]};
}

}  // namespace

int main() {
   check_source_and_advection();
   const auto coarse = mapped_tangent(0.1, 1e-3);
   const auto fine_mesh = mapped_tangent(0.05, 1e-3);
   std::cout << "classical mapped tangent RMS: " << coarse.approximate_rms
             << " full mapped: " << coarse.full_rms
             << " error: " << coarse.error_rms
             << " half-step error: " << coarse.half_step_error_rms
             << " half-step noise floor: " << coarse.half_step_noise_floor
             << " step-error ratio: " << coarse.step_error_ratio
             << " refined-mesh error: " << fine_mesh.error_rms
             << " mesh-error ratio: " << coarse.error_rms / fine_mesh.error_rms
             << " coarse noise floor: " << coarse.noise_floor
             << " refined-mesh noise floor: " << fine_mesh.noise_floor << '\n';
   if (coarse.error_rms > coarse.noise_floor &&
       fine_mesh.error_rms > fine_mesh.noise_floor) {
      require(fine_mesh.error_rms < coarse.error_rms,
              "resolved tangent discrepancy did not decrease on mesh refinement");
   }
   return 0;
}
