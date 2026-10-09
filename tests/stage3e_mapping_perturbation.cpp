#include <gplspec/All>
#include <gplspec/src/Pseudospectral_Matrix_Wrapper3D.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Complex = std::complex<double>;
using Operator = MatrixReplacement3D<Complex>;

struct SmoothRadialMap {
   double amplitude;
   double body_radius;

   auto RadialMapping(int) const {
      const double amplitude_value = amplitude;
      const double body_radius_value = body_radius;
      return [amplitude_value, body_radius_value](double radius,
                                                  double theta, double) {
         const double normalized = radius / body_radius_value;
         const double angular =
             std::cos(theta) + 0.25 * (3.0 * std::cos(theta) *
                                           std::cos(theta) -
                                       1.0);
         return amplitude_value * radius * (1.0 - normalized) * angular;
      };
   }
};

struct ZeroRadialMap {
   auto RadialMapping(int) const {
      return [](double, double, double) { return 0.0; };
   }
};

void require(bool condition, const std::string &message) {
   if (!condition) throw std::runtime_error(message);
}

std::size_t scalar_index(int degree, int order) {
   return static_cast<std::size_t>(degree * degree + degree + order);
}

double norm_limit(double relative, double absolute_scale, double scale) {
   return relative * scale + absolute_scale * std::max(1.0, scale);
}

void require_vector_close(const Eigen::VectorXcd &actual,
                          const Eigen::VectorXcd &expected,
                          double relative, double absolute_scale,
                          const std::string &label) {
   require(actual.size() == expected.size(), label + ": vector size mismatch");
   const double error = (actual - expected).norm();
   const double scale = expected.norm();
   if (!actual.allFinite() || !expected.allFinite() ||
       error > norm_limit(relative, absolute_scale, scale)) {
      throw std::runtime_error(label + ": error=" + std::to_string(error) +
                               " scale=" + std::to_string(scale));
   }
}

Eigen::VectorXcd surface_field(const GeneralEarthModels::Density3D &model,
                               const Eigen::VectorXcd &solution) {
   const auto &mesh = model.Node_InformationP();
   const auto &grid = model.GSH_Grid();
   const std::size_t stride =
       static_cast<std::size_t>(grid.MaxDegree() + 1) *
       (grid.MaxDegree() + 1);
   const double radius = mesh.PlanetRadius();
   std::size_t surface_radial = mesh.NumberOfGlobalRadialNodes();
   for (std::size_t element = 0; element < mesh.NumberOfElements(); ++element) {
      for (std::size_t node = 0; node < mesh.NumNodesWithinElement(); ++node) {
         if (std::abs(mesh.NodeRadius(element, node) - radius) <=
             1e-13 * std::max(1.0, radius)) {
            surface_radial = mesh.GlobalNodeIndex(element, node);
         }
      }
   }
   require(surface_radial < mesh.NumberOfGlobalRadialNodes(),
           "body-surface radial node was not found");

   std::vector<Complex> coefficients(stride),
       spatial(grid.NumberOfCoLatitudes() * grid.NumberOfLongitudes(),
               Complex{});
   for (std::size_t harmonic = 0; harmonic < stride; ++harmonic) {
      coefficients[harmonic] = solution[surface_radial * stride + harmonic];
   }
   grid.InverseTransformation(grid.MaxDegree(), 0, coefficients, spatial);

   Eigen::VectorXcd values(static_cast<Eigen::Index>(spatial.size()));
   for (std::size_t point = 0; point < spatial.size(); ++point) {
      values[static_cast<Eigen::Index>(point)] = spatial[point];
   }
   return values;
}

double angular_rms(const GeneralEarthModels::Density3D &model,
                   const Eigen::VectorXcd &field) {
   const auto &grid = model.GSH_Grid();
   double weighted_norm2 = 0.0;
   double weight_sum = 0.0;
   std::size_t point = 0;
   for (std::size_t theta = 0; theta < grid.NumberOfCoLatitudes(); ++theta) {
      for (std::size_t phi = 0; phi < grid.NumberOfLongitudes(); ++phi) {
         const double weight = grid.CoLatitudeWeights()[theta] *
                               grid.LongitudeWeights()[phi];
         weighted_norm2 += weight * std::norm(field[point++]);
         weight_sum += weight;
      }
   }
   return model.PotentialNorm() * std::sqrt(weighted_norm2 / weight_sum);
}

Eigen::VectorXcd central_difference(const Eigen::VectorXcd &plus,
                                    const Eigen::VectorXcd &minus,
                                    double amplitude) {
   return (plus - minus) / (2.0 * amplitude);
}

}  // namespace

int main() {
   constexpr int lmax = 2;
   constexpr double max_step = 0.5;
   constexpr double outer_radius = 1.2;
   constexpr double solver_tolerance = 1e-13;
   const auto quadrature =
       GaussQuad::GaussLobattoLegendreQuadrature1D<double>(3);
   GeneralEarthModels::RadialEarthModel<> material(
       GPLSpec::ModelInput::ReadPremDeck(
           std::string(GPLSPEC_SOURCE_DIR) + "/modeldata/PREM500"));
   GeneralEarthModels::Density3D base =
       GeneralEarthModels::Density3D::OneDimensionalPlanetFromModel(
           material, 2, lmax, max_step, outer_radius);
   const double body_radius = base.Node_InformationP().PlanetRadius();
   const SmoothRadialMap unit_map{1.0, body_radius};
   const double pi = std::numbers::pi_v<double>;

   const Eigen::VectorXcd source = Gravity_Tools::FindForce(base);
   Operator base_operator(base, false);
   const std::size_t dimension = static_cast<std::size_t>(base_operator.cols());
   Eigen::VectorXcd input(static_cast<Eigen::Index>(dimension));
   for (Eigen::Index i = 0; i < input.size(); ++i) {
      input[i] = Complex(0.125 * (i + 1), -0.0625 * (i + 2));
   }
   const Eigen::VectorXcd base_action = base_operator.Apply(input);

   GeneralEarthModels::MappingPerturbation perturbation(base, unit_map);
   Operator delta_operator(base, perturbation, false);
   const Eigen::VectorXcd delta_action = delta_operator.Apply(input);
   require(delta_action.norm() >
               1e-8 * std::max(1.0, base_action.norm()),
           "mapping perturbation operator signal is too small");

   for (double amplitude : {1e-3, 5e-4}) {
      const SmoothRadialMap plus_map{amplitude, body_radius};
      const SmoothRadialMap minus_map{-amplitude, body_radius};
      GeneralEarthModels::Density3D plus_model(
          material, PlanetaryModel::TomographyZeroModel(), plus_map, 2, lmax,
          max_step, outer_radius);
      GeneralEarthModels::Density3D minus_model(
          material, PlanetaryModel::TomographyZeroModel(), minus_map, 2, lmax,
          max_step, outer_radius);
      Operator plus_operator(plus_model, false);
      Operator minus_operator(minus_model, false);

      require_vector_close(Gravity_Tools::FindForce(plus_model), source,
                           1e-11, 1e-13,
                           "plus mapped model source invariance");
      require_vector_close(Gravity_Tools::FindForce(minus_model), source,
                           1e-11, 1e-13,
                           "minus mapped model source invariance");

      const Eigen::VectorXcd action_fd = central_difference(
          plus_operator.Apply(input), minus_operator.Apply(input), amplitude);
      const double action_error = (action_fd - delta_action).norm();
      const double action_bound =
          2e-5 * delta_action.norm() +
          1e-11 * std::max(1.0, base_action.norm());
      require(action_error <= action_bound,
              "mapped operator tensor finite difference: error=" +
                  std::to_string(action_error) + " bound=" +
                  std::to_string(action_bound));
   }

   GeneralEarthModels::MappingPerturbation zero_perturbation(base, ZeroRadialMap{});
   const auto &zero_tensor = zero_perturbation.delta_laplace_tensor();
   for (const auto &element : zero_tensor)
      for (const auto &node : element)
         for (const auto &tensor : node)
            require(tensor.norm() == 0.0,
                    "zero mapping produced a nonzero tensor");
   Operator zero_operator(base, zero_perturbation, false);
   require(zero_operator.Apply(input).norm() == 0.0,
           "zero mapping produced a nonzero action");

   using CG = Eigen::ConjugateGradient<
       Operator, Eigen::Lower | Eigen::Upper,
       Eigen::SphericalGeometryPreconditioner<Complex>>;
   const auto preconditioner =
       Gravity_Tools::detail::SphericalPreconditionerMatrix(base, lmax);
   Operator base_solve_operator(base, true);
   Gravity_Tools::detail::PreparedSolver<CG, Operator> base_solver(
       base_solve_operator, preconditioner);
   const auto base_result =
       base_solver.Solve(source, solver_tolerance, true, true);
   require(base_result.status == Eigen::Success,
           "base mapped solve did not converge: status=" +
               std::to_string(base_result.status) + " iterations=" +
               std::to_string(base_result.iterations) + " recursive_error=" +
               std::to_string(base_result.recursive_error) +
               " relative_residual=" +
               std::to_string(base_result.relative_residual));

   Eigen::VectorXcd delta_rhs = -delta_operator.Apply(base_result.solution);
   const auto tangent_result =
       base_solver.Solve(delta_rhs, solver_tolerance, true, true);
   require(tangent_result.status == Eigen::Success,
           "mapping tangent solve did not converge");
   const Eigen::VectorXcd base_surface = surface_field(base, base_result.solution);
   const Eigen::VectorXcd tangent_surface =
       surface_field(base, tangent_result.solution);
   const double base_surface_rms = angular_rms(base, base_surface);
   const double tangent_surface_rms = angular_rms(base, tangent_surface);
   require(tangent_surface_rms > 1e-8 * base_surface_rms,
           "mapping perturbation field signal is too small");

   const Eigen::Index surface_stride = (lmax + 1) * (lmax + 1);
   std::size_t surface_radial = base.Node_InformationP().NumberOfGlobalRadialNodes();
   for (std::size_t element = 0;
        element < base.Node_InformationP().NumberOfElements(); ++element) {
      for (std::size_t node = 0;
           node < base.Node_InformationP().NumNodesWithinElement(); ++node) {
         if (std::abs(base.Node_InformationP().NodeRadius(element, node) -
                      body_radius) <= 1e-13 * std::max(1.0, body_radius)) {
            surface_radial = base.Node_InformationP().GlobalNodeIndex(element, node);
         }
      }
   }
   for (int degree : {1, 2}) {
      double harmonic_norm2 = 0.0;
      for (int order = -degree; order <= degree; ++order) {
         const auto index = static_cast<Eigen::Index>(surface_radial) *
                                surface_stride +
                            static_cast<Eigen::Index>(scalar_index(degree, order));
         harmonic_norm2 += std::norm(tangent_result.solution[index]);
      }
      const double physical_degree_rms =
          base.PotentialNorm() * std::sqrt(harmonic_norm2 / (4.0 * pi));
      require(physical_degree_rms > 1e-8 * base_surface_rms,
              "mapping tangent lost a declared l=1 or l=2 signal");
   }

   std::vector<double> amplitudes{1e-3, 5e-4};
   std::vector<double> field_errors;
   std::vector<Eigen::VectorXcd> field_derivatives;
   for (double amplitude : amplitudes) {
      const SmoothRadialMap plus_map{amplitude, body_radius};
      const SmoothRadialMap minus_map{-amplitude, body_radius};
      GeneralEarthModels::Density3D plus_model(
          material, PlanetaryModel::TomographyZeroModel(), plus_map, 2, lmax,
          max_step, outer_radius);
      GeneralEarthModels::Density3D minus_model(
          material, PlanetaryModel::TomographyZeroModel(), minus_map, 2, lmax,
          max_step, outer_radius);
      Operator plus_operator(plus_model, true);
      Operator minus_operator(minus_model, true);
      Gravity_Tools::detail::PreparedSolver<CG, Operator> plus_solver(
          plus_operator, preconditioner);
      Gravity_Tools::detail::PreparedSolver<CG, Operator> minus_solver(
          minus_operator, preconditioner);
      const auto plus_result =
          plus_solver.Solve(source, solver_tolerance, true, true);
      const auto minus_result =
          minus_solver.Solve(source, solver_tolerance, true, true);
      require(plus_result.status == Eigen::Success &&
                  minus_result.status == Eigen::Success,
              "finite-difference mapped solve did not converge");
      const auto field_fd = central_difference(
          surface_field(plus_model, plus_result.solution),
          surface_field(minus_model, minus_result.solution), amplitude);
      const auto field_error = field_fd - tangent_surface;
      const double error_rms = angular_rms(base, field_error);
      const double bound = 2e-4 * tangent_surface_rms +
                           1e-10 * std::max(1.0, base_surface_rms);
      require(error_rms <= bound,
              "physical surface potential tangent finite difference failed");
      require(amplitude * tangent_surface_rms >
                  100.0 * solver_tolerance * base_surface_rms,
              "field finite-difference signal is below solver noise scale");
      field_errors.push_back(error_rms);
      field_derivatives.push_back(field_fd);
   }

   const auto derivative_difference = field_derivatives[0] - field_derivatives[1];
   const double derivative_difference_rms =
       angular_rms(base, derivative_difference);
   const auto solve_noise_floor = [&](double amplitude) {
      return std::max(
          1e-10 * std::max(1.0, base_surface_rms),
          100.0 * solver_tolerance * base_surface_rms / amplitude);
   };
   const double coarse_noise_floor = solve_noise_floor(amplitudes[0]);
   const double fine_noise_floor = solve_noise_floor(amplitudes[1]);
   const double agreement_limit = 1e-4 * tangent_surface_rms +
                                  coarse_noise_floor + fine_noise_floor;
   require(derivative_difference_rms <= agreement_limit,
           "successive physical-field derivative estimates disagree: difference=" +
               std::to_string(derivative_difference_rms) + " bound=" +
               std::to_string(agreement_limit));
   const double refinement_ratio = field_errors[0] / field_errors[1];
   const bool refinement_above_noise =
       field_errors[0] > coarse_noise_floor &&
       field_errors[1] > fine_noise_floor;
   if (refinement_above_noise) {
      require(refinement_ratio >= 3.0 && refinement_ratio <= 5.0,
              "central-difference error does not show the declared second-order refinement: ratio=" +
                  std::to_string(refinement_ratio));
   }

   std::cout << "Stage 3E.2 mapping tangent: operator central differences at "
                "1e-3/5e-4 and physical surface RMS errors "
             << field_errors[0] << ", " << field_errors[1]
             << " m^2/s^2 per unit amplitude; successive derivative difference "
             << derivative_difference_rms << " (limit " << agreement_limit
             << "); second-order error ratio " << refinement_ratio
             << (refinement_above_noise ? " (resolved)" : " (noise-limited)")
             << "; successive difference/base RMS "
             << derivative_difference_rms / base_surface_rms
             << "; tangent/base RMS "
             << tangent_surface_rms / base_surface_rms
             << "; source and zero-map checks passed\n";
}
