#include <gplspec/All>

int consumer_has_model() {
  auto model = GeneralEarthModels::Density3D::SphericalHomogeneousPlanet(
      1.0, 1.0, 6371000.0, 3600.0, 5.972e24, 0.5, 1.2);
  const auto result = Gravity_Tools::SolvePotential(model, 1e-6);
  const auto field = Gravity_Tools::ToHarmonicField(model, result);
  return model.Num_Elements() > 0 && model.Density().size() > 0 &&
         result.solves.size() == 1 &&
         result.final_solve().status == Eigen::Success &&
         result.solution.allFinite() && !field.empty();
}
