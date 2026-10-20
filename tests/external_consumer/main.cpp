#include <gplspec/All>

int consumer_has_model();

int main() {
  auto model = GeneralEarthModels::Density3D::SphericalHomogeneousPlanet(
      1.0, 1.0, 6371000.0, 3600.0, 5.972e24, 0.5, 1.2);
  const auto result = Gravity_Tools::SolvePotential(model, 1e-6);
  const auto field = Gravity_Tools::ToHarmonicField(model, result);
  const bool main_translation_unit_uses_gplspec =
      model.Num_Elements() > 0 && !model.Density().empty() &&
      result.solves.size() == 1 &&
      result.final_solve().status == Eigen::Success &&
      result.solution.allFinite() && !field.empty();
  return main_translation_unit_uses_gplspec && consumer_has_model() ? 0 : 1;
}
