#include <gplspec/All>

int consumer_has_model();

int main() {
  auto model = GeneralEarthModels::Density3D::SphericalHomogeneousPlanet(
      1.0, 1.0, 6371000.0, 3600.0, 5.972e24, 0.5, 1.2);
  const bool main_translation_unit_uses_gplspec =
      model.Num_Elements() > 0 && !model.Density().empty();
  return main_translation_unit_uses_gplspec && consumer_has_model() ? 0 : 1;
}
