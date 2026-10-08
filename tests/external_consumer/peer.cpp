#include <gplspec/All>

int consumer_has_model() {
  auto model = GeneralEarthModels::Density3D::SphericalHomogeneousPlanet(
      1.0, 1.0, 6371000.0, 3600.0, 5.972e24, 0.5, 1.2);
  return model.Num_Elements() > 0 && model.Density().size() > 0;
}
