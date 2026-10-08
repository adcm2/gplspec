#include <gplspec/GeneralModels>
#include <iomanip>
#include <iostream>

int main(int argc, char **argv) {
  const char *path = argc > 1 ? argv[1] : "modeldata/PREM500";
  TestTools::EarthModel<> model(path);
  std::cout << std::setprecision(17);
  std::cout << "META," << model.NumberOfLayers() << ',' << model.OuterRadius()
            << ',' << model.LengthNorm() << ',' << model.DensityNorm() << '\n';
  for (int layer = 0; layer < model.NumberOfLayers(); ++layer) {
    const auto radii = model.LayerRadii(layer);
    std::cout << "L," << layer << ',' << model.LayerLowerIndex(layer) << ','
              << model.LayerUpperIndex(layer) << ',' << model.IsSolid(layer)
              << ',' << radii.size() << '\n';
    for (std::size_t i = 0; i < radii.size(); ++i) {
      const double r = radii[i];
      std::cout << "N," << layer << ',' << i << ',' << r << ','
                << model.Density(layer)(r) << ',' << model.VPV(layer)(r)
                << ',' << model.VSV(layer)(r) << ',' << model.VPH(layer)(r)
                << ',' << model.VSH(layer)(r) << ',' << model.Eta(layer)(r)
                << '\n';
    }
  }
}
