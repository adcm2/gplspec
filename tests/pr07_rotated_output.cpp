#include "../gplspec/All"

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition, const std::string &message) {
   if (!condition) throw std::runtime_error(message);
}

struct TwoLayerModel {
   int NumberOfLayers() const { return 2; }
   double LowerRadius(int layer) const { return layer == 0 ? 0.0 : 0.5 * layer; }
   double UpperRadius(int layer) const { return 0.5 * (layer + 1); }
   double OuterRadius() const { return 1.0; }
   double LengthNorm() const { return 2.0; }
   double MassNorm() const { return 64.0; }
   double TimeNorm() const { return 1.0; }
   auto Density(int layer) const {
      const double rho = layer == 0 ? 4.0 : 10.0;
      return [rho](double) { return rho; };
   }
};

struct RadialScaleMapping {
   auto RadialMapping(int) const {
      return [](double radius, double, double) { return 0.1 * radius; };
   }
};

std::vector<std::vector<double>> read_rows(const std::filesystem::path &path) {
   std::ifstream input(path);
   require(static_cast<bool>(input), "could not open rotated output");
   std::vector<std::vector<double>> rows;
   std::string line;
   while (std::getline(input, line)) {
      std::replace(line.begin(), line.end(), ';', ' ');
      std::istringstream stream(line);
      std::vector<double> values;
      double value;
      while (stream >> value) values.push_back(value);
      rows.push_back(std::move(values));
   }
   return rows;
}

void near(double actual, double expected, const std::string &message) {
   require(std::abs(actual - expected) <= 1.0e-10 * std::max(1.0, std::abs(expected)),
           message);
}
}  // namespace

int main(int argc, char **argv) {
   using GeneralEarthModels::Density3D;
   require(argc == 2, "usage: pr07_rotated_output OUTPUT_DIR");
   const std::filesystem::path output_dir(argv[1]);
   std::filesystem::create_directories(output_dir);

   const TwoLayerModel source;
   const PlanetaryModel::TomographyZeroModel tomography;
   Density3D model(source, tomography, RadialScaleMapping{}, 2, 2, 0.25, 1.2);
   near(model.DensityNorm(), 8.0, "fixture density norm changed");
   near(model.PotentialNorm(), 4.0, "fixture potential norm changed");

   std::vector<double> p1{0.0, 0.0}, p2{1.5707963267948966, 0.0};
   const auto elements = static_cast<std::size_t>(model.Num_Elements());
   const auto poly = static_cast<std::size_t>(model.Poly_Order());
   std::vector<std::vector<std::vector<std::complex<double>>>> discontinuous(
       elements,
       std::vector<std::vector<std::complex<double>>>(
           poly + 1, std::vector<std::complex<double>>(9)));
   std::vector<std::vector<std::vector<std::complex<double>>>> continuous = discontinuous;
   for (std::size_t elem = 0; elem < elements; ++elem) {
      for (std::size_t node = 0; node <= poly; ++node) {
         discontinuous[elem][node][0] = {100.0 * elem + node + 1.0, 0.0};
         continuous[elem][node][0] = {3.0, 0.0};
      }
   }
   model.ReferentialOutputRotated((output_dir / "referential.csv").string(),
                                  p1, p2, discontinuous);
   model.PhysicalOutputRotated((output_dir / "physical.csv").string(), p1, p2,
                               discontinuous);
   model.ReferentialOutputRotated((output_dir / "continuous.csv").string(), p1,
                                  p2, continuous);
   model.PhysicalOutputRotated((output_dir / "continuous-physical.csv").string(),
                               p1, p2, continuous);
   model.ModelDensityOutputRotated(
       (output_dir / "density-referential.csv").string(), p1, p2, false);
   model.ModelDensityOutputRotated(
       (output_dir / "density-physical.csv").string(), p1, p2, true);

   const auto ref_rows = read_rows(output_dir / "referential.csv");
   const auto physical_rows = read_rows(output_dir / "physical.csv");
   const auto continuous_rows = read_rows(output_dir / "continuous.csv");
   const auto continuous_physical_rows =
       read_rows(output_dir / "continuous-physical.csv");
   require(ref_rows.size() == elements + 1 && physical_rows.size() == elements + 1 &&
               continuous_rows.size() == elements + 1 &&
               continuous_physical_rows.size() == elements + 1,
           "rotated potential writers must retain first and final endpoints");
   constexpr double harmonic_scale = 0.28209479177387814;  // 1 / sqrt(4 pi)
   for (std::size_t elem = 0; elem < elements; ++elem) {
      require(ref_rows[elem].size() >= 5 && physical_rows[elem].size() >= 5 &&
                  continuous_rows[elem].size() >= 5 &&
                  continuous_physical_rows[elem].size() >= 5,
              "rotated writer emitted an incomplete row");
      const double trace = elem == 0 ? discontinuous[0][0][0].real()
                                     : discontinuous[elem - 1][poly][0].real();
      near(ref_rows[elem][3], trace * model.PotentialNorm() * harmonic_scale,
           "referential writer selected the wrong radial trace or norm");
      near(physical_rows[elem][3], trace * model.PotentialNorm() * harmonic_scale,
           "physical writer selected the wrong radial trace or norm");
      near(continuous_rows[elem][3], 3.0 * model.PotentialNorm() * harmonic_scale,
           "continuous referential potential changed at a duplicated-radius boundary");
      near(continuous_physical_rows[elem][3], 3.0 * model.PotentialNorm() * harmonic_scale,
           "continuous physical potential changed at a duplicated-radius boundary");
   }
   require(ref_rows.back().size() >= 5 && physical_rows.back().size() >= 5,
           "rotated writer emitted an incomplete final endpoint");
   near(ref_rows.back()[3], discontinuous.back()[poly][0].real() *
                                  model.PotentialNorm() * harmonic_scale,
        "referential writer changed the final endpoint");
   near(physical_rows.back()[3], discontinuous.back()[poly][0].real() *
                                      model.PotentialNorm() * harmonic_scale,
        "physical writer changed the final endpoint");

   const auto density_ref = read_rows(output_dir / "density-referential.csv");
   const auto density_phys = read_rows(output_dir / "density-physical.csv");
   require(density_ref.size() == elements + 1 && density_phys.size() == elements + 1,
           "density writer must retain first and final endpoints");
   bool saw_jump = false;
   bool saw_nonunit_jacobian = false;
   for (std::size_t elem = 0; elem < elements; ++elem) {
      require(density_ref[elem].size() >= 5 && density_phys[elem].size() >= 5,
              "density writer emitted an incomplete row");
      const int selected_elem = elem == 0 ? 0 : static_cast<int>(elem - 1);
      const int selected_node = elem == 0 ? 0 : static_cast<int>(poly);
      const double rho_ref = model.Density_Point(selected_elem, selected_node, 0);
      const double jacobian = model.Jacobian_Point(selected_elem, selected_node, 0);
      const double rho_phys = rho_ref / jacobian;
      if (std::abs(jacobian - 1.0) > 1.0e-8) saw_nonunit_jacobian = true;
      near(density_ref[elem][3], rho_ref * model.DensityNorm(),
           "density writer selected the wrong inner-side trace or density norm");
      near(density_phys[elem][3], rho_phys * model.DensityNorm(),
           "physical-density writer violated rho_ref = J rho_phys");
      if (elem > 0 &&
          std::abs(density_ref[elem][3] - density_ref[elem - 1][3]) > 1.0e-8)
         saw_jump = true;
   }
   require(density_ref.back().size() >= 5 && density_phys.back().size() >= 5,
           "density writer emitted an incomplete final endpoint");
   require(saw_jump, "fixture did not expose a discontinuous density trace");
   require(saw_nonunit_jacobian, "fixture must exercise a nonunit Jacobian");
   const double final_density = model.Density_Point(elements - 1, poly, 0);
   const double final_jacobian = model.Jacobian_Point(elements - 1, poly, 0);
   near(density_ref.back()[3], final_density * model.DensityNorm(),
        "density writer changed its final endpoint");
   near(density_phys.back()[3], final_density / final_jacobian * model.DensityNorm(),
        "physical-density writer changed its final endpoint scaling");
   require(model.Density_Point(static_cast<int>(elements - 1), 0, 0) == 0.0 &&
               density_ref[density_ref.size() - 2][3] > 0.0 &&
               density_phys[density_phys.size() - 2][3] > 0.0 && final_density == 0.0,
           "fixture must expose the inner body trace and outer vacuum endpoint");
}
