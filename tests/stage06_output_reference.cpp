#include "../gplspec/All"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct SmoothMapping {
   auto RadialMapping(int) const {
      return [](double radius, double theta, double phi) {
         return 0.005 * radius *
                (1.0 + 0.2 * std::sin(theta) * std::cos(phi));
      };
   }
};
}  // namespace

int main(int argc, char **argv) {
   if (argc != 2) {
      throw std::runtime_error("usage: stage06_output_reference OUTPUT_DIR");
   }
   const std::filesystem::path output_dir(argv[1]);
   std::filesystem::create_directories(output_dir / "referential");
   std::filesystem::create_directories(output_dir / "physical");

   auto model = GeneralEarthModels::Density3D(
       0.8, 2.0, 0.8, SmoothMapping{}, 2, 3, 1.0, 1.0, 1.0, 0.25, 1.2);
   const int lMax = model.GSH_Grid().MaxDegree();
   const std::size_t coefficients = static_cast<std::size_t>((lMax + 1) * (lMax + 1));
   std::vector<std::vector<std::vector<std::complex<double>>>> potential(
       model.Num_Elements(),
       std::vector<std::vector<std::complex<double>>>(
           model.Poly_Order() + 1, std::vector<std::complex<double>>(coefficients)));
   for (std::size_t elem = 0; elem < potential.size(); ++elem) {
      for (std::size_t node = 0; node < potential[elem].size(); ++node) {
         for (std::size_t lm = 0; lm < potential[elem][node].size(); ++lm) {
            potential[elem][node][lm] = {
                0.001 * static_cast<double>(1 + elem + 2 * node + lm),
                -0.0003 * static_cast<double>(2 + elem + node + lm)};
         }
      }
   }
   std::vector<double> p1{0.7, 0.2};
   std::vector<double> p2{1.1, 2.0};
   const auto path = [&output_dir](const char *name) {
      return (output_dir / name).string();
   };

   model.ReferentialOutputRotated(path("referential-rotated.csv"), p1, p2,
                                  potential);
   model.PhysicalOutputRotated(path("physical-rotated.csv"), p1, p2,
                               potential);
   model.ModelDensityOutputRotated(path("density-referential-rotated.csv"),
                                   p1, p2, false);
   model.ModelDensityOutputRotated(path("density-physical-rotated.csv"), p1,
                                   p2, true);
   model.ReferentialOutputAtElement((output_dir / "referential").string(),
                                    potential);
   model.PhysicalOutputAtElement((output_dir / "physical").string(),
                                 potential);

   const auto slice = model.RotateSliceToEquator(p1, p2, potential);
   std::ofstream slice_file(output_dir / "slice-rotation.csv");
   for (const auto &element : slice) {
      for (const auto &node : element) {
         for (const auto &coefficient : node) {
            slice_file << std::setprecision(17) << coefficient.real() << ';'
                       << coefficient.imag() << '\n';
         }
      }
   }
   return 0;
}
