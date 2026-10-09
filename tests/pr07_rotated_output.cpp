#include "../gplspec/All"

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char *message) {
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
}  // namespace

int main(int argc, char **argv) {
   using GeneralEarthModels::Density3D;
   if (argc != 2) throw std::runtime_error("usage: pr07_rotated_output OUTPUT_DIR");
   std::filesystem::create_directories(argv[1]);

   const TwoLayerModel source;
   const PlanetaryModel::TomographyZeroModel tomography;
   Density3D model(source, tomography, RadialScaleMapping{}, 2, 2, 0.25, 1.2);
   require(std::abs(model.DensityNorm() - 8.0) < 1.0e-14 &&
               std::abs(model.PotentialNorm() - 4.0) < 1.0e-14,
           "fixture normalization changed");

   const int lmax = model.GSH_Grid().MaxDegree();
   const std::size_t stride = static_cast<std::size_t>((lmax + 1) * (lmax + 1));
   Eigen::VectorXcd solution(model.Node_Information().NumberOfGlobalRadialNodes() *
                             static_cast<Eigen::Index>(stride));
   solution.setZero();
   for (Eigen::Index node = 0; node < solution.size() /
                                         static_cast<Eigen::Index>(stride); ++node)
      solution[node * static_cast<Eigen::Index>(stride)] = 3.0;

   using namespace GPLSpec::FieldEvaluation;
   std::vector<Sample> outer;
   ForEachAngularGridSample(
       model, solution, Quantity::Potential,
       [&](const Sample &sample) { outer.push_back(sample); }, nullptr,
       TraceSelection::Outer);
   const std::size_t angular_count = model.GSH_Grid().Longitudes().size() *
                                     model.GSH_Grid().CoLatitudes().size();
   require(outer.size() == (static_cast<std::size_t>(model.Num_Elements()) + 1) *
                               angular_count,
           "single-side angular samples must retain each element boundary");
   for (const auto &sample : outer) {
      require(sample.side == TraceSide::None || sample.side == TraceSide::Outer,
              "outer cohort selected the wrong one-sided trace");
      require(std::isfinite(sample.value.real()) &&
                  std::isfinite(sample.physical_radius_m),
              "rotated field evaluator produced non-finite output");
   }

   const std::vector<double> azimuths{0.0, std::numbers::pi / 2.0};
   const Rotation rotation{0.0, std::numbers::pi / 2.0, 0.0};
   std::vector<Sample> slice;
   ForEachRotatedSliceSample(
       model, solution, Quantity::Potential, rotation, azimuths,
       [&](const Sample &sample) { slice.push_back(sample); });
   require(slice.size() == static_cast<std::size_t>(model.Num_Elements()) *
                               2 * azimuths.size(),
           "rotated slices must keep both radial endpoint traces");
   for (const auto &sample : slice) {
      require(sample.representation == Representation::RotatedSlice &&
                  std::isfinite(sample.value.real()) &&
                  sample.theta_rad == std::numbers::pi / 2.0,
              "rotated slice sample metadata is invalid");
   }

   std::vector<Sample> ref_density, physical_density;
   ForEachAngularGridSample(model, std::span<const std::complex<double>>{},
       Quantity::ReferentialDensity,
       [&](const Sample &sample) { ref_density.push_back(sample); });
   ForEachAngularGridSample(model, std::span<const std::complex<double>>{},
       Quantity::PhysicalDensity,
       [&](const Sample &sample) { physical_density.push_back(sample); });
   require(ref_density.size() == physical_density.size(),
           "density trace layouts must agree");
   bool saw_jump = false;
   for (std::size_t i = 0; i < ref_density.size(); ++i) {
      const auto &ref = ref_density[i];
      const auto &physical = physical_density[i];
      const double jacobian = model.Jacobian_Point(
          static_cast<int>(ref.element), ref.radial_node, 0);
      require(std::abs(physical.value.real() * jacobian -
                       ref.value.real()) < 1.0e-11,
              "physical and referential density lost their Jacobian relation");
      if (ref.boundary == BoundaryKind::Material) saw_jump = true;
   }
   require(saw_jump, "fixture must exercise a material density jump");
   return 0;
}
