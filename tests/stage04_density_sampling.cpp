#include "../gplspec/src/GeneralModels/Model_Construction_Utilities.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition, const std::string &message) {
   if (!condition) {
      throw std::runtime_error(message);
   }
}

struct MockGrid {
   const std::vector<double> &CoLatitudes() const { return colatitudes; }
   const std::vector<double> &Longitudes() const { return longitudes; }
   std::vector<double> colatitudes{0.1, 0.9, 2.2};
   std::vector<double> longitudes{0.0, 1.7};
};

struct MockTomography {
   const std::vector<double> &GetDepths() const { return depths; }
   double GetValueAt(double depth, double longitude, double latitude) const {
      ++samples;
      return 1.0e-5 * depth + 2.0e-3 * longitude - 3.0e-3 * latitude;
   }
   std::vector<double> depths{0.0, 1500.0, 3000.0};
   mutable std::size_t samples = 0;
};

std::vector<double> original_density_update(double density_1d,
                                            std::size_t spatial_count,
                                            double radius,
                                            double depth_reference_radius,
                                            double length_norm,
                                            const MockTomography &tomography,
                                            const MockGrid &grid) {
   std::vector<double> density(spatial_count, density_1d);
   auto depth = depth_reference_radius - radius;
   depth *= length_norm / 1000.0;
   if (depth > tomography.GetDepths()[0] &&
       depth < tomography.GetDepths().back()) {
      int idxspatial = 0;
      for (auto idxt : grid.CoLatitudes()) {
         double pi_db = 3.1415926535897932384626433;
         double latitude = pi_db / 2.0 - idxt;
         latitude *= 180.0 / pi_db;
         for (auto idxp : grid.Longitudes()) {
            double multfact = 180.0 / 3.1415926535897932384626433;
            auto longitude = multfact * idxp;
            density[idxspatial] *=
                (1.0 + 0.005 * tomography.GetValueAt(depth, longitude, latitude));
            ++idxspatial;
         }
      }
   }
   return density;
}

void check_case(double radius, double depth_reference_radius,
                const MockGrid &grid) {
   constexpr double initial_density = 2.5;
   constexpr double length_norm = 1000.0;
   const std::size_t spatial_count =
       grid.colatitudes.size() * grid.longitudes.size();
   const MockTomography original_tomography, shared_tomography;
   const auto expected = original_density_update(
       initial_density, spatial_count, radius, depth_reference_radius,
       length_norm, original_tomography, grid);
   std::vector<double> actual(spatial_count, initial_density);
   auto depth = depth_reference_radius - radius;
   depth *= length_norm / 1000.0;
   GPLSpec::model_detail::ApplyReferentialTomographyDensityVariation(
       actual, depth, shared_tomography, grid);
   require(actual == expected, "referential tomography density differs from original loop");
   require(shared_tomography.samples == original_tomography.samples,
           "tomography lookup count changed");
}
}  // namespace

int main() {
   const MockGrid grid;
   check_case(0.4, 1.6, grid);
   check_case(0.4, 1.0, grid);
   check_case(1.6, 1.6, grid);  // strict lower-depth boundary
   check_case(-1.4, 1.6, grid); // strict upper-depth boundary
   check_case(2.0, 1.0, grid);  // outside the sampled depth interval
}
