#include "../gplspec/src/GeneralModels/Model_Construction_Utilities.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using ScalarField = std::vector<std::vector<std::vector<double>>>;

void require(bool condition, const std::string &message) {
   if (!condition) {
      throw std::runtime_error(message);
   }
}

struct MockRadialMesh {
   int LayerNumber(int element) const { return element; }
   double NodeRadius(int element, int node) const {
      return 0.35 + 0.45 * element + 0.3 * node;
   }
   double PlanetRadius() const { return 1.0; }
   double OuterRadius() const { return 1.6; }
};

struct MockGrid {
   const std::vector<double> &CoLatitudes() const { return colatitudes; }
   const std::vector<double> &Longitudes() const { return longitudes; }
   std::vector<double> colatitudes{0.2, 1.1};
   std::vector<double> longitudes{0.0, 1.3, 2.2};
};

struct MockMapping {
   auto RadialMapping(int layer) const {
      ++requests;
      return [layer](double radius, double colatitude, double longitude) {
         return 0.5 * (layer + 1) * radius + 0.25 * colatitude -
                0.125 * longitude;
      };
   }
   mutable std::size_t requests = 0;
};

void original_mapping_geometry(const MockRadialMesh &node_data,
                               const MockGrid &grid, const MockMapping &mapping,
                               int num_layers, int polynomial_order,
                               ScalarField &displacement) {
   for (int idxelem = 0; idxelem < num_layers; ++idxelem) {
      int laynum = node_data.LayerNumber(idxelem);
      for (int idxnode = 0; idxnode < polynomial_order + 1; ++idxnode) {
         auto radr = node_data.NodeRadius(idxelem, idxnode);
         auto multfact = 1.0;
         auto raduse = radr;
         if (radr > node_data.PlanetRadius()) {
            raduse = node_data.PlanetRadius();
            multfact = (node_data.OuterRadius() - radr) /
                       (node_data.OuterRadius() - node_data.PlanetRadius());
         }
         int idxspatial = 0;
         for (auto it : grid.CoLatitudes()) {
            for (auto ip : grid.Longitudes()) {
               displacement[idxelem][idxnode][idxspatial] =
                   mapping.RadialMapping(laynum)(raduse, it, ip) * multfact;
               ++idxspatial;
            }
         }
      }
   }
}
}  // namespace

int main() {
   constexpr int layers = 2;
   constexpr int polynomial_order = 2;
   const MockRadialMesh radial_mesh;
   const MockGrid grid;
   const MockMapping original_mapping, shared_mapping;
   auto original = GPLSpec::model_detail::InitializeScalarFieldStorage(
       layers, polynomial_order + 1, grid.colatitudes.size() * grid.longitudes.size(), 0.0);
   auto actual = original;
   original_mapping_geometry(radial_mesh, grid, original_mapping, layers,
                             polynomial_order, original);
   GPLSpec::model_detail::PopulateRadialMappingGeometry(
       radial_mesh, grid, shared_mapping, layers, polynomial_order, actual);
   require(actual == original, "radial mapping geometry differs from original loops");
   require(shared_mapping.requests == original_mapping.requests,
           "radial mapping call count/order changed");
   require(shared_mapping.requests == layers * (polynomial_order + 1) *
                                          grid.colatitudes.size() *
                                          grid.longitudes.size(),
           "radial mapping was not sampled once per spatial node");
}
