#ifndef GPLSPEC_MODEL_CONSTRUCTION_UTILITIES_H
#define GPLSPEC_MODEL_CONSTRUCTION_UTILITIES_H

#include <cstddef>
#include <vector>

namespace GPLSpec::model_detail {

// Allocate an element/node/spatial scalar field with the requested initial value.
inline std::vector<std::vector<std::vector<double>>> InitializeScalarFieldStorage(
    std::size_t element_count, std::size_t node_count,
    std::size_t spatial_count, double initial_value) {
   return std::vector<std::vector<std::vector<double>>>(
       element_count,
       std::vector<std::vector<double>>(
           node_count, std::vector<double>(spatial_count, initial_value)));
}

template <class RadialMesh, class TransformGrid, class Mapping, class ScalarField>
void PopulateRadialMappingGeometry(const RadialMesh &node_data,
                                   const TransformGrid &grid,
                                   const Mapping &mapping, int num_layers,
                                   int polynomial_order, ScalarField &displacement) {
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

}  // namespace GPLSpec::model_detail

#endif
