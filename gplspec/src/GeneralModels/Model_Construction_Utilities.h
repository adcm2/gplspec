#ifndef GPLSPEC_MODEL_CONSTRUCTION_UTILITIES_H
#define GPLSPEC_MODEL_CONSTRUCTION_UTILITIES_H

#include <cstddef>
#include <GaussQuad/All>
#include <vector>

namespace GPLSpec::model_detail {

template <class TransformGrid>
void InitializeModelDiscretization(
    GaussQuad::Quadrature1D<double> &quadrature, TransformGrid &grid,
    int &polynomial_order, int requested_polynomial_order, int l_max) {
   quadrature =
       GaussQuad::GaussLobattoLegendreQuadrature1D<double>(
           requested_polynomial_order + 1);
   grid = TransformGrid(l_max, 2);
   polynomial_order = quadrature.N() - 1;
}

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

template <class Tomography, class TransformGrid>
void ApplyReferentialTomographyDensityVariation(
    std::vector<double> &density, double depth, const Tomography &tomography,
    const TransformGrid &grid) {
   if (depth > tomography.GetDepths()[0] &&
       depth < tomography.GetDepths().back()) {
      int idxspatial = 0;
      for (auto idxt : grid.CoLatitudes()) {
         double pi_db = 3.1415926535897932384626433;
         auto latitude = pi_db / 2.0 - idxt;
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
}

}  // namespace GPLSpec::model_detail

#endif
