#ifndef SPHERICAL_INTEGRATOR_H
#define SPHERICAL_INTEGRATOR_H

// Compatibility surface for the protected ellipticity example. The former
// spherical solver implementations have been retired; these mesh utilities
// are owned by Radial_Tools and forwarded here only until that external
// consumer can migrate independently.
#include "Radial_Tools.h"

namespace GravityFunctions {

template <typename FLOAT,
          template <typename, typename> class SphericalModel>
   requires PlanetaryModel::SphericalElasticModel<
       SphericalModel<FLOAT, int>>
auto Radial_Node(const SphericalModel<FLOAT, int> &model,
                 const FLOAT max_element_size) {
   return Radial_Tools::Radial_Node(model, max_element_size);
}

template <typename FLOAT,
          template <typename, typename> class SphericalModel>
   requires PlanetaryModel::SphericalElasticModel<
       SphericalModel<FLOAT, int>>
auto Radial_Node(const SphericalModel<FLOAT, int> &model,
                 const FLOAT max_element_size, const FLOAT outer_radius) {
   return Radial_Tools::Radial_Node(model, max_element_size, outer_radius);
}

template <typename FLOAT>
auto All_Node(const std::vector<FLOAT> &node_radii,
              const GaussQuad::Quadrature1D<FLOAT> &quadrature) {
   return Radial_Tools::All_Node(node_radii, quadrature);
}

}   // namespace GravityFunctions

#endif
