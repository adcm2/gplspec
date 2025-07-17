#ifndef SIMPLE_MODEL_CONSTRUCTOR_DEFINITION_H
#define SIMPLE_MODEL_CONSTRUCTOR_DEFINITION_H

namespace GeneralEarthModels {
namespace SimpleModels {

// factory function that makes a homogeneous sphere
spherical_model
spherical_model::HomogeneousSphere(
    double radius, double density,
    double length_norm = EarthModels::EarthConstants<double>().LengthNorm(),
    double time_norm = EarthModels::EarthConstants<double>().TimeNorm(),
    double mass_norm = EarthModels::EarthConstants<double>().MassNorm()) {
   auto vec_radii = std::vector<double>(1, radius);
   auto vec_density = std::vector<double>(1, density);
   return spherical_model(vec_radii, vec_density, length_norm, time_norm,
                          mass_norm);
};

// factory function that makes a layered sphere
spherical_model
spherical_model::HomogeneousLayers(
    std::vector<double> &vec_layer_boundaries,
    std::vector<double> &vec_layer_densities,
    double length_norm = EarthModels::EarthConstants<double>().LengthNorm(),
    double time_norm = EarthModels::EarthConstants<double>().TimeNorm(),
    double mass_norm = EarthModels::EarthConstants<double>().MassNorm()) {
   return spherical_model(vec_layer_boundaries, vec_layer_densities,
                          length_norm, time_norm, mass_norm);
};

// full constructor
spherical_model::spherical_model(std::vector<double> &vec_layer_boundaries,
                                 std::vector<double> &vec_layer_densities,
                                 double length_norm, double time_norm,
                                 double mass_norm)
    : _length_norm(length_norm), _time_norm(time_norm), _mass_norm(mass_norm),
      _vec_layer_densities{vec_layer_densities} {

   // define things
   _number_of_layers = vec_layer_densities.size();
   assert(_number_of_layers == vec_layer_boundaries.size() &&
          "Incorrect dimensions");

   // layer boundaries
   _vec_layer_boundaries.resize(_number_of_layers + 1);
   _vec_layer_boundaries[0] = 0.0;
   for (int idx = 1; idx < _number_of_layers + 1; ++idx) {
      _vec_layer_boundaries[idx] = vec_layer_boundaries[idx - 1];
   }
};

}   // namespace SimpleModels
}   // namespace GeneralEarthModels

#endif