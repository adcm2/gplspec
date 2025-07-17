#ifndef SIMPLE_MODEL_RETURN_DEFINITION_H
#define SIMPLE_MODEL_RETURN_DEFINITION_H

namespace GeneralEarthModels {
namespace SimpleModels {

double
spherical_model::LengthNorm() const {
   return _length_norm;
};

double
spherical_model::MassNorm() const {
   return _mass_norm;
};
double
spherical_model::TimeNorm() const {
   return _time_norm;
}
double
spherical_model::DensityNorm() const {
   return _mass_norm / std::pow(_length_norm, 3.0);
};
double
spherical_model::InertiaNorm() const {
   return _mass_norm * std::pow(_length_norm, 2.0);
};
double
spherical_model::VelocityNorm() const {
   return _length_norm / _time_norm;
};
double
spherical_model::AccelerationNorm() const {
   return _length_norm / std::pow(_time_norm, 2.0);
};
double
spherical_model::ForceNorm() const {
   return _mass_norm * _length_norm / std::pow(_time_norm, 2.0);
};
double
spherical_model::StressNorm() const {
   return _mass_norm / (std::pow(_time_norm, 2.0) * _length_norm);
};

int
spherical_model::NumberOfLayers() const {
   return _number_of_layers;
};
auto
spherical_model::LowerRadius(int i) const {
   assert(i < _number_of_layers && "Not in model");
   return _vec_layer_boundaries[i];
};
auto
spherical_model::UpperRadius(int i) const {
   assert(i < _number_of_layers && "Not in model");
   return _vec_layer_boundaries[i + 1];
};
auto
spherical_model::OuterRadius() const {
   return _vec_layer_boundaries.back();
};

auto
spherical_model::Density(int i) const {
   assert(i < _number_of_layers && "Not in model");
   auto densitylambda = [intdensity = this->_vec_layer_densities[i]](
                            double eval_rad) { return intdensity; };
   return densitylambda;
};
}   // namespace SimpleModels
}   // namespace GeneralEarthModels

#endif