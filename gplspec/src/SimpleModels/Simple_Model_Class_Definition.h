#ifndef SIMPLE_MODEL_CLASS_DEFINITION_H
#define SIMPLE_MODEL_CLASS_DEFINITION_H

#include <PlanetaryModel/All>

namespace GeneralEarthModels {
namespace SimpleModels {
class spherical_model {
 public:
   // default
   spherical_model() {};

   // spherical_model(double, double, double, double, double);
   static spherical_model HomogeneousSphere(double, double, double, double,
                                            double);
   static spherical_model HomogeneousLayers(std::vector<double> &,
                                            std::vector<double> &, double,
                                            double, double);

   // return functions
   double LengthNorm() const;
   double MassNorm() const;
   double TimeNorm() const;
   double DensityNorm() const;
   double InertiaNorm() const;
   double VelocityNorm() const;
   double AccelerationNorm() const;
   double ForceNorm() const;
   double StressNorm() const;
   int NumberOfLayers() const;
   auto LowerRadius(int i) const;
   auto UpperRadius(int i) const;
   auto OuterRadius() const;
   auto Density(int i) const;

 private:
   double _length_norm, _time_norm, _mass_norm;
   int _number_of_layers = 1;
   std::vector<double> _vec_layer_boundaries, _vec_layer_densities;

   // general constructor
   spherical_model(std::vector<double> &, std::vector<double> &, double, double,
                   double);
};

}   // namespace SimpleModels

}   // namespace  GeneralEarthModels
#endif