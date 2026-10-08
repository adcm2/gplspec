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
   inline static spherical_model HomogeneousSphere(double, double, double, double,
                                            double);
   inline static spherical_model HomogeneousLayers(std::vector<double> &,
                                            std::vector<double> &, double,
                                            double, double);

   // return functions
   inline double LengthNorm() const;
   inline double MassNorm() const;
   inline double TimeNorm() const;
   inline double DensityNorm() const;
   inline double InertiaNorm() const;
   inline double VelocityNorm() const;
   inline double AccelerationNorm() const;
   inline double ForceNorm() const;
   inline double StressNorm() const;
   inline int NumberOfLayers() const;
   inline auto LowerRadius(int i) const;
   inline auto UpperRadius(int i) const;
   inline auto OuterRadius() const;
   inline auto Density(int i) const;

 private:
   double _length_norm, _time_norm, _mass_norm;
   int _number_of_layers = 1;
   std::vector<double> _vec_layer_boundaries, _vec_layer_densities;

   // general constructor
   inline spherical_model(std::vector<double> &, std::vector<double> &, double, double,
                   double);
};

}   // namespace SimpleModels

}   // namespace  GeneralEarthModels
#endif