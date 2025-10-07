---
title: "Tutorial 4: Topography  "
description: "Gravitational potential calculation for aspherical"
layout: default
---

# Surface topography
## Outline
This tutorial demonstrates how to construct a density model that contains a mapping $$\xi$$. We choose to construct a model that has a non-zero mapping but is of a spherical, homogeneous planet. We specify the physical properties of the planet, such as density and radius, as well as the mapping. The mapping needs to have a function called 
```cpp
RadialMapping(int i)
```
This specifies the radial mapping within the $i$-th layer. It must return a function that can take in three doubles, for the radius, colatitude and longitude. An example of the radial mapping class is
```cpp
class inp_map {
    public:
      inp_map() {};
      inp_map(const double h) : _h{h} {};
      auto RadialMapping(int i) const {
         auto lambdamap = [hmult = _h](double r, double theta, double phi) {
            return hmult * r;
         };
         return lambdamap;
      }

    private:
      double _h = 0.0;
   };
```


The model that is generated will have the following properties
* Physically it will be a homogeneous sphere with density and radius as specified
* The reference model will not necessarily have spherical internal boundaries 
* The mapping will be that specified by the user


### Source code
```cpp
#include <gplspec/All>
int
main() {

   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;
   using namespace SimpleModels;
   using FLOAT = double;

   // parameters:
   double maxstep = 0.1, ballrad = 1.4, radius = 1.0, density = 1.0;
   int npoly = 5, lmax = 2;
   double h = 0.2;
   double lengthnorm = EarthModels::EarthConstants<FLOAT>().LengthNorm();
   double timenorm = EarthModels::EarthConstants<FLOAT>().TimeNorm();
   double massnorm = EarthModels::EarthConstants<FLOAT>().MassNorm();

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////
   // declare model and find the potential
   // mapping class
   class inp_map {
    public:
      inp_map() {};
      inp_map(const double h) : _h{h} {};
      auto RadialMapping(int i) const {
         auto lambdamap = [hmult = _h](double r, double theta, double phi) {
            return hmult * r;
         };
         return lambdamap;
      }

    private:
      double _h = 0.0;
   };

   // construct radial displacement
   inp_map map1(h);

   // sphere
   spherical_model inp_sphere = spherical_model::HomogeneousSphere(
       radius, density, lengthnorm, timenorm, massnorm);

   // construct sphere that has a non-zero radial displacement
   Density3D homogeneous_sphere_with_map(inp_sphere,
                                         PlanetaryModel::TomographyZeroModel(),
                                         map1, npoly, lmax, maxstep, ballrad);

   // find potential
   auto stdvec_potsol = FindGravitationalPotential(homogeneous_sphere_with_map,
                                                   std::pow(10.0, -12.0));
   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // benchmarking using an equivalent spherical homogeneu
   //  equivalent spherical model
   Density3D equivalent_sphere = Density3D::SphericalHomogeneousPlanet(
       radius * (1.0 + h), density / std::pow(1.0 + h, 3.0), lengthnorm,
       timenorm, massnorm, maxstep, ballrad);

   // find radii at which the potential is evaluated
   auto inp_radii = homogeneous_sphere_with_map.PhysicalRadius_Line(0);

   // find the solution to the equivalent problem
   auto vec_exactsol = HomogeneousSphereIntegral(equivalent_sphere, inp_radii);

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
   std::string pathtofolder = "./work/Bench3";
   homogeneous_sphere_with_map.PhysicalOutputAtElement(pathtofolder,
                                                       stdvec_potsol);
   equivalent_sphere.OutputAtElement(pathtofolder, inp_radii, vec_exactsol);

   return 0;
}


```





