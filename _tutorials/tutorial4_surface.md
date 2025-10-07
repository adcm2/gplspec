---
title: "Tutorial 4: Topography  "
description: "Gravitational potential calculation for aspherical"
layout: default
---

# Surface topography

This tutorial demonstrates how to construct a density model that contains a mapping $$\xi$$. We choose to construct a model that has a non-zero mapping but is of a spherical, homogeneous planet. 

## Code
### Model
To construct a 3D model we will use the constructor 
```cpp 
SphericalThreeDimensionalPlanetFromFile
```
This requires a path to a radial model, a path to the tomography model as well as radial parameters. Optionally one can specify norms, but the "default" values for the norms are those used in PREM so in this case we do not need to change these. 

### Sensitivity kernel
The sensitivity kernel functionality allows one to find the sensitivity kernel for a specific combination of spherical harmonic components of the gravitational potential on the ball surrounding the planet. This is given by 

$$Q(\zeta) = b^{-2} \int_{\partial \mathcal{B}} \sum_{lm} q_{lm}\overline{Y}_{lm} \zeta dS,$$

where the $q_{lm}$ are the relative components of the sensitivity kernel. To specify this, one simply needs to provide three vectors, with the degree l, the order m and the coefficient $q_{lm}$. A single function call will then give the sensitivity kernel in the standard output format.

### Source code
In this code we use PREM as the 1D reference model and S40RTS as the tomography model. 
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

## Results
### Output Files

The benchmark generates several output files in the `./work/Bench6/` directory:

- **3D Solution**: Spectral element method results
- **Integral Solution**: Spherical integration method results







