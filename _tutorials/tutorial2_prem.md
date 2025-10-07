---
title: "Tutorial 2: PREM "
description: "Gravitational potential calculation for PREM"
layout: default
---

# Homogeneous Sphere 

This tutorial demonstrates how to construct a density model from PREM and the pseudospectral and radial integral methods for potential calculation.


## Construction
### Outline
To construct a radial model we will use the named constructor ```cpp OneDimensionalPlanetFromFile ```. This requires a path to a radial model, as well as parameters for the mesh. Optionally one can specify norms, but the ``default'' values for the norms are those used in PREM so in this case we do not need to change these. 

### Source code

```cpp

#include <gplspec/All>

int
main() {

   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;

   // parameters
   double maxstep = 0.01, ballrad = 1.2;
   int npoly = 5,  lmax = 2;

   // Path to the PREM model data file.
   std::string pathtoprem = "modeldata/prem.200";

   // Declare model and find the potential
   Density3D testprem = Density3D::OneDimensionalPlanetFromFile(
       pathtoprem, npoly, lmax, maxstep, ballrad);

   // Pseudospectral method
   auto stdvec_potsol = FindGravitationalPotential(testprem);

   // Integral method
   auto vec_integral_potential = GravitationalSphericalIntegral(testprem);


   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////
   // outputting results
   std::string pathtofolder = "./work/Bench2";
   testprem.ReferentialOutputAtElement(pathtofolder, stdvec_potsol);
   testprem.OutputAtElement(pathtofolder, vec_integral_potential);

   return 0;
}
```

## Results
### Output Files

The benchmark generates several output files in the `./work/Bench2/` directory:

- **3D Solution**: Spectral element method results
- **Integral Solution**: Spherical integration method results  







