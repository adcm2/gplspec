---
title: "Tutorial 1: Homogeneous Sphere "
description: "Gravitational potential calculation for a homogeneous sphere"
layout: default
---

# Homogeneous Sphere 

This tutorial demonstrates how to construct a homogeneous sphere model and then the different ways to calculate its potential.


## Construction
### Outline
The class which is used throughout the code base is the Density3D class. It possesses multiple constructors for different situations. To construct a homogeneous sphere one needs to specify the norms associated with the model as well as its normalised radius and density. There are several methods available to find the gravitational potential, discussed below. There are also functions available to output the results. 

### Source code

```cpp

#include <gplspec/All>

int
main() {
   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;

   // Physical and model parameters
   // Normalised radius, density, length, time and mass normalisations, the
   // maximum mesh step size and the ball radius
   double nrad = 1.0, nrho = 1.0, lengthnorm = 6371000.0, timenorm = 3600.0,
          massnorm = 5.972e24, maxstep = 0.01, ballrad = 1.2;     

   // Construct homogeneous sphere model
   Density3D testsphere = Density3D::SphericalHomogeneousPlanet(
       nrad, nrho, lengthnorm, timenorm, massnorm, maxstep, ballrad);

   // Find potential via the pseudospectral method
   auto stdvec_potsol = FindGravitationalPotential(testsphere);

   // Compute potential using spherical integration method
   auto vec_integral_potential = GravitationalSphericalIntegral(testsphere);

   // Compute exact solution for homogeneous sphere
   auto vec_exactsol = HomogeneousSphereIntegral(testsphere);

   // Output results to files
   std::string pathtofolder = "./work/Bench1";
   testsphere.ReferentialOutputAtElement(pathtofolder, stdvec_potsol);
   testsphere.OutputAtElement(pathtofolder, vec_integral_potential);
   testsphere.OutputAtElement(pathtofolder, vec_exactsol);

   return 0;
}
```

## Results
### Output Files

The benchmark generates several output files in the `./work/Bench1/` directory:

- **3D Solution**: Spectral element method results
- **Integral Solution**: Spherical integration method results  
- **Exact Solution**: Analytical homogeneous sphere solution







