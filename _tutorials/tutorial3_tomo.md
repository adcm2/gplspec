---
title: "Tutorial 3: S40RTS  "
description: "Gravitational potential calculation for PREM"
layout: default
---

# Tomography models and sensitivity kernels

This tutorial demonstrates how to construct a density model from PREM and how to add 3D variation via a tomography model. In addition the calculation of a sensitivity kernel is shown.

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

   // declaring variables:
   double maxstep = 0.05,  ballrad = 1.2;
   int npoly = 8, lmax = 20;
   std::string pathtoprem = "modeldata/prem.200";
   std::string pathtotomo = "modeldata/S40RTS_dvs.nc";


   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////
   // make model
   Density3D testtomo = Density3D::SphericalThreeDimensionalPlanetFromFile(
       pathtoprem, pathtotomo, npoly, lmax, maxstep, ballrad);
   
   // pseudospectral solution
   auto stdvec_potsol =
       FindGravitationalPotential(testtomo, std::pow(10.0, -14.0));

   // spherical integration
   auto vec_integral_potential = GravitationalSphericalIntegral(testtomo);


   // sensitivity kernel
   std::vector<int> vec_l{0, 1, 2, 3, 4}, vec_m{0, 0, 2, 3, 0};
   std::vector<double> multval{1.0, 3.0, 5.0, 7.0, 9.0};
   auto stdvec_senskernel =
       SphericalHarmonicSensitivityKernel(testtomo, vec_l, vec_m, multval);

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output
   std::string pathtofolder = "./work/Bench6";
   std::string pathtofolderlm = "./work/Bench6/lm";
   testtomo.PhysicalOutputAtElement(pathtofolder, stdvec_potsol);
   testtomo.PhysicalOutputAtElement(pathtofolder, vec_integral_potential);
   testtomo.ReferentialOutputAtElement(pathtofolderlm, stdvec_potsol);
   testtomo.OutputAtElement(pathtofolderlm, vec_integral_potential);
   testtomo.ReferentialOutputAtElement(pathtofolderlm, stdvec_senskernel, true);
   
   return 0;
}

```

## Results
### Output Files

The benchmark generates several output files in the `./work/Bench6/` directory:

- **3D Solution**: Spectral element method results
- **Integral Solution**: Spherical integration method results







