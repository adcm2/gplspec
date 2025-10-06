---
title: "Homogeneous Sphere Benchmark"
description: "Gravitational potential calculation benchmark for a homogeneous sphere"
layout: default
---

# Homogeneous Sphere Benchmark

This benchmark demonstrates the calculation of gravitational potential for a homogeneous sphere using multiple methods to validate accuracy and performance.

## Overview

- **Purpose**: Validate gravitational potential calculations against analytical solutions
- **Methods**: 
  - Direct 3D spectral element calculation
  - Spherical integration method
  - Exact analytical solution
- **Model**: Homogeneous sphere with Earth-like parameters

## Simplified first example

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
### Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `nrad` | 1.0 | Normalized planetary radius |
| `nrho` | 1.0 | Normalized density |
| `lengthnorm` | 6,371,000 m | Earth radius normalization |
| `massnorm` | 5.972×10²⁴ kg | Earth mass normalization |
| `maxstep` | 0.01 | Maximum radial step size |
| `ballrad` | 1.2 | Computational domain radius |

## Results
### Output Files

The benchmark generates several output files in the `./work/Bench1/` directory:

- **3D Solution**: Spectral element method results
- **Integral Solution**: Spherical integration method results  
- **Exact Solution**: Analytical homogeneous sphere solution

These files can be used to:
- Compare numerical accuracy between methods
- Validate implementation correctness
- Analyze computational performance

### Comparison of methods

The benchmark produces gravitational potential solutions that can be compared for accuracy validation:

<div class="figure">
<img src="{{ '/benchfigures/Bench1.png' | relative_url }}" alt="Benchmark 1 Results" class="figure-small">
</div>

*Figure: Comparison of gravitational potential calculations for a homogeneous sphere using different numerical methods.*






