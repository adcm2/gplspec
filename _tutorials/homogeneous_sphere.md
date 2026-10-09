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
#include <filesystem>
#include <utility>

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
   auto solve_result = SolvePotential(testsphere);
   auto potential = ToHarmonicField(testsphere, solve_result);

   // Compute potential using spherical integration method
   auto vec_integral_potential = GravitationalSphericalIntegral(testsphere);

   // Compute exact solution for homogeneous sphere
   auto vec_exactsol = HomogeneousSphereIntegral(testsphere);

   // Write the solved field as a schema-v1 harmonic bundle.
   auto metadata = GPLSpec::FieldBundle::WithSolve(
       solve_result, "homogeneous_sphere", "homogeneous_sphere", "CG",
       solve_result.final_solve().requested_tolerance);
   metadata.parameters = {{"normalized_radius", GPLSpec::FieldBundle::ParameterNumber(nrad)},
       {"normalized_density", GPLSpec::FieldBundle::ParameterNumber(nrho)},
       {"length_norm_m", GPLSpec::FieldBundle::ParameterNumber(lengthnorm)},
       {"time_norm_s", GPLSpec::FieldBundle::ParameterNumber(timenorm)},
       {"mass_norm_kg", GPLSpec::FieldBundle::ParameterNumber(massnorm)},
       {"max_step_norm", GPLSpec::FieldBundle::ParameterNumber(maxstep)},
       {"exterior_radius_norm", GPLSpec::FieldBundle::ParameterNumber(ballrad)}};
   const auto output = std::filesystem::path("work/Bench1/field_bundle");
   GPLSpec::FieldBundle::WriteRadialHarmonics(
       output, testsphere, solve_result.solution,
       GPLSpec::FieldEvaluation::Quantity::Potential, std::move(metadata));

   return 0;
}
```

## Results
### Output Files

Run this program with its working directory set to an isolated run directory.
It writes a schema-v1 `work/Bench1/field_bundle` containing the solved
potential's radial harmonic coefficients. The integral and exact calculations
above are computed in memory; this short example does not serialize them.

- **Field bundle**: Spectral element potential coefficients






