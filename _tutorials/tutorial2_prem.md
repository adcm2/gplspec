---
title: "Tutorial 2: PREM "
description: "Gravitational potential calculation for PREM"
layout: default
---

# PREM

This tutorial demonstrates how to construct a density model from PREM and the pseudospectral and radial integral methods for potential calculation.


## Construction
### Outline
To construct a radial model we will use the named constructor 
```cpp 
OneDimensionalPlanetFromFile 
```
This requires a path to a radial model, as well as parameters for the mesh. Optionally one can specify norms, but the ``default'' values for the norms are those used in PREM so in this case we do not need to change these. 

### Source code

```cpp

#include <gplspec/All>
#include <filesystem>
#include <utility>

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
   auto solve_result = SolvePotential(testprem);
   auto potential = ToHarmonicField(testprem, solve_result);

   // Integral method
   auto vec_integral_potential = GravitationalSphericalIntegral(testprem);


   auto field_metadata = GPLSpec::FieldBundle::WithSolve(
       solve_result, "clean_bench_2", "prem_200", "CG",
       solve_result.final_solve().requested_tolerance);
   field_metadata.inputs.push_back(GPLSpec::FieldBundle::InputPath(pathtoprem));
   field_metadata.parameters = {{"polynomial_order", std::to_string(npoly)},
       {"degree_max", std::to_string(lmax)},
       {"max_step_norm", GPLSpec::FieldBundle::ParameterNumber(maxstep)},
       {"exterior_radius_norm", GPLSpec::FieldBundle::ParameterNumber(ballrad)}};
   const std::filesystem::path output{"work/Bench2"};
   GPLSpec::FieldBundle::Writer(output / "field_bundle", testprem,
       field_metadata).Write([&](auto &&emit) {
      GPLSpec::FieldEvaluation::ForEachAngularGridSample(
          testprem, solve_result.solution,
          GPLSpec::FieldEvaluation::Quantity::Potential, emit);
      GPLSpec::FieldEvaluation::ForEachRadialHarmonicSample(
          testprem, solve_result.solution,
          GPLSpec::FieldEvaluation::Quantity::Potential, emit);
   });
   auto integral_metadata = GPLSpec::FieldBundle::RunMetadata{};
   integral_metadata.workflow = "clean_bench_2";
   integral_metadata.dataset = "prem200_spherical_integral";
   integral_metadata.inputs = field_metadata.inputs;
   integral_metadata.parameters = field_metadata.parameters;
   GPLSpec::FieldBundle::Writer(output / "integral_bundle", testprem,
       std::move(integral_metadata)).Write([&](auto &&emit) {
      GPLSpec::FieldEvaluation::ForEachElementBoundaryHarmonics(
          testprem, vec_integral_potential,
          GPLSpec::FieldEvaluation::Quantity::Potential, emit);
   });

   return 0;
}
```

## Results
### Output Files

Run from an isolated working directory containing the unchanged input at
`modeldata/prem.200`. The program writes schema-v1 field and integral bundles
under `work/Bench2/`.

- **3D Solution**: Spectral element method results
- **Integral Solution**: Spherical integration method results  






