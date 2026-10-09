---
title: "Tutorial 5: Phobos  "
description: "Gravitational potential calculation for Phobos"
layout: default
---

# Phobos
## Outline
This tutorial demonstrates how to construct a model of Phobos. Phobos' density variation is not well constrained. However, its surface topography is relatively well known. The assumption made in the constructor for Phobos can be adapted to a more complex density distribution, but it assumes that the density is homogeneous by default. 

## Implementation
### Information
The model data is taken from Willner 2014. The surface radius is given in terms of a real spherical harmonic decomposition. The constructor used here assumes that the model data is in this format. The constructor needs the density, a path to the file and some more parameters. 

### Rotation
The compiled `phobos_gravity` workflow also writes original and rotated products;
the isolated producer and plotting commands are documented in `work/plots.md`.
The self-contained source below writes its original-frame potential and degree-power bundle.

### Source code
Run this example from an isolated working directory containing
`modeldata/TABLEA1.DAT`; it writes schema-v1 bundles under `work/Phobos/`.

```cpp
#include <cmath>
#include <filesystem>
#include <string>
#include <utility>
#include <gplspec/All>

int
main() {
   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;
   using namespace SimpleModels;
   using FLOAT = double;
   std::string pathtofile = "modeldata/TABLEA1.DAT";

   // parameters
   double physdensity = 1860.0, lengthnorm = 10000.0,  usedensity = 1342.0,
    massnorm = std::pow(lengthnorm, 3.0) * usedensity, timenorm = 3600.0,
     scaledensity = physdensity / massnorm * std::pow(lengthnorm, 3.0);
   int lMax = 100;
   int npoly = 5;

   // constructing model
   const auto phobos_data =
       GPLSpec::ModelInput::ReadPhobosHarmonics(pathtofile);
   Density3D phobos(scaledensity, phobos_data, npoly, lMax, lengthnorm, timenorm,
                    massnorm, 0.1, 1.5);

   // get gravitational field
   auto solve_result = SolvePotential(phobos, std::pow(10.0, -12.0));
   auto field_metadata = GPLSpec::FieldBundle::WithSolve(
       solve_result, "phobos_gravity", "potential", "CG",
       solve_result.final_solve().requested_tolerance);
   field_metadata.inputs.push_back(
       GPLSpec::FieldBundle::InputPath(pathtofile));
   field_metadata.parameters = {{"density_kg_m-3", GPLSpec::FieldBundle::ParameterNumber(physdensity)},
       {"length_norm_m", GPLSpec::FieldBundle::ParameterNumber(lengthnorm)},
       {"mass_norm_kg", GPLSpec::FieldBundle::ParameterNumber(massnorm)},
       {"time_norm_s", GPLSpec::FieldBundle::ParameterNumber(timenorm)},
       {"polynomial_order", std::to_string(npoly)}, {"degree_max", std::to_string(lMax)},
       {"max_step_norm", "0.1"}, {"exterior_radius_norm", "1.5"}};
   const std::filesystem::path output{"work/Phobos/Spherical"};
   GPLSpec::FieldBundle::Writer(output / "field_bundle", phobos,
       field_metadata).Write([&](auto &&emit) {
      GPLSpec::FieldEvaluation::ForEachAngularGridSample(
          phobos, solve_result.solution,
          GPLSpec::FieldEvaluation::Quantity::Potential, emit);
      GPLSpec::FieldEvaluation::ForEachRadialHarmonicSample(
          phobos, solve_result.solution,
          GPLSpec::FieldEvaluation::Quantity::DegreePower, emit);
   });
       
 
   return 0;
}


```

