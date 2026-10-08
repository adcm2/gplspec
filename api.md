---
title: "API"
layout: default
permalink: /api/
---

# Supported API overview

GPLSpec currently exposes a header-only C++ interface through the CMake `gplspec` target. It requires C++23. There is no install/export package or promise that every internal header is a stable API; use the umbrella headers below and treat direct `gplspec/src/` includes as implementation details.

## Umbrella headers

| Header | Intended contents |
|---|---|
| [`gplspec/All`](gplspec/All) | General and simple models, mapping, matrix wrappers, spectral-element tools, preconditioning, integration, and timer utilities. |
| [`gplspec/GeneralModels`](gplspec/GeneralModels) | General density/model constructors, model return and perturbation tools, and gravity calculations. |
| [`gplspec/SimpleModels`](gplspec/SimpleModels) | Simple model definition, construction, and return helpers. |
| [`gplspec/Core`](gplspec/Core) | Mapping, pseudospectral matrix, preconditioner, integration, and timer components. |
| [`gplspec/Test`](gplspec/Test) | Model, matrix, radial, and test utilities used by the repository's examples and checks. |
| [`gplspec/Timer`](gplspec/Timer) | Timer utility. |

A linked target should include a supported umbrella, for example:

```cpp
#include <gplspec/GeneralModels>
```

The build and integration contract is documented in [Getting started](getting-started.md). The regression checks include separate compile probes for the umbrella headers, but those probes establish compilation only.

## Model construction and gravity routines

The model constructors in [`Density_Model_Constructors.h`](gplspec/src/GeneralModels/Density_Model_Constructors.h) include `Density3D::SphericalHomogeneousPlanet`, `Density3D::OneDimensionalPlanetFromFile`, and `Density3D::SphericalThreeDimensionalPlanetFromFile`. Their model layouts, units, file formats, and parameter conventions are defined by the corresponding declarations and implementations in `gplspec/src/GeneralModels/` and its pinned model dependencies.

[`Gravity_Tools.h`](gplspec/src/GeneralModels/Gravity_Tools.h) contains the public gravity operations used by the examples, including `FindGravitationalPotential` overloads, `GravitationalSphericalIntegral`, `HomogeneousSphereIntegral`, `SphericalHarmonicSensitivityKernel`, `FindGravitationalPotentialPerturbation`, and `FindGravitationalPotentialClassicalPerturbation`. The solver and perturbation routines have overload-specific tolerance and truncation arguments; consult their declarations at the linked source before relying on defaults. Example parameter choices and generated file contracts are mapped in [Scientific workflows](benchmark.md).

These named interfaces are a source-oriented overview, not a guarantee that all internal types or overloads are stable. Phase 3 recommendations for narrowing or retiring overlapping APIs will be recorded separately; no API retirement is implemented in Phase 2.
