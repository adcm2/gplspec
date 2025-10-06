---
layout: default
title: "gplspec Documentation"
---

# gplspec: Gravity and Planetary Spectra Solver

Welcome to the documentation for `gplspec`, a C++ library for solving problems related to gravity and planetary spectra.

## Quick Navigation

- **[Getting Started]({{ '/getting-started/' | relative_url }})**: Installation and first steps
- **[Tutorials]({{ '/tutorials/' | relative_url }})**: Learn the core concepts
- **[Examples]({{ '/examples/' | relative_url }})**: Code examples and benchmarks
- **[API Reference]({{ '/api/' | relative_url }})**: Complete class and function reference

## Quick Example

Here's a simple example of using gplspec to compute gravitational potential:

```cpp
#include <gplspec/All>
using namespace GeneralEarthModels;
using namespace Gravity_Tools;

int main() {
    // Create a density model from PREM
    Density3D model = Density3D::OneDimensionalPlanetFromFile(
        "prem.200", 5, 2, 0.01, 1.2);
    
    // Compute gravitational potential
    auto potential = FindGravitationalPotential(model);
    
    return 0;
}
```