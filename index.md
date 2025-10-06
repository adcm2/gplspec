---
title: "Home"
---

# gplspec: Gravity and Planetary Spectra Solver

Welcome to the documentation for `gplspec`, a C++ library for solving problems related to gravity and planetary spectra using advanced numerical methods.

## Overview

`gplspec` provides efficient tools for:

- **Gravitational potential calculations** using spectral-element methods
- **Planetary density modeling** with radial basis functions
- **Spherical harmonic analysis** for global field representations
- **High-performance computing** with optimized numerical algorithms

## Quick Start

Get up and running with gplspec in just a few steps:

```cpp
#include <gplspec/All>
using namespace GeneralEarthModels;
using namespace Gravity_Tools;

int main() {
    // Load a planetary density model
    Density3D model = Density3D::OneDimensionalPlanetFromFile(
        "prem.200", 5, 2, 0.01, 1.2);
    
    // Compute gravitational potential
    auto potential = FindGravitationalPotential(model);
    
    return 0;
}
```

## Documentation Structure

This documentation is organized into several sections:

- **[Getting Started]({{ '/getting-started/' | relative_url }})** - Installation and basic setup
- **[Tutorials]({{ '/tutorials/' | relative_url }})** - Step-by-step guides to core concepts
- **[Examples]({{ '/examples/' | relative_url }})** - Complete working examples with explanations
- **[API Reference]({{ '/api/' | relative_url }})** - Detailed class and function documentation

## Key Features

### Spectral-Element Methods
High-order accuracy for gravitational potential calculations with efficient sparse matrix operations.

### Planetary Modeling
Support for 1D radial models (like PREM) extended to full 3D representations using spherical harmonics.

### Performance Optimization
Optimized algorithms designed for modern C++ with parallel computing support.

## Getting Help

- Browse the [tutorials]({{ '/tutorials/' | relative_url }}) for guided learning
- Check out [examples]({{ '/examples/' | relative_url }}) for practical applications
- Refer to the [API documentation]({{ '/api/' | relative_url }}) for detailed function references

## License

This project is open source. See the repository for license details.