---
title: "Home"
---

# gplspec: Gravitational potential via pseudospectral particle relabelling and radial spectral elements

Welcome to the documentation for `gplspec`, a C++ library for finding the gravitational potential of arbitrary, heterogeneous planets. The original theory paper for this method is available at https://academic.oup.com/gji/article/219/2/1043/5541065 whilst further theoretical developments and implementation details for gplspec are discussed in https://arxiv.org/abs/2508.07910.

## Overview

`gplspec` implements a novel numerical method that combines the **particle-relabelling transformation** with **radial spectral elements** to solve Poisson's equation for gravitational potential. This approach provides several key advantages:

### Core Methodology

- **Particle-Relabelling Transformation**: Maps aspherical planets onto a spherical reference planet. Only requirement is concentric, non-intersecting layers
- **Radial Spectral Elements**: High-order Gauss-Lobatto-Legendre basis functions provide exponential convergence for smooth solutions
- **Spherical Harmonic Decomposition**: Spherical geometry of reference planet enables usage of spherical harmonic decomposition
- **Matrix-Free Implementation**: Avoids storage of large system matrices through efficient operator applications

### Technical Features

- **Arbitrary Heterogeneity**: Handles complex 3D density variations including discontinuities and rapid spatial changes
- **Spectral Accuracy**: Achieves machine precision for smooth density distributions with relatively few degrees of freedom
- **Computational Efficiency**: Scales as O(L^3) for maximum spherical harmonic degree L

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

## Performance Characteristics

The particle-relabelling method offers significant computational advantages:

- **Convergence Rate**: Exponential convergence for smooth density distributions
- **Memory Efficiency**: Matrix-free implementation reduces memory requirements by orders of magnitude
- **Accuracy**: Machine precision achievable for well-resolved problems

## Getting Help

- Browse the [tutorials]({{ '/tutorials/' | relative_url }}) for guided learning
- Check out [examples]({{ '/examples/' | relative_url }}) for practical applications
- Refer to the [API documentation]({{ '/api/' | relative_url }}) for detailed function references

## Benchmarks

Explore our performance benchmarks and validation tests:

{% for benchmark in site.benchmarks limit:3 %}
- [{{ benchmark.title }}]({{ benchmark.url | relative_url }})
{% endfor %}

[View all benchmarks →]({{ '/benchmarks/' | relative_url }})

## Citation

If you use gplspec in your research, please cite:

- Myhill, A. D., Maitra, M. A., & Al-Attar, D. (2025). Forward and adjoint calculations of gravitational potential in heterogeneous, aspherical planets. arXiv preprint arXiv:2508.07910.
- Maitra, M., & Al-Attar, D. (2019). A non-perturbative method for gravitational potential calculations within heterogeneous and aspherical planets. Geophysical Journal International, 219(2), 1043-1055.

## License

This project is open source. See the repository for license details.