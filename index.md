---
title: "Home"
---

# gplspec: Gravitational potential via pseudospectral particle relabelling and radial spectral elements

Welcome to the documentation for `gplspec`, a C++ library for finding the gravitational potential of arbitrary, heterogeneous planets. The original theory paper for this method is available at https://academic.oup.com/gji/article/219/2/1043/5541065 whilst further theoretical developments and implementation details for gplspec are discussed in https://arxiv.org/abs/2508.07910.

## Overview

`gplspec` implements a novel numerical method that combines the **particle-relabelling transformation** with **radial spectral elements** to solve Poisson's equation for gravitational potential. 

### Core Methodology

- **Particle-Relabelling Transformation**: Maps aspherical planets onto a spherical reference planet. Only requirement is concentric, non-intersecting layers
- **Radial Spectral Elements**: High-order Gauss-Lobatto-Legendre basis functions provide exponential convergence for smooth models
- **Spherical Harmonic Decomposition**: Spherical geometry of reference planet enables usage of spherical harmonic decomposition
- **Matrix-Free Implementation**: Avoids storage of large system matrices through efficient operator applications

### Technical features and performance

- **Arbitrary Heterogeneity**: Handles complex 3D density variations including radial discontinuities 
- **Computational Efficiency**: Scales as O(L^3) for maximum spherical harmonic degree L
- **Memory Efficiency**: Matrix-free implementation reduces memory requirements by orders of magnitude
- **Accuracy**: Machine precision achievable for well-resolved problems


## Documentation Structure

This documentation is organized into several sections:

- **[Getting Started]({{ '/getting-started/' | relative_url }})** - Installation and basic setup
- **[Tutorials]({{ '/tutorials/' | relative_url }})** - Step-by-step guides to core concepts
- **[Benchmarks]({{ '/benchmark/' | relative_url }})** - Benchmarks of code accuracy and performance


## Benchmarks

Explore our performance benchmarks and validation tests:

{% for benchmark in site.benchmarks limit:3 %}
- [{{ benchmark.title }}]({{ benchmark.url | relative_url }})
{% endfor %}

[View all benchmarks →]({{ '/benchmark/' | relative_url }})

## Citation

If you use gplspec in your research, please cite:

- Myhill, A. D., Maitra, M. A., & Al-Attar, D. (2025). Forward and adjoint calculations of gravitational potential in heterogeneous, aspherical planets. arXiv preprint arXiv:2508.07910.
- Maitra, M., & Al-Attar, D. (2019). A non-perturbative method for gravitational potential calculations within heterogeneous and aspherical planets. Geophysical Journal International, 219(2), 1043-1055.

## License

This project is open source. See the repository for license details.