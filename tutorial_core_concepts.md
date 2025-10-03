/** @page tutorial_core_concepts_page Core Concepts Tutorial */

# Tutorial: Core Concepts

This tutorial introduces the fundamental classes and concepts in `gplspec`.

## The Density3D Class

The `gplspec::Density3D` class is the central component for representing a planet's density distribution. It is built upon a spectral-element discretization in radius and spherical harmonic expansions for angular variations.

You can initialize a `Density3D` model from a 1D profile file, such as PREM, like this:
```cpp
Density3D testprem = Density3D::OneDimensionalPlanetFromFile(
    pathtoprem, npoly, lmax, maxstep, ballrad);
```

## Gravitational Potential Solvers

Once you have a density model, you can compute the corresponding gravitational potential. The library provides multiple methods.

### The Spectral-Element Method
This is the primary method, which solves Poisson's equation using the same spectral-element basis as the density model.
```cpp
auto stdvec_potsol = FindGravitationalPotential(testprem);
```

### The Spherical Integration Method
For comparison and validation, a direct integration method is also available.
```cpp
auto vec_integral_potential = GravitationalSphericalIntegral(testprem);
```

See the @ref clean_bench_2.cpp example for a full benchmark comparison.

