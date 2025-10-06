---
title: "Core Concepts"
layout: default
---

# Core Concepts Tutorial

This tutorial introduces the fundamental classes and concepts in `gplspec`.

## The Density3D Class

The `Density3D` class is the central component for representing a planet's density distribution. It uses spectral-element discretization in radius and spherical harmonic expansions for angular variations.

### Creating a Density Model

You can initialize a `Density3D` model from a 1D profile file like PREM:

```cpp
Density3D testprem = Density3D::OneDimensionalPlanetFromFile(
    "prem.200", npoly, lmax, maxstep, ballrad);
```

## Gravitational Potential Solvers

Once you have a density model, you can compute the gravitational potential using different methods:

### Spectral-Element Method
```cpp
auto potential = FindGravitationalPotential(testprem);
```

### Spherical Integration Method
```cpp
auto potential = GravitationalSphericalIntegral(testprem);
```

## Performance Considerations

The spectral-element method is typically faster for complex 3D models, while the spherical integration method is useful for validation and simple models.