---
title: "Homogeneous sphere"
---

# Homogeneous sphere

This tutorial introduces the fundamental concepts and classes in gplspec, helping you understand how to effectively use the library for gravitational modeling.

## Overview

gplspec is built around several key concepts:

1. **Density Models** - Representing planetary density distributions
2. **Spectral Elements** - High-order numerical methods for accuracy
3. **Spherical Harmonics** - Global field representations
4. **Potential Solvers** - Computing gravitational fields from density

## The Density3D Class

The `Density3D` class is the cornerstone of gplspec, representing a three-dimensional density distribution of a planetary body.

### Creating a Density Model

The most common way to create a density model is from a 1D radial profile:

```cpp
#include <gplspec/All>
using namespace GeneralEarthModels;

// Parameters for the model
double maxstep = 0.01;    // Maximum radial step size
double ballrad = 1.2;     // Radius of computational domain
int npoly = 5;            // Polynomial degree for spectral elements
int lmax = 2;             // Maximum spherical harmonic degree
std::string modelfile = "modeldata/prem.200";

// Create the 3D model from 1D profile
Density3D earth_model = Density3D::OneDimensionalPlanetFromFile(
    modelfile, npoly, lmax, maxstep, ballrad);
```

### Understanding the Parameters

- **`maxstep`**: Controls the radial resolution of the spectral elements
- **`ballrad`**: Defines the size of the computational domain (normalized units)
- **`npoly`**: Higher values give more accuracy but increase computational cost
- **`lmax`**: Maximum degree of spherical harmonic expansion for lateral variations