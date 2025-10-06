---
title: "Gravitational Potential Benchmark"
layout: default
---

# Gravitational Potential Benchmark

This example compares two methods for calculating gravitational potential and measures their performance.

## Overview

The benchmark demonstrates:
- Loading a planetary model from PREM data
- Computing potential with spectral-element method
- Computing potential with spherical integration
- Performance comparison

## Source Code

```cpp
#include <gplspec/All>
#include <gplspec/Timer>

int main() {
    using namespace GeneralEarthModels;
    using namespace Gravity_Tools;
    
    // Setup parameters
    double maxstep = 0.01;
    double ballrad = 1.2;
    int npoly = 5;
    int lmax = 2;
    std::string pathtoprem = "modeldata/prem.200";
    
    Timer timer;
    
    // Load PREM model
    Density3D testprem = Density3D::OneDimensionalPlanetFromFile(
        pathtoprem, npoly, lmax, maxstep, ballrad);
    
    // Method 1: Spectral-element
    timer.start();
    auto potential_3d = FindGravitationalPotential(testprem);
    timer.stop("3D Spectral Method");
    
    // Method 2: Spherical integration
    timer.start();
    auto potential_integral = GravitationalSphericalIntegral(testprem);
    timer.stop("Spherical Integration");
    
    return 0;
}
```

## Expected Output

The program will display timing information comparing both methods, typically showing the spectral-element method is faster for complex models.