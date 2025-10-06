---
title: "Gravitational Potential Benchmark"
---

# Gravitational Potential Benchmark

This example demonstrates how to benchmark different methods for computing gravitational potential and compare their performance and accuracy.

## Overview

This benchmark compares two approaches:

1. **Spectral-Element Method**: High-order finite element approach
2. **Spherical Integration**: Direct numerical integration

The example uses the PREM (Preliminary Reference Earth Model) to represent Earth's density structure.

## Complete Source Code

```cpp
#include <GaussQuad/All>
#include <gplspec/All>
#include <gplspec/Test>
#include <gplspec/Timer>
#include <PlanetaryModel/All>
#include <TomographyModels/All>
#include <filesystem>
#include <iomanip>
#include <iostream>

int main() {
    using namespace GeneralEarthModels;
    using namespace Gravity_Tools;

    // Step 1: Setup Model Parameters
    double maxstep = 0.01;              // Maximum radial step size
    double ballrad = 1.2;               // Computational domain radius
    int npoly = 5;                      // Polynomial degree
    int lmax = 2;                       // Max spherical harmonic degree
    std::string pathtoprem = "modeldata/prem.200";  // PREM data file
    
    Timer timer;

    std::cout << "=== Gravitational Potential Benchmark ===" << std::endl;
    std::cout << "Model parameters:" << std::endl;
    std::cout << "  Max step size: " << maxstep << std::endl;
    std::cout << "  Ball radius: " << ballrad << std::endl;
    std::cout << "  Polynomial degree: " << npoly << std::endl;
    std::cout << "  Max l: " << lmax << std::endl;
    std::cout << std::endl;

    // Step 2: Load PREM Model and Create 3D Density
    std::cout << "Loading PREM model..." << std::endl;
    timer.start();
    
    Density3D testprem = Density3D::OneDimensionalPlanetFromFile(
        pathtoprem, npoly, lmax, maxstep, ballrad);
    
    timer.stop("Model loading");
    std::cout << std::endl;

    // Step 3: Calculate Potential with Spectral-Element Method
    std::cout << "Computing potential with spectral-element method..." << std::endl;
    timer.start();
    
    auto stdvec_potsol = FindGravitationalPotential(testprem);
    
    timer.stop("Spectral-element method");
    std::cout << std::endl;

    // Step 4: Calculate Potential with Spherical Integration
    std::cout << "Computing potential with spherical integration..." << std::endl;
    timer.start();
    
    auto vec_integral_potential = GravitationalSphericalIntegral(testprem);
    
    timer.stop("Spherical integration method");
    std::cout << std::endl;

    // Step 5: Output Results for Analysis
    std::string output_path = "./work/Bench2";
    std::cout << "Saving results to " << output_path << "..." << std::endl;
    
    // Create output directory if it doesn't exist
    std::filesystem::create_directories(output_path);
    
    // Save spectral-element results
    testprem.ReferentialOutputAtElement(output_path, stdvec_potsol);
    
    // Save integration results
    testprem.OutputAtElement(output_path, vec_integral_potential);
    
    std::cout << "Benchmark complete!" << std::endl;

    return 0;
}
```

## Understanding the Code

### Step 1: Parameter Setup

The benchmark begins by defining key parameters that control the accuracy and computational cost:

```cpp
double maxstep = 0.01;    // Smaller values = higher radial resolution
double ballrad = 1.2;     // Computational domain size
int npoly = 5;            // Higher values = more accurate elements
int lmax = 2;             // Spherical harmonic expansion degree
```

### Step 2: Model Creation

The PREM model is loaded and converted to a 3D representation:

```cpp
Density3D testprem = Density3D::OneDimensionalPlanetFromFile(
    pathtoprem, npoly, lmax, maxstep, ballrad);
```

This creates a spectral-element mesh with spherical harmonic basis functions.

### Step 3: Spectral-Element Solution

The primary solver uses the same basis functions as the density representation:

```cpp
auto stdvec_potsol = FindGravitationalPotential(testprem);
```

This method is typically faster and more accurate for complex 3D models.

### Step 4: Integration Method

For validation, a direct spherical integration is performed:

```cpp
auto vec_integral_potential = GravitationalSphericalIntegral(testprem);
```

This method is slower but provides an independent verification of results.

## Expected Output

When you run this benchmark, you should see output similar to:

```
=== Gravitational Potential Benchmark ===
Model parameters:
  Max step size: 0.01
  Ball radius: 1.2
  Polynomial degree: 5
  Max l: 2

Loading PREM model...
Model loading: 0.15 seconds

Computing potential with spectral-element method...
Spectral-element method: 2.34 seconds

Computing potential with spherical integration...
Spherical integration method: 8.71 seconds

Saving results to ./work/Bench2...
Benchmark complete!
```

## Performance Analysis

### Typical Performance Characteristics

- **Spectral-element method**: Usually 3-5x faster than integration
- **Memory usage**: Scales with `npoly³` and number of elements
- **Accuracy**: Both methods should agree to within numerical precision

### Parameter Impact on Performance

| Parameter | Effect on Speed | Effect on Accuracy |
|-----------|-----------------|-------------------|
| `maxstep` ↓ | Slower (more elements) | Higher accuracy |
| `npoly` ↑ | Slower (more DOF per element) | Higher accuracy |
| `lmax` ↑ | Slower (more spherical modes) | More lateral detail |

## Analyzing Results

The output files contain the gravitational potential at each computational point. You can analyze them using visualization tools or compare them numerically:

```cpp
// Example: Simple comparison
double max_difference = 0.0;
for (size_t i = 0; i < stdvec_potsol.size(); ++i) {
    double diff = std::abs(stdvec_potsol[i] - vec_integral_potential[i]);
    max_difference = std::max(max_difference, diff);
}
std::cout << "Maximum difference: " << max_difference << std::endl;
```

## Modifying the Benchmark

### Testing Different Parameters

Try different parameter combinations to see their effect:

```cpp
// High accuracy test
double maxstep = 0.005;
int npoly = 7;
int lmax = 4;

// Quick test
double maxstep = 0.02;
int npoly = 3;
int lmax = 1;
```

### Using Different Models

You can replace PREM with other 1D models or create synthetic models:

```cpp
// Example: Homogeneous sphere
Density3D sphere = Density3D::SphericalHomogeneousPlanet(
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
```

## Next Steps

- Try the [Advanced Modeling Tutorial]({{ '/tutorials/advanced-modeling/' | relative_url }})
- Explore [other examples]({{ '/examples/' | relative_url }})
- Check the [API Reference]({{ '/api/' | relative_url }}) for more functions