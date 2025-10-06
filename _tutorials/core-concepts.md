---
title: "Core Concepts"
---

# Core Concepts Tutorial

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

## Gravitational Potential Calculation

Once you have a density model, you can compute the gravitational potential using different methods.

### Method 1: Spectral-Element Solver

This is the primary method, optimized for high accuracy and performance:

```cpp
#include <gplspec/Timer>
using namespace Gravity_Tools;

Timer timer;
timer.start();

// Compute potential using spectral-element method
auto potential_spectral = FindGravitationalPotential(earth_model);

timer.stop("Spectral-element method");
```

### Method 2: Spherical Integration

For validation and comparison, a direct integration method is available:

```cpp
timer.start();

// Compute potential using spherical integration
auto potential_integration = GravitationalSphericalIntegral(earth_model);

timer.stop("Spherical integration method");
```

## Working with Results

The potential solvers return vectors containing the gravitational potential at each computational point.

### Output and Visualization

You can save results for analysis:

```cpp
std::string output_directory = "./results";

// Output spectral-element results
earth_model.ReferentialOutputAtElement(output_directory, potential_spectral);

// Output integration results for comparison
earth_model.OutputAtElement(output_directory, potential_integration);
```

## Performance Considerations

### Choosing Parameters

The computational cost scales with your parameter choices:

- **Higher `npoly`**: More accurate but slower
- **Smaller `maxstep`**: More radial resolution but more elements
- **Higher `lmax`**: More lateral detail but increased memory usage

### Typical Parameter Ranges

For different applications:

```cpp
// Quick testing
int npoly = 3; double maxstep = 0.02; int lmax = 1;

// Standard accuracy
int npoly = 5; double maxstep = 0.01; int lmax = 2;

// High precision
int npoly = 7; double maxstep = 0.005; int lmax = 4;
```

## Complete Example

Here's a complete example putting these concepts together:

```cpp
#include <iostream>
#include <gplspec/All>
#include <gplspec/Timer>

int main() {
    using namespace GeneralEarthModels;
    using namespace Gravity_Tools;
    
    // Model parameters
    double maxstep = 0.01;
    double ballrad = 1.2;
    int npoly = 5;
    int lmax = 2;
    std::string prem_file = "modeldata/prem.200";
    
    Timer timer;
    
    std::cout << "Loading planetary model..." << std::endl;
    
    // Create density model
    Density3D model = Density3D::OneDimensionalPlanetFromFile(
        prem_file, npoly, lmax, maxstep, ballrad);
    
    std::cout << "Computing gravitational potential..." << std::endl;
    
    // Compute potential with both methods
    timer.start();
    auto potential_spectral = FindGravitationalPotential(model);
    timer.stop("Spectral-element method");
    
    timer.start();
    auto potential_integration = GravitationalSphericalIntegral(model);
    timer.stop("Integration method");
    
    // Save results
    std::string output_dir = "./tutorial_output";
    model.ReferentialOutputAtElement(output_dir, potential_spectral);
    model.OutputAtElement(output_dir, potential_integration);
    
    std::cout << "Results saved to " << output_dir << std::endl;
    
    return 0;
}
```

## Next Steps

Now that you understand the core concepts:

1. **Explore the [Benchmarks]({{ '/benchmark/' | relative_url }})** to see these concepts in action
2. **Learn about advanced features** in the specialized tutorials
3. **Consult the [API Reference]({{ '/api/' | relative_url }})** for detailed function documentation

The next tutorial covers [Advanced Modeling Techniques]({{ '/tutorials/advanced-modeling/' | relative_url }}) for more complex scenarios.