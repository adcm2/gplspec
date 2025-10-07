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

The `Density3D` class is the main class within gplspec, representing a three-dimensional density distribution of a planetary body. There are other classes which can be used for one-dimensional or spherical bodies, but there is little to be gained in most practical scenarios.

### Creating a Density Model

There are multiple ways in which to construct a density model. These are discussed in more detail through the rest of the tutorials. In the below example a one dimensional model is constructed using the data from PREM. 

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
- **`npoly`**: Maximum polynomial degree of the radial spectral elements
- **`lmax`**: Maximum degree of spherical harmonic expansion for lateral variations

## Gravitational Potential Calculation

Once you have a density model, you can compute the gravitational potential using different methods.

### Method 1: Spectral-Element Solver

This is the primary method, optimized for high accuracy and performance:

```cpp
#include <gplspec/All>
using namespace Gravity_Tools;

// Compute potential using spectral-element method
auto potential_spectral = FindGravitationalPotential(earth_model);

```

### Method 2: Spherical Integration

For validation and comparison, a direct integration method is available:

```cpp
// Compute potential using spherical integration
auto potential_integration = GravitationalSphericalIntegral(earth_model);

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

- **Higher `npoly`**: Generally increases accuracy
- **Smaller `maxstep`**: Decreases the maximum element size
- **Higher `lmax`**: More lateral detail but increased memory usage


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

1. **Explore the rest of the [Tutorials]({{ '/tutorials/' | relative_url }})** to see these concepts in action
2. **See the [Benchmarks]({{ '/benchmark/' | relative_url }})** for validation of the model framework
