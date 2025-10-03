# gplspec: Gravity solver 

`gplspec` is a C++ library for finding the gravitational potential of isolated, aspherical, heterogeneous planets. This documentation provides an overview of the library, its API, and usage examples.

## Getting Started

This section describes how to build and run the project.

### Prerequisites

*   A C++ compiler (e.g., g++) that meets the C++20 standard 
*   CMake
*   Doxygen (for documentation)
*   NetCDF

### Building the Project

1.  Clone the repository
2.  Create a build directory 
3.  Compile the source code 

Building can be performed via:
```bash
git clone https://github.com/adcm2/gplspec.git
cmake -S . -B build
cmake --build build/
```
Be careful to ensure that the compiler used by CMake is C++20 compliant. If you have multiple compilers installed this can be enforced with the option -DCMAKE_CXX_COMPILER=path_to_compiler.

### Generating Documentation

To generate this documentation, run Doxygen from the project root:
```bash
doxygen Doxyfile
```
The output will be in the `docs/` directory.

## Usage
Here is a basic example of how to use the library:

```cpp
#include <iostream>
#include "gplspec/All"

int
main() {
   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;


   // Construct homogeneous sphere model and compute potential
   Density3D testsphere = Density3D::SphericalHomogeneousPlanet(
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
   auto potential = FindGravitationalPotential(testsphere);

   return 0;
}
```