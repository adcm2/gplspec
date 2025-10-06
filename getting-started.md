---
title: "Getting Started"
permalink: /getting-started/
---

# Getting started with gplspec
gplspec is a header only library, consequently it does not require installation. Indeed, one only needs the header files to be able to include it. It is, however, dependent upon other libraries that we have developed, and consequently it is set up so that it can be included via CMake for ease. Although there is no requirement to use CMake the rest of this guide will describe how to include via CMake.

## Prerequisites

Before including gplspec via CMake, ensure you have:

- **C++ Compiler**: C++20 support
- **CMake**: Version 3.15 or later
- **Git**: For cloning the repository

### Optional Dependencies

- **Doxygen**: For generating API documentation
- **Graphviz**: For generating dependency diagrams in documentation

## Testing gplspec
If you wish to test gplspec before including it in larger projects you can clone the repository and run the examples. Steps to do this using CMake are given below.

### 1. Cloning and building from GitHub

```bash
git clone https://github.com/adcm2/gplspec.git
cd gplspec
```

### 2. Build the Project

Create a build directory and configure with CMake, for example:

```bash
cmake -S . -B build
```

Compile the project:

```bash
cmake --build build/
```

### 3. Run the Examples

Test your installation by running one of the included examples:

```bash
./examples/clean_bench_1
```

### Project Structure

After building, the gplspec directory will look like this:

```
gplspec/
├── src/                 # Source code
├── examples/            # Example programs
├── build/               # Build artifacts
├── CMakeLists.txt       # CMake configuration
└── README.md
```

## Including in larger projects

### Fetching the content
To include it in a program via CMake it is recommended that one fetches the repository from Git, ie include the following in your CMakeLists.txt file:
```cmake 
include(FetchContent)
FetchContent_Declare(
  gplspec
  GIT_REPOSITORY https://github.com/adcm2/gplspec.git
  GIT_TAG main
)
FetchContent_MakeAvailable(gplspec)
```
### A simple program
Once you have included the library and appropriately linked it to your executable you can include different parts or all of the library. For example to include all of the library:

```cpp
// file: hello_gplspec.cpp
#include <iostream>
#include <gplspec/All>

int main() {
    std::cout << "Hello from gplspec!" << std::endl;
    
    // Create a simple test to verify the library is working
    using namespace GeneralEarthModels;
    
    std::cout << "gplspec library loaded successfully!" << std::endl;
    return 0;
}
```



## Build Options

### Debug Build

For development and debugging:

```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
cmake --build build/
```

### Release Build

For optimized performance:

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build build/
```

## Troubleshooting

### Common Issues


**Compilation errors:**
- Ensure you're using a C++20 compatible compiler
- Check that all dependencies are properly installed

**Runtime errors:**
- Verify that input data files (like PREM models) are in the correct location
- Check file permissions for output directories

## Next Steps

Now that you have gplspec installed:

1. **Learn the basics**: Read the [Core Concepts Tutorial]({{ '/tutorials/core-concepts/' | relative_url }})
2. **See it in action**: Explore the [Benchmarks]({{ '/benchmark/' | relative_url }})
3. **Deep dive**: Check the [API Reference]({{ '/api/' | relative_url }})

## Getting Help

If you encounter issues:

- Review the [API documentation]({{ '/api/' | relative_url }})
- Open an issue on the GitHub repository