---
title: "Getting Started"
permalink: /getting-started/
---

# Getting Started with gplspec

This guide will help you install, build, and run your first gplspec program.

## Prerequisites

Before building gplspec, ensure you have:

- **C++ Compiler**: C++20 support
- **CMake**: Version 3.15 or later
- **Git**: For cloning the repository

### Optional Dependencies

- **Doxygen**: For generating API documentation
- **Graphviz**: For generating dependency diagrams in documentation

## Installation

### 1. Clone the Repository

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

## Your First Program

Create a simple program to verify everything is working:

### Step 1: Create a new file

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

### Step 2: Compile and run

```bash
g++ -std=c++20 -I/path/to/gplspec/src hello_gplspec.cpp -o hello_gplspec
./hello_gplspec
```

## Project Structure

After building, your project directory will look like this:

```
gplspec/
├── src/                 # Source code
├── examples/            # Example programs
├── build/               # Build artifacts
├── CMakeLists.txt       # CMake configuration
└── README.md
```

## Build Options

### Debug Build

For development and debugging:

```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```

### Release Build

For optimized performance:

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

### Building Documentation

If you have Doxygen installed:

```bash
cmake -DBUILD_DOCS=ON ..
make docs
```

## Troubleshooting

### Common Issues

**CMake can't find dependencies:**
```bash
# Make sure you have the required packages
sudo apt-get install build-essential cmake git
```

**Compilation errors:**
- Ensure you're using a C++20 compatible compiler
- Check that all dependencies are properly installed

**Runtime errors:**
- Verify that input data files (like PREM models) are in the correct location
- Check file permissions for output directories

## Next Steps

Now that you have gplspec installed:

1. **Learn the basics**: Read the [Core Concepts Tutorial]({{ '/tutorials/core-concepts/' | relative_url }})
2. **See it in action**: Explore the [Examples]({{ '/examples/' | relative_url }})
3. **Deep dive**: Check the [API Reference]({{ '/api/' | relative_url }})

## Getting Help

If you encounter issues:

- Check the [Examples]({{ '/examples/' | relative_url }}) for similar use cases
- Review the [API documentation]({{ '/api/' | relative_url }})
- Open an issue on the GitHub repository