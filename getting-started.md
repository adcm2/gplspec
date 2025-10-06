---
layout: default
title: "Getting Started"
permalink: /getting-started/
---

# Getting Started with gplspec

This guide covers the basic steps to compile and use the `gplspec` library.

## Prerequisites

- A modern C++ compiler (supporting C++20 or later)
- CMake (version 3.10 or later)
- Doxygen (optional, for building documentation)
- Graphviz (optional, for generating diagrams)

## Building the Project

1. **Clone the repository:**
   ```bash
   git clone https://github.com/adcm2/gplspec.git
   cd gplspec
   ```

2. **Configure and compile with CMake:**
   ```bash
   cmake -S . -B build
   cmake --build build/
   ```
   
The executables for the examples will be located in the `build/examples` directory.

## Your First Program

Create a simple program to test the installation:

```cpp
#include <iostream>
#include <gplspec/All>

int main() {
    std::cout << "gplspec library loaded successfully!" << std::endl;
    return 0;
}
```

## Next Steps

- Read the [Core Concepts Tutorial]({{ '/tutorials/core-concepts/' | relative_url }})
- Explore the [Examples]({{ '/examples/' | relative_url }})
- Browse the [API Reference]({{ '/api/' | relative_url }})