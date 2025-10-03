/** @page getting_started_page Getting Started */

# Getting Started with gplspec

This guide covers the basic steps to compile and use the `gplspec` library.

## Prerequisites

- A modern C++ compiler (supporting C++20 or later)
- CMake (version 3.10 or later)
- Doxygen (optional, for building documentation)
- Graphviz (optional, for generating diagrams)

## Building the Project

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/adcm2/gplspec.git
    cd gplspec
    ```

2.  **Configure and compile with CMake:**
    ```bash
    cmake -S . -B build
    cmake --build build/
    ```
    The executables for the examples will be located in the `build/examples` directory.

## Generating Documentation

To generate the documentation yourself, run Doxygen from the project's root directory:
```bash
doxygen Doxyfile
```
The output will be generated in the `docs/` folder. Open `docs/index.html` to view it.
