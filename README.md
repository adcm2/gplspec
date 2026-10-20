# gplspec

`gplspec` is a header-only C++ library for gravitational-potential calculations in aspherical bodies using particle relabelling and radial spectral elements. Its supported CMake target is `gplspec`; it supplies the include paths, pinned header dependencies, and C++23 requirement to linked consumers.

Current documentation: [GPLSpec guides](https://adcm2.github.io/SpecRelab/gplspec/). The legacy GPLSpec Pages address remains a compatibility entry point.

- [Build and use the library](getting-started.md)
- [Supported headers and API overview](api.md)
- [clean_bench and Phobos workflows, outputs, and validation limits](benchmark.md)
- [Scientific tutorials](tutorials.md)

The cleanup campaign records its behavior-preserving baseline and evidence in [docs/cleanup/campaign.md](docs/cleanup/campaign.md). Content under `experimental/` is excluded from production builds and has not been numerically validated. The repository contains no install/export package configuration; consumers currently add the source tree as a CMake subdirectory.
