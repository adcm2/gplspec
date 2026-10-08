---
title: "Home"
---

# gplspec

`gplspec` is a header-only C++ library for gravitational-potential calculations in heterogeneous, aspherical planetary models. It combines particle relabelling with radial spectral elements and spherical-harmonic representation. The implementation and current user workflow are documented in this repository; the pages do not assert accuracy or performance beyond the listed tests and evidence.

## Documentation

- [Getting started](getting-started.md): requirements, build, regression checks, and CMake consumer setup.
- [Supported API overview](api.md): umbrella headers and principal model/gravity interfaces.
- [Scientific workflows and validation](benchmark.md): `clean_bench`, Phobos, files, parameters, output paths, plot consumers, and evidence limits.
- [Tutorials](tutorials.md): model and method explanations.

## Papers

- Myhill, A. D., Maitra, M. A., & Al-Attar, D. (2025). *Forward and adjoint calculations of gravitational potential in heterogeneous, aspherical planets*. [arXiv:2508.07910](https://arxiv.org/abs/2508.07910).
- Maitra, M., & Al-Attar, D. (2019). *A non-perturbative method for gravitational potential calculations within heterogeneous and aspherical planets*. [Geophysical Journal International, 219(2), 1043–1055](https://academic.oup.com/gji/article/219/2/1043/5541065).

## License

See the repository license file.
