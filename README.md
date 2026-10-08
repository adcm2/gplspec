# gplspec

`gplspec` is a header-only C++ library for gravitational-potential calculations in aspherical bodies using particle relabelling. The project uses CMake to fetch its pinned source dependencies and to build examples and regression checks. See [Getting started](getting-started.md) for prerequisites, build instructions, and public-consumer setup.

The cleanup campaign records its accepted behavior-preserving baseline, validation evidence, and known limitations in [docs/cleanup/campaign.md](docs/cleanup/campaign.md). The historical experiments in `experimental/` are excluded from production builds and have not been numerically validated.
