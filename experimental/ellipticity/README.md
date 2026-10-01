# Experimental ellipticity comparison

Historical experiment copied unchanged from Gravitational_Field at commit
`50a0264427d6de7f1fef4ef340932b3c2f361241` on 2026-10-01.
File hashes are recorded in `provenance.json`.

## Contents and possible use

- `examples/ellipticity_test.cpp`: compares a finite-difference ellipticity solution, spectral-element formulations of Clairaut's equation, and ellipticity recovered from a potential-based spectral-element calculation.
- `work/ellipticityplot.py`: compares the profiles and relative differences from `gravitycheck.out` and `ellipticity2.out`.
- `Gravitational_Field/TestEllipticity` and `Gravitational_Field/src/ellipticitytools.h`: historical supporting headers. The `FDEllipticitySolver` class is only a placeholder; the calculation is implemented in the example.

This could provide an independent cross-check for future equilibrium-figure or ellipticity work. Its presence does not establish numerical correctness.

## Status and future integration

This is a source-preservation copy, not an enabled GPLSpec example. It has not been compiled or run as part of this transfer. No CMake targets or scientific outputs were added, and no numerical expressions were changed.

The original directory structure and legacy `Gravitational_Field` includes are deliberately retained. Before compiling, adapt the includes to the current `gplspec` headers, review the experimental header arrangement and current dependency APIs, and add an opt-in build target. This directory alone is not a complete copy of the old library.

The example expects `modeldata/prem.200` relative to its working directory and writes `work/gravitycheck.out` and `work/ellipticity2.out`. GPLSpec already contains an identical copy of that input model, so it was not duplicated here. For any future run, use a separate output directory and adjust the plotting script's input paths accordingly; do not overwrite existing research results.

Validate boundary/interface treatment, normalization, convergence and agreement between the formulations before treating this as a supported method.
