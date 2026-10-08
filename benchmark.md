---
title: "Scientific workflows and validation"
layout: default
permalink: /benchmark/
---

# Scientific workflows and validation

This page maps the tracked `clean_bench` and Phobos programs to their inputs, effective parameters, outputs, and plotting consumers. It describes what a user can run; it does not claim that every workflow was executed or that every model is scientifically validated. All relative paths below are resolved from the program's current working directory. Programs do not create output directories. Run each in a disposable directory with the required `work/...` directories prepared. Do not run these workflows in a checkout containing saved outputs you need to preserve.

Build the executables as described in [Getting started](getting-started.md). Executables are in `build/bin/`. To plot a benchmark, run its script from the temporary `work/` directory containing the named `BenchN` outputs; scripts in the source tree are under `work/`.

## clean_bench programs

All programs use fixed source parameters except where an stdin prompt is listed. “No file input” means the source's declared path is unused or the model constructor is commented out. Model inputs are paths relative to the process working directory.

| Program | Inputs and effective setup | Output paths and consumers |
|---|---|---|
| `clean_bench_1` | No files or prompts. One homogeneous-material sphere discretized using max step 0.1; length norm 6,371,000, time norm 3,600, mass norm 5.972e24; max step 0.1, ball radius 1.2. Compares pseudospectral solution, spherical integral, and homogeneous-sphere analytic result. | `work/Bench1/{MatrixSolution,IntegralSolution,ExactSolution}.out`; `work/CleanBench1Plot.py` reads these. |
| `clean_bench_2` | Reads `modeldata/prem.200`; max step 0.1, ball radius 1.2, polynomial order 5, harmonic degree 2. Compares pseudospectral solution and spherical integral. | `work/Bench2/{MatrixSolution,IntegralSolution}.out`; `work/CleanBench2Plot.py`. |
| `clean_bench_3` | No effective file input; homogeneous unit-radius/unit-density sphere, mapped radius `h*r`, with `h=0.2`; polynomial order 5, degree 2, max step 0.1, ball radius 1.4. Compares mapped solution with the equivalent sphere (radius 1.2, density `1/(1.2^3)`); tolerance 1e-12. | `work/Bench3/MatrixSolution.out` and `work/Bench3/ExactSolution.out`; `work/CleanBench3Plot.py`. |
| `clean_bench_4` | No effective file input; unit sphere, map `h*r*(1-r/R)`, `h=0.02`; polynomial order 5, degree 2, max step 0.1, ball radius 1.2; tolerance 1e-12. | `work/Bench4/MatrixSolution.out` and `work/Bench4/ExactSolution.out`; `work/CleanBench4Plot.py`. |
| `clean_bench_5` | No files. Unit sphere, map `h*sin(theta)*cos(theta)*sin(phi)*r*(1-r/R)`, `h=0.2`; polynomial order 8, ball radius 1.4. Prompts for `lmax2` (initial value 100) and `maxstep` (initial value 0.05); declared degree 2 is not the effective prompt default. Tolerance 1e-12. | `work/Bench5/MatrixSolution.out` and `work/Bench5/ExactSolution.out`; `work/CleanBench5Plot.py`. |
| `clean_bench_6` | Reads `modeldata/prem.200` and `modeldata/S40RTS_dvs.nc`; max step 0.05, ball radius 1.2, polynomial order 12, degree 20. Pseudospectral solve tolerance 1e-14, spherical integral, and sensitivity for `(l,m)=(0,0),(1,0),(2,2),(3,3),(4,0)` with multipliers 1,3,5,7,9. | `work/Bench6/` matrix and integral outputs; `work/Bench6/lm/` physical/referential and sensitivity outputs. `CleanBench6Plot.py`, `CleanBench6Plot1.py`, and `CleanBench6Plot_SK.py` consume the corresponding products. |
| `clean_bench_7` | No effective file input; homogeneous normalized sphere, zero tomography, max step 0.1, ball radius 1.5, polynomial order 5, degree 20. Prompts for `h`; map `h*r`. Compares mapped and baseline solutions and sweeps 100 values from 0 to below 0.2 for a zero-radius product. Tolerance 1e-12. | `work/Bench7/exact/MatrixSolution.out`, `perturbed/MatrixSolution.out`, `unperturbed/MatrixSolution.out`, and `zeroradius/MatrixSolution.out` beneath `work/Bench7/`. `CleanBench7Plot.py` and `CleanBench7Plot1.py`. |
| `clean_bench_8` | No effective file input; same base model/resolution as 7. Prompts for `h`; maps `h*cos(theta)*r*(1-r/R)` and a 2h variant; tolerance 1e-12. | `work/Bench8/{exact,perturbed,unperturbed}/MatrixSolution.out`; `work/CleanBench8Plot.py`. |
| `clean_bench_9` | No effective file input; same base model/resolution as 7. Prompts for `h`; map `h*r`; compares direct/classical perturbation; tolerance 1e-12. | `work/Bench9/{exact,perturbed,unperturbed}/MatrixSolution.out`; `work/CleanBench9Plot.py`. |
| `clean_bench_10` | No effective file input; same base model, degree 30. Prompts for `h`; map `h*cos(theta)*r`; compares exact, classical, and mapping perturbations, each at 1e-12. | `work/Bench10/{exact,perturbed,unperturbed,perturbed2}/MatrixSolution.out`; `work/CleanBench10Plot.py`. |

The exact output spelling and the plotting scripts are part of the workflow contract. Some scripts read additional named or alternative products; inspect the script before running a partial producer workflow. `work/` currently contains tracked plot scripts and saved products/figures; this map does not authorize replacing them.

## Phobos workflows

The Phobos model constructor and `OutputZiheng` read the tracked `modeldata/TABLEA1.DAT` through the typed `GPLSpec::ModelInput::ReadPhobosHarmonics` input function. Common source parameters are physical density 1860, length norm 10,000, used density 1342, mass norm `lengthnorm^3 * useddensity`, time norm 3600, scaled density `physical_density / massnorm * lengthnorm^3`, degree 100, polynomial order 5, max step 0.1, and ball radius 1.5. They write below `work/Phobos/` or `work/Ziheng/`; create those directories inside a disposable run tree first.

| Program | Differences and outputs |
|---|---|
| `phobos_gravity` | Solves with tolerances 1e-12 and 1e-3 and computes sensitivity for degree 2, order 0. Writes Cartesian products under `work/Phobos/Cartesian/`, spherical products and rotated files under `work/Phobos/Spherical/`, and sensitivity products under `work/Phobos/Sensitivity/`. |
| `phobos_heterogeneous` | Uses heterogeneous density (`nhomo=true`), solves at 1e-5 and 1e-3, computes degree-2/order-0 sensitivity and a power spectrum. Writes the common coordinate products and `work/Phobos/Spherical/PowerSpectrum.out`. |
| `phobos_gravity_perturbation` | Prompts for `lsmall` (initial value 0); computes full and low-degree models at 1e-6, a mapping perturbation truncated from `lsmall+1` through 100, and the perturbed solve at 1e-6. Writes Cartesian, spherical, rotated, and perturbed products under `work/Phobos/`. |
| `OutputZiheng` | Reads `modeldata/TABLEA1.DAT` for Phobos and separately reads `modeldata/prem.200.no` through `OneDimensionalPlanetFromFile(path, 5, 2, 0.001, 1.2)` to produce a PREM spherical-integral comparison. Uses the common Phobos model scales, solves at 1e-12, and applies the source coefficient transforms. Writes `work/Ziheng/PhobosRadius.out`, `work/Ziheng/PHOBOS_POTENTIAL_OUTER.out`, and `work/Ziheng/PREM200.NO.INTEGRAL.out`. The other coordinate-output block in this source is commented out. |

The tracked Phobos plot/analysis scripts read these relative paths from their working directory:

| Consumer | Products read under `work/` |
|---|---|
| `PhobosFullSlice.py` | `Phobos/Cartesian/MatrixSolution.out`; `Phobos/Spherical/MatrixSolutionReferentialRotated.out`, `MatrixSolutionPhysicalRotated.out`, `ReferentialDensitySolutionRotated.out`, `PhysicalDensitySolutionRotated.out`. |
| `PhobosPerturbation.py` | `Phobos/Cartesian/MatrixSolution.out`, `Phobos/Spherical/MatrixSolution.out`, `Phobos/Rotated/MatrixSolution.out`, `Phobos/Perturbed/MatrixSolution.out`. |
| `PhobosSensitivityKernel.py` | Cartesian, spherical, and rotated `MatrixSolution.out`; `Phobos/Sensitivity/SensitivityKernel.out`. |
| `PhobosSensitivityKernel2.py` | Cartesian `MatrixSolution.out`; `Phobos/Sensitivity/SensitivityKernel.out` and `SensitivityKernelSlice.out`. |
| `PhobosSliceDensity.py` | `Phobos/Cartesian/MatrixSolution.out` and `Phobos/Spherical/MatrixSolutionReferentialRotated.out`. |
| `PhobosSlicePlot.py` | `Phobos/Cartesian/MatrixSolution.out`, `Phobos/Sensitivity/SensitivityKernel.out`, and `Phobos/Spherical/MatrixSolutionSlice.out`. |
| `PhobosSlicePlot2.py` | `Phobos/Cartesian/MatrixSolution.out`, `Phobos/Sensitivity/SensitivityKernel.out`, `Phobos/Spherical/MatrixSolutionSlice.out`, and `Phobos/Spherical/MatrixSolutionPhysical.out`. The source assigns `path6` several times; the final active value is the physical-output file. |
| `PhobosSurfacePotential.py` | Cartesian and spherical `MatrixSolution.out`. |
| `PhobosSurfacePotentialRotated.py` | Cartesian, spherical, and rotated `MatrixSolution.out`. |
| `PhobosSurfaceTopography.py` | Cartesian and spherical `MatrixSolution.out`. |
| `PowerSpectrum.py` | `Phobos/Spherical/PowerSpectrum.out`. |

These consumers cover more than one producer/output variant; inspect a script and select a producer that actually writes every file it reads. For example, the rotated `Phobos/Rotated/MatrixSolution.out` is written by the perturbation program, while `phobos_gravity` writes its rotated coordinate products into `Phobos/Spherical/`. Preserve each producer-consumer pair, model input, relative path, and saved product.

## What validation establishes

The Stage 2B Release build compiled all configured targets. Its 14 CTests passed, the Stage 00 runner matched three frozen 20,700-record references with zero changed components and byte-identical representative output, and the standalone header-link smoke check exited successfully. The embedded CMake consumer also configured, built, and ran.

These checks do not execute the full `clean_bench_1..10` or Phobos workflows, do not exercise the interactive prompt values as a campaign, and do not validate plots or establish paper-level reproducibility. The bounded tests and saved evidence are linked from [the cleanup campaign record](docs/cleanup/campaign.md) and [Phase 2 evidence](docs/core-cleanup/evidence/phase2b-validation.log). Keep the frozen references and production inputs unchanged.
