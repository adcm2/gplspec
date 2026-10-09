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

## Reproducing the reconstructed workflows

Build the maintained API and registered quantity gates with the normal Release configuration:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j2
ctest --test-dir build --output-on-failure
```

CTest runs the Stage 00 two-capture scientific comparison and the Stage 05 1-D/file-backed scientific comparisons. To see those comparisons directly, run `bash tests/run_stage00.sh build/bin/stage00_reference`; run a single Stage 05 gate with `bash tests/run_stage05_scientific_reference.sh stage05_1d build/bin/stage05_operator_reference tests/reference/stage05-1d-original.csv` or the corresponding `stage05_file` executable and frozen CSV. Each runner also prints the unchanged strict comparison as a diagnostic. The Python contract checks for missing/extra keys, non-finite values, and out-of-bound field/action changes run with `python3 tests/test_scientific_reference_gate.py`.

The benchmark executables use paths relative to their current working directory and expect output directories to exist. Run them in a disposable tree, linking the tracked inputs for the readers and preparing every output folder first. The symlink is a path convenience, not a filesystem write-protection mechanism; the programs read model inputs and write only below `work/`. For example:

```sh
repo=$PWD
run=$(mktemp -d)
ln -s "$repo/modeldata" "$run/modeldata"
for n in 1 2 3 4 5 6; do mkdir -p "$run/work/Bench$n"; done
for path in Bench6/lm Bench7/exact Bench7/perturbed Bench7/unperturbed \
            Bench7/zeroradius Bench8/exact Bench8/perturbed \
            Bench8/unperturbed Bench9/exact Bench9/perturbed \
            Bench9/unperturbed Bench10/exact Bench10/perturbed \
            Bench10/unperturbed Bench10/perturbed2; do
  mkdir -p "$run/work/$path"
done
for path in Cartesian Spherical Rotated Perturbed Sensitivity; do
  mkdir -p "$run/work/Phobos/$path"
done
mkdir -p "$run/work/Ziheng"
(cd "$run" && "$repo/build/bin/clean_bench_3")
(cd "$run" && printf '0.001\n' | "$repo/build/bin/clean_bench_9")
(cd "$run" && printf '0.001\n' | "$repo/build/bin/clean_bench_10")
```

The table above gives each program's source parameters and numerical settings; preserve those settings when selecting a different benchmark. `clean_bench_5` reads `lmax` and `maxstep`; `clean_bench_7` through `_10` read the mapping amplitude `h`. `phobos_gravity_perturbation` reads `lsmall`. Supply the intended values on stdin explicitly and record them with the run. The other listed programs use their tracked inputs and source settings without prompts. `OutputZiheng` consumes both `modeldata/TABLEA1.DAT` and `modeldata/prem.200.no`. Existing plot scripts can consume the current output layout from `work/`; their plotting code and input names are retained.

These commands reproduce supported source workflows. The mapping from those runs to particular paper figures or ambiguous Phobos `lsmall` selections remains unverified unless separately documented; a successful run is not itself a paper-level claim.

## What validation establishes

### Current Phase 3 validation

The fresh disconnected Release configure and full build completed successfully; all 25 registered CTests passed, and the standalone multi-header/two-translation-unit link check exited 0. The Stage 00 and Stage 05 scientific gates passed their retained-key and quantity criteria. Unchanged strict comparators are still run and their nonzero diagnostic results are preserved separately.

Ten supported original-settings workflows (`clean_bench_1..4`, `_6`, `_9`, `_10`, `phobos_gravity`, `phobos_heterogeneous`, and `OutputZiheng`) exited 0 and produced all 40 expected files. Thirty-eight files met the initially sufficient conservative physical bound. The Phobos physical-radius output also passes its already-declared geometry tolerance, yielding 39 of 40 passing scientific comparisons. The remaining file is the clean_bench_10 classical `perturbed/MatrixSolution.out`; its `461813.485791 m^2/s^2` maximum difference is retained as a baseline Grid-view lifetime defect attribution, not described as parity or accepted by a relaxed tolerance. The other Phobos radius exact-check difference and all strict exact/bytewise diagnostic results remain visible in the evidence record.

Matched public API timings, including the candidate's independent residual calculation, were 14.026254 ms mapped and 342.357401 ms file-backed versus 13.893254 ms and 333.042666 ms for the old API lifecycle. Candidate medians changed by +0.96% and +2.80%; median RSS was 16,792 KiB versus 16,784 KiB mapped and 393,060 KiB versus 416,476 KiB file-backed. Candidate independent residuals were `2.69037115e-12` and `6.39413501e-10`; the solver algorithm and requested settings remain unchanged.

The exact commands, limits, output hashes, comparison scripts, and build/performance manifests are in [Phase 3F validation evidence](docs/core-cleanup/evidence/phase3f/README.md). These runs reproduce the listed source workflows and do not resolve ambiguous paper labels or interactive Phobos `lsmall` figure mapping.

### Historical Phase 2 validation

The following records the earlier Phase 2 evidence and remains historical.

The Stage 2B Release build compiled all configured targets. Its 14 CTests passed, the Stage 00 runner matched three frozen 20,700-record references with zero changed components and byte-identical representative output, and the standalone header-link smoke check exited successfully. The embedded CMake consumer also configured, built, and ran.

These checks do not execute the full `clean_bench_1..10` or Phobos workflows, do not exercise the interactive prompt values as a campaign, and do not validate plots or establish paper-level reproducibility. The bounded tests and saved evidence are linked from [the cleanup campaign record](docs/cleanup/campaign.md) and [Phase 2 evidence](docs/core-cleanup/evidence/phase2b-validation.log). Keep the frozen references and production inputs unchanged.
