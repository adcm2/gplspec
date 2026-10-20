---
title: "Scientific workflows and validation"
layout: default
permalink: /benchmark/
---

# Scientific workflows and validation

This page records the scientific inputs and effective settings for the tracked `clean_bench` and Phobos workflows. It describes what users can run; it does not claim that every workflow was executed or that every model is scientifically validated. Relative input paths are resolved from the process working directory. Current producers write schema-version-1 field bundles under `work/`; the shared bundle writer creates output directories as needed. Run workflows in isolated temporary directories so repository data and saved outputs remain untouched.

Build executables as described in [Getting started](getting-started.md). For reproducible producer commands, current bundle paths, and plotting commands using `work/plot_fields.py`, follow the [field-bundle plotting guide](work/plots.md). That guide is authoritative for the current solve-to-plot path.

## clean_bench programs

All programs use fixed source parameters except where an stdin prompt is listed. “No file input” means the source's declared path is unused or the model constructor is commented out. Model inputs are paths relative to the process working directory.

| Program | Inputs and effective setup | Output paths and consumers |
|---|---|---|
| `clean_bench_1` | No files or prompts. One homogeneous-material sphere discretized using max step 0.1; length norm 6,371,000, time norm 3,600, mass norm 5.972e24; max step 0.1, ball radius 1.2. Compares pseudospectral solution, spherical integral, and homogeneous-sphere analytic result. | Bundles: `work/Bench1/{field_bundle,integral_bundle,exact_bundle}`. See [the bundle guide](work/plots.md) for plotting. |
| `clean_bench_2` | Reads `modeldata/prem.200`; max step 0.1, ball radius 1.2, polynomial order 5, harmonic degree 2. Compares pseudospectral solution and spherical integral. | Bundles: `work/Bench2/{field_bundle,integral_bundle}`. See [the bundle guide](work/plots.md) for plotting. |
| `clean_bench_3` | No effective file input; homogeneous unit-radius/unit-density sphere, mapped radius `h*r`, with `h=0.2`; polynomial order 5, degree 2, max step 0.1, ball radius 1.4. Compares mapped solution with the equivalent sphere (radius 1.2, density `1/(1.2^3)`); tolerance 1e-12. | Bundles: `work/Bench3/{field_bundle,exact_bundle}`. See [the bundle guide](work/plots.md) for plotting. |
| `clean_bench_4` | No effective file input; unit sphere, map `h*r*(1-r/R)`, `h=0.02`; polynomial order 5, degree 2, max step 0.1, ball radius 1.2; tolerance 1e-12. | Bundles: `work/Bench4/{field_bundle,exact_bundle}`. See [the bundle guide](work/plots.md) for plotting. |
| `clean_bench_5` | No files. Unit sphere, map `h*sin(theta)*cos(theta)*sin(phi)*r*(1-r/R)`, `h=0.2`; polynomial order 8, ball radius 1.4. Prompts for `lmax2` (initial value 100) and `maxstep` (initial value 0.05); declared degree 2 is not the effective prompt default. Tolerance 1e-12. | Bundles: `work/Bench5/{field_bundle,exact_bundle}`. See [the bundle guide](work/plots.md) for plotting. |
| `clean_bench_6` | Reads `modeldata/prem.200` and `modeldata/S40RTS_dvs.nc`; max step 0.05, ball radius 1.2, polynomial order 12, degree 20. Pseudospectral solve tolerance 1e-14, spherical integral, and sensitivity for `(l,m)=(0,0),(1,0),(2,2),(3,3),(4,0)` with multipliers 1,3,5,7,9. | Bundles: `work/Bench6/field_bundle` and `work/Bench6/lm/{harmonic_bundle,integral_bundle,sensitivity_bundle}`. The retired plot roles are mapped to `work/plot_fields.py` in [the bundle guide](work/plots.md). |
| `clean_bench_7` | No effective file input; homogeneous normalized sphere, zero tomography, max step 0.1, ball radius 1.5, polynomial order 5, degree 20. Prompts for `h`; map `h*r`. Compares mapped and baseline solutions and sweeps 100 values from 0 to below 0.2 for a zero-radius product. Tolerance 1e-12. | Bundles: `work/Bench7/{exact,perturbed,unperturbed}/field_bundle`; the zero-radius sweep remains a diagnostic table at `work/Bench7/zeroradius/MatrixSolution.out`. Plot the field bundles using [the bundle guide](work/plots.md). |
| `clean_bench_8` | No effective file input; same base model/resolution as 7. Prompts for `h`; maps `h*cos(theta)*r*(1-r/R)` and a 2h variant; tolerance 1e-12. | Bundles: `work/Bench8/{exact,perturbed,unperturbed}/field_bundle`. See [the bundle guide](work/plots.md) for plotting. |
| `clean_bench_9` | No effective file input; same base model/resolution as 7. Prompts for `h`; map `h*r`; compares direct/classical perturbation; tolerance 1e-12. | Bundles: `work/Bench9/{exact,perturbed,unperturbed}/field_bundle`. See [the bundle guide](work/plots.md) for plotting. |
| `clean_bench_10` | No effective file input; same base model, degree 30. Prompts for `h`; map `h*cos(theta)*r`; compares exact, classical, and mapping perturbations, each at 1e-12. | Bundles: `work/Bench10/{exact,perturbed,unperturbed,perturbed2}/field_bundle`. See [the bundle guide](work/plots.md) for plotting. |

Current bundle names and plot commands are documented in [work/plots.md](work/plots.md). Retired solver `.out` products and plot-script paths are historical; the retained Bench 7 `work/Bench7/zeroradius/MatrixSolution.out` diagnostic table is an input to `work/plot_fields.py sweep`, as shown in the bundle guide. `work/` also contains protected saved products and figures; use an isolated run directory and do not replace them.

## Phobos workflows

The Phobos model constructor and `OutputZiheng` read the tracked `modeldata/TABLEA1.DAT` through the typed `GPLSpec::ModelInput::ReadPhobosHarmonics` input function. Common source parameters are physical density 1860, length norm 10,000, used density 1342, mass norm `lengthnorm^3 * useddensity`, time norm 3600, scaled density `physical_density / massnorm * lengthnorm^3`, degree 100, polynomial order 5, max step 0.1, and ball radius 1.5. Their bundle products are listed in [work/plots.md](work/plots.md); bundle writers create the output directories. Run them in an isolated temporary tree.

| Program | Differences and outputs |
|---|---|
| `phobos_gravity` | Solves with tolerances 1e-12 and 1e-3 and computes sensitivity for degree 2, order 0. Writes field bundles under `work/Phobos/Spherical/` and `work/Phobos/Sensitivity/`; exact bundle names and plot commands are in [work/plots.md](work/plots.md). |
| `phobos_heterogeneous` | Uses heterogeneous density (`nhomo=true`), solves at 1e-5 and 1e-3, computes degree-2/order-0 sensitivity and a power spectrum. Writes field and degree-power bundles under `work/Phobos/Spherical/` and sensitivity bundles under `work/Phobos/Sensitivity/`; see [work/plots.md](work/plots.md). |
| `phobos_gravity_perturbation` | Prompts for `lsmall` (initial value 0); computes full and low-degree models at 1e-6, a mapping perturbation truncated from `lsmall+1` through 100, and the perturbed solve at 1e-6. Writes field bundles under `work/Phobos/Spherical/`, `work/Phobos/Rotated/`, and `work/Phobos/Perturbed/`; see [work/plots.md](work/plots.md) for supported plot commands. |
| `OutputZiheng` | Reads `modeldata/TABLEA1.DAT` for Phobos and separately reads `modeldata/prem.200.no` through `OneDimensionalPlanetFromFile(path, 5, 2, 0.001, 1.2)` to produce a PREM spherical-integral comparison. Uses the common Phobos model scales, solves at 1e-12, and applies the source coefficient transforms. Writes `work/Ziheng/PhobosRadius.bundle`, `work/Ziheng/PHOBOS_POTENTIAL_OUTER.bundle`, and `work/Ziheng/PREM200.NO.INTEGRAL.bundle`; see [work/plots.md](work/plots.md) for current plots. |

Current plot consumers use the field-bundle reader and `work/plot_fields.py`; the [bundle plotting guide](work/plots.md) maps the retired script roles to bundle quantities and executable commands. It includes isolated Phobos and topography recipes. The retained Bench 7 sweep table is a current `sweep` input; other historical solver `.out` product filenames remain in frozen evidence and are not current consumer inputs.

## Reproducing the reconstructed workflows

Build the maintained API and registered quantity gates with the normal Release configuration:

```sh
python3.12 -m venv build/python-test-env
build/python-test-env/bin/python -m pip install -r tests/requirements.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON \
  -DMY_PROJECT_BUILD_EXAMPLES=ON \
  -DGPLSPEC_BUILD_BASELINE_HARNESS=ON \
  -DPYTHON_EXECUTABLE="$PWD/build/python-test-env/bin/python"
cmake --build build -j2
ctest --test-dir build --output-on-failure --no-tests=error
```

The Python environment supplies the NumPy and Matplotlib dependencies used by the registered field-bundle/plotting tests and plot commands. C++ consumers linking only `gplspec` do not use these packages.

CTest runs the Stage 00 two-capture scientific comparison and the Stage 05 1-D/file-backed scientific comparisons. To see those comparisons directly, run `bash tests/run_stage00.sh build/bin/stage00_reference`; run a single Stage 05 gate with `bash tests/run_stage05_scientific_reference.sh stage05_1d build/bin/stage05_operator_reference tests/reference/stage05-1d-original.csv` or the corresponding `stage05_file` executable and frozen CSV. Each runner also prints the unchanged strict comparison as a diagnostic. The Python contract checks for missing/extra keys, non-finite values, and out-of-bound field/action changes run with `python3 tests/test_scientific_reference_gate.py`.

The benchmark executables resolve input paths relative to their current working directory. Bundle writers create output directories. Run them in a disposable tree, linking the tracked inputs for the readers. The symlink is a path convenience, not a filesystem write-protection mechanism; the programs read model inputs and write only below `work/`. For example:

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

The table above gives each program's source parameters and numerical settings; preserve those settings when selecting a different benchmark. `clean_bench_5` reads `lmax` and `maxstep`; `clean_bench_7` through `_10` read the mapping amplitude `h`. `phobos_gravity_perturbation` reads `lsmall`. Supply the intended values on stdin explicitly and record them with the run. The other listed programs use their tracked inputs and source settings without prompts. `OutputZiheng` consumes both `modeldata/TABLEA1.DAT` and `modeldata/prem.200.no`. Plot current bundles with `work/plot_fields.py` as shown in [the bundle plotting guide](work/plots.md).

These commands reproduce supported source workflows. The mapping from those runs to particular paper figures or ambiguous Phobos `lsmall` selections remains unverified unless separately documented; a successful run is not itself a paper-level claim.

## What validation establishes

### Phase 3 historical validation

The current Phase 4 validation evidence and its scope limits are recorded in the [Phase 4 handoff](docs/core-cleanup/phase-4-handoff.md). The following records the earlier Phase 3 validation.

The fresh disconnected Release configure and full build completed successfully; all 25 registered CTests passed, and the standalone multi-header/two-translation-unit link check exited 0. The Stage 00 and Stage 05 scientific gates passed their retained-key and quantity criteria. Unchanged strict comparators are still run and their nonzero diagnostic results are preserved separately.

Ten supported original-settings workflows (`clean_bench_1..4`, `_6`, `_9`, `_10`, `phobos_gravity`, `phobos_heterogeneous`, and `OutputZiheng`) exited 0 and produced all 40 expected files. Thirty-eight files met the initially sufficient conservative physical bound. The Phobos physical-radius output also passes its already-declared geometry tolerance, yielding 39 of 40 passing scientific comparisons. The remaining file is the clean_bench_10 classical `perturbed/MatrixSolution.out`; its `461813.485791 m^2/s^2` maximum difference is retained as a baseline Grid-view lifetime defect attribution, not described as parity or accepted by a relaxed tolerance. The other Phobos radius exact-check difference and all strict exact/bytewise diagnostic results remain visible in the evidence record.

Matched public API timings, including the candidate's independent residual calculation, were 14.026254 ms mapped and 342.357401 ms file-backed versus 13.893254 ms and 333.042666 ms for the old API lifecycle. Candidate medians changed by +0.96% and +2.80%; median RSS was 16,792 KiB versus 16,784 KiB mapped and 393,060 KiB versus 416,476 KiB file-backed. Candidate independent residuals were `2.69037115e-12` and `6.39413501e-10`; the solver algorithm and requested settings remain unchanged.

The exact commands, limits, output hashes, comparison scripts, and build/performance manifests are in [Phase 3F validation evidence](docs/core-cleanup/evidence/phase3f/README.md). These runs reproduce the listed source workflows and do not resolve ambiguous paper labels or interactive Phobos `lsmall` figure mapping.

### Historical Phase 2 validation

The following records the earlier Phase 2 evidence and remains historical.

The Stage 2B Release build compiled all configured targets. Its 14 CTests passed, the Stage 00 runner matched three frozen 20,700-record references with zero changed components and byte-identical representative output, and the standalone header-link smoke check exited successfully. The embedded CMake consumer also configured, built, and ran.

These checks do not execute the full `clean_bench_1..10` or Phobos workflows, do not exercise the interactive prompt values as a campaign, and do not validate plots or establish paper-level reproducibility. The bounded tests and saved evidence are linked from [the cleanup campaign record](docs/cleanup/campaign.md) and [Phase 2 evidence](docs/core-cleanup/evidence/phase2b-validation.log). Keep the frozen references and production inputs unchanged.
