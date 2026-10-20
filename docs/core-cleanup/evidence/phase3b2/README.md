# Phase 3B.2 model and workflow evidence

## Owned model values

`run-prem-model-comparison.sh` recompiles the 3B.1 baseline probe against detached baseline `15582edab8112dbb63e0101be5b26f5135298644` and the post-migration probe in this directory against the candidate. Both runs use `-O0 -DNDEBUG -std=gnu++23 -fopenmp` and `OMP_NUM_THREADS=1`. It writes all four raw profile CSVs under `/tmp/gplspec-phase3b2-prem-repro/` and requires `cmp` equality. The baseline and current candidate profile hashes match for all PREM500 and `prem.200` normalized layer radii/material interpolants. Input, probe, output, and script hashes are in `manifest.json`.

`stage3b_model_input` checks the parsed PREM and TABLEA1 known values, one-sided PREM layers, scale/radius formulas, malformed input, owned-value lifetime after input-row clearing, and move-only radial-model safety. A compile-time concept check accepts `double` and rejects `float`, matching the model’s double-valued layer storage and spline iterator type. The moved radial model is evaluated after its source owner leaves scope. Returned spline and derived material callables borrow/capture the model; discard them before moving or destroying their owner. `stage05_file_perturbation_original_reference` remains an exact old-source comparison for the shared typed-row perturbation path.

## Build and tests

Before the final double-only scalar constraint and OutputZiheng timer-placement correction, `cmake --build /tmp/gplspec-phase3-candidate-release -j2` completed all configured targets, including examples/tutorials, public-header and two-translation-unit probes, and external consumer. Full CTest passed 16/16 in 10.69 s, and `tests/run_stage00.sh /tmp/gplspec-phase3-candidate-release/bin/stage00_reference` passed three 20,700-record comparisons with zero changed components and a byte-identical representative output. These are pre-final-correction results, not validation of the exact committed source. After those corrections, the focused three-target rebuild and 2/2 CTest passed; see `review-fix-build.log` and `review-fix-ctest.log`.

## Captured workflows

The coordinator ran five fixed workflows at one thread, from candidate source `71ee3f366c31951a39c1ff2eac015c9d077b9fd2` plus the initial 3B.2 implementation, before the final move-only owner, double-only scalar constraint, and timer-placement corrections. All five workflows exited 0 and all compared values were finite. Their output hashes and mismatch lists are in `workflow-manifest.json`; per-file numeric deltas and tolerance failures are in the adjacent JSONs. These are pre-final-correction source captures, not validation of every line in committed HEAD `0912fd62a2ecf1ebab6ab0c2b4fbc05ac9b14c45`; the current manifest separately records the later focused validation. No arithmetic change is implied by the timing-boundary correction.

`clean_bench_2` and `clean_bench_6` output files were byte-identical. In `OutputZiheng`, the Phobos radius and PREM integral were byte-identical; the outer potential differed by at most `3.552713678800501e-15` (47 components changed, zero failures under its existing numeric tolerance). Phobos gravity and heterogeneous outputs were not byte-identical and retain strict componentwise potential-field failures: gravity reached 74,231 failing components and heterogeneous reached 216,724 in referential-density output. Their outputs were finite; sensitivity and power files had no strict comparison failures. These are comparison results, not a blanket equivalence claim. The prior 3B.1 untouched-run variability remains context only; it does not explain every mismatch.

The tracked `experimental/ellipticity/examples/ellipticity_test.cpp` remains unchanged preserved experimental source and is not an active target. Higher-level `Density3D`/discretization ownership coupled to the solver layout is deferred with 3C. No solver accuracy or full Phase 3 completion is claimed.


## Luna review correction checkpoint

The first independent review marked this checkpoint **BLOCKED** on two interface/evidence defects. The radial owner now explicitly supports `double` only, matching its stored `vector<double>` profiles and spline iterators; the compile-time checks in `stage3b_model_input` accept `double` and reject `float`. Its returned spline callables borrow layer storage, while derived material functions capture the owner, so discard returned functions before moving or destroying the model. `OutputZiheng` now starts its timer before `ReadPhobosHarmonics`, including TABLEA1 parsing in the measured workflow. No formulas, solver behavior, or scientific outputs were changed.

Focused rebuild command: `cmake --build /tmp/gplspec-phase3-candidate-release --target stage3b_model_input OutputZiheng stage05_file_perturbation_reference -j2`. All three targets built. Focused CTest ran `stage3b_model_input` and `stage05_file_perturbation_original_reference`; both passed (2/2). Build and test transcripts, including hashes in the manifest, are `review-fix-build.log` and `review-fix-ctest.log`. The earlier full-workflow source hashes and comparisons above remain a historical pre-review checkpoint; they do not represent the timer-placement/scalar-constraint source revision. Independent Luna re-review: **PASS**, no implementation blockers. Stage 3B.2 was accepted pending coordinator commit at review time and was subsequently committed as `0912fd62a2ecf1ebab6ab0c2b4fbc05ac9b14c45`.


The cumulative Sol-high partial-checkpoint verdict and stop boundary are recorded in [`../../phase-3-handoff.md`](../../phase-3-handoff.md); its exact report is archived at `../phase3-final/sol-high-cumulative-review.md`.
