# Phase 3F final validation

This evidence belongs to the candidate based on accepted Phase 3 commit `9a4f6a11449716560d52638442d0f0cd2acf766d` and baseline source `15582edab8112dbb63e0101be5b26f5135298644`. The source/test freeze hashes are in `code-freeze.json`; `protected-manifest.json` records the 121 protected paths, including frozen references, model inputs, saved products, and `experimental/`. The protected manifest is rooted at resume commit `f06ed749a0060f91e048103013b371006aa5aacb`; the user data and frozen references were not changed.

`code-freeze.json` describes the source used for validation before the post-review explanation-only floor correction in `tests/scientific_reference_gate.py`; `post-review-comment.json` maps that validated hash to the current hash and records that the executable Python AST is unchanged.

## Build and focused gates

`validation-manifest.json` records the fresh disconnected Release configuration, full `cmake --build ... -j2` (exit 0, 204.5 s), complete CTest (25/25, exit 0, 12.23 s), and standalone multi-header/two-TU link check (exit 0). The registered Stage 00 and Stage 05 scientific gates pass. The Stage 00 gate retains 18,648 of 20,700 frozen rows per capture and retires only the 2,052 `dxi`, `dxilm`, and `df` intermediate rows. Stage 05 1-D retains all 1,178 rows; Stage 05 file perturbation retains 1,766 of 3,818 rows and retires the same three intermediate categories. The gate contract's missing/extra, non-finite, excessive action, and excessive physical-potential rejection checks pass 6/6. The fixed potential floor `6.25631e-5 m^2/s^2` is conservative because fixture mass is at least `0.8*(4*pi/3)*M_norm`, greater than `1.2*M_norm`, and the modeled outer radius is at most `1.2*R_norm`.

Strict exact comparisons remain unchanged diagnostics and their mismatches are preserved: Stage 00 strict comparison and byte comparison each report differences for both captures, while two current Stage 00 captures compare strictly equal to each other; Stage 05 1-D and file strict comparisons also report differences. The scientific gates pass independently. No reference was regenerated or tolerance changed.

## Reproducible workflows

`workflow-manifest.json` gives the ten original-settings executions, inputs, exit codes, output paths, and hashes. The runs are `clean_bench_1..4`, `_6`, `_9`, `_10`, `phobos_gravity`, `phobos_heterogeneous`, and `OutputZiheng`; all ten exited 0 and produced 40 output files with no missing or extra files. The comparison settings and category criteria are in `workflow-comparison-contract.json`; the full per-file comparison is in `workflow-comparison.json` and the concise result log is `workflow-comparison.log`.

Thirty-eight files pass the initial conservative sufficient bound. The computed physical-radius array in `phobos_gravity/work/Phobos/Spherical/MatrixSolutionPhysical.out` differs at 36,369 points by at most `1.0914e-11 m`; the angles are exact and maximum potential difference is `1.2790e-13 m^2/s^2`. The initial contract's `10,000 m` radius scale covers Phobos surface-radius quantities such as `OutputZiheng` data. This computed volume-coordinate file extends to `R_out=1.5*10,000 m=15,000 m`, so the same predeclared geometry formula uses `atol=1e-14*R_out=1.5e-10 m` for this quantity. This is a scale application, not a changed tolerance; the measured difference passes both scales. `coordinate-comparison.json` preserves the original exact-check failure and the criterion-based result. The sole remaining workflow comparison exception is `clean_bench_10/work/Bench10/perturbed/MatrixSolution.out`, with maximum potential difference `461813.485791 m^2/s^2`. It is not reported as parity or a tolerance pass: a controlled baseline-only `GSH_Grid()` const-reference fix reproduced the accepted 3D parent result, attributing this difference to a temporary borrowed grid view in the original baseline path. The historical exact mismatch remains in the comparison evidence, and its bounded source attribution is in [`phase3e3/README.md`](../phase3e3/README.md).

All runs reproduce the supported source workflows and their preserved inputs/parameters. They do not establish ambiguous paper labels or an interactive Phobos `lsmall` figure mapping. See the current commands and isolated-output recipe in [benchmark.md](../../../../benchmark.md).

## Public API performance

The matched baseline timer includes the old `FindGravitationalPotential` nested-field conversion. The candidate timer covers `SolvePotential`, including its independent residual calculation, and stops before the explicit `ToHarmonicField` conversion used later by output workflows. With one thread and three process samples, median-of-process-medians were:

| Fixture | Baseline | Candidate | Change | Median RSS baseline → candidate |
|---|---:|---:|---:|---:|
| mapped | 13.893254 ms | 14.026254 ms | +0.96% | 16,784 → 16,792 KiB |
| file-backed | 333.042666 ms | 342.357401 ms | +2.80% | 416,476 → 393,060 KiB |

Candidate residuals were `2.69037115e-12` mapped and `6.39413501e-10` file-backed. A separate single-call capture confirms status `Eigen::Success`; the mapped solve reports 13 iterations, while the file-backed Eigen solve reports 0 displayed iterations. Their recursive errors were `4.94505538e-13` and `4.88398849e-21`; independent relative residuals are shown above. These values do not change the accepted solver algorithms or the workflow-specific residual interpretation.

The public API timing sources, baseline/candidate build manifests, per-process summaries, all six raw timing CSVs and RSS files for each side, comparison, and single-call diagnostics are kept here. `ctest.log` records the full 25-test run; `strict-diagnostic-manifest.json` and the short `strict-*.log` files preserve exact-comparison outcomes and first differences without copying the multi-megabyte raw capture trees. The timing comparison includes output-conversion scope differences noted above and does not reuse the older manual Eigen-lifecycle timing as an API measurement.
