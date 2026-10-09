# Phase 4 final-validation handoff

Phase 4 automated validation is complete. Independent Luna and final cumulative Sol reviews PASS; the user's final human checkpoint remains open, so Phase 4 is not yet finally accepted. The reviewed implementation candidate is parent `637e28f7f2b99266c037ca5b127c5c0ae1274b50` plus the 18-file closeout freeze in `/tmp/gplspec-phase4/4f-closeout-freeze.json`. The reviewed 4F implementation is committed locally as `156c9cce7673b09ea7b91d50453404901c4c4fa3`. Phase 4E was accepted by independent Luna and committed locally as `637e28f7f2b99266c037ca5b127c5c0ae1274b50`.

The reviewed implementation is parent commit `637e28f7f2b99266c037ca5b127c5c0ae1274b50` plus the validated source/test working diff across `examples/OutputZiheng.cpp`, `examples/clean_bench_9.cpp`, `examples/clean_bench_10.cpp`, `gplspec/src/GeneralModels/Field_Bundle.h`, `tests/stage4c_bundle_fixture.cpp`, and `tests/test_field_bundle.py`. The complete eight-file 4F candidate diff at validation freeze had SHA-256 `b8fc442e2fd48c2819b9131bb43067d64df67731d6bbb060c31ace08193e36c1`; the other two files in that freeze were `docs/core-cleanup/phase-4.md` and `implementation_status.md`. This handoff and status reconciliation are documentation updates after that code freeze.

## Validation evidence

| Check | Result | Evidence |
|---|---|---|
| Fresh Release build of configured targets and header consumers at commit `637e28f` (before bounded 4F fixes) | PASS; evidence reused | `/tmp/gplspec-phase4/final-validation-manifest.json` |
| CTest and standalone header link at commit `637e28f` (before bounded 4F fixes) | 26/26 PASS; PASS; evidence reused | `/tmp/gplspec-phase4/final-ctest.log`, final validation manifest |
| Embedded parent-Eigen configure, build, and link at commit `637e28f` (before bounded 4F fixes) | PASS; evidence reused | `/tmp/gplspec-phase4/final-embedded-results.json` |
| Ten unchanged-settings workflows | 10/10 PASS across combined evidence; three affected workflows rerun after fixes | `/tmp/gplspec-phase4/final-workflows/manifest.json` |
| Forty historical-product scientific comparisons | 40/40 PASS at predeclared criteria across combined workflow evidence | `/tmp/gplspec-phase4/final-workflow-comparison.json` |
| Representative plot outputs | 9 rendered; coordinator visual checks passed | `/tmp/gplspec-phase4/final-plots/manifest.json` |
| Protected-file audit | All 97 protected hashes unchanged | `/tmp/gplspec-phase4/final-protected.json` |

After the bounded 4F fixes, `OutputZiheng`, `clean_bench_9`, `clean_bench_10`, and `stage4c_bundle_fixture` were rebuilt in the final build tree. The surface-row strict bundle regression passed 1/1. The three affected workflows were rerun; other final-validation build/test evidence and unaffected workflow captures above were reused rather than rerun on the post-fix working tree.

The coordinator visually checked the layered trace plot (both density jumps remain separated), Bench 3 analytic overlay, Phobos rotated surface and density slices, degree-power plot, and actual OutputZiheng topography. OutputZiheng's first isolated run failed on its retired `.out` writers; that capture is preserved. The writers were removed, optional `radial_node` now serializes blank for its `-1` sentinel, the surface bundle carries the model reference-surface radius, and a surface-row regression passes the strict reader. The three affected workflows were rerun after these corrections.

Supplemental review-response validation is separate from the ten captured workflows, 40-product comparison, and nine plots above. Two additional isolated synthetic workflows, `clean_bench_7` and `clean_bench_8`, completed with stdin `0.001` and no `modeldata` symlink; their synthetic bundle input lists are empty. The Bench 7 100-row sweep table is byte-identical to a fresh historical executable capture from the immutable baseline archive (SHA-256 `4ef58c634501c9a75a10be006b717818b7cdda19e17f53b670145c45e65b2108`), with identical `h` values and zero maximum physical-value difference; the supplementary scientific comparison passed. The actual plot rendered to `/tmp/gplspec-phase4/sweep-candidate/bench7_sweep.png` and was visually checked. This adds one supplementary table and one plot; it does not alter or regenerate any of the frozen 40 products or original nine plots. Evidence: `/tmp/gplspec-phase4/sweep-candidate/manifest.json`, `comparison.json`, and `/tmp/gplspec-phase4/sweep-baseline/manifest.json`.

The 40-product comparison matches the historical single-trace physical values and coordinates under the predeclared criteria. The exact both-side sample IDs/counts and continuity behavior are separately covered by `tests/stage4b_field_evaluation.cpp`, `tests/pr07_rotated_output.cpp`, and `tests/stage4c_bundle_fixture.cpp`; the 40-product comparison is not the evidence for those additions. The workflow manifest top-level source hashes belong to the original seven reused captures. Each of the three rerun records has its own current source map and executable hash, and `/tmp/gplspec-phase4/final-evidence-identity.json` records the current snapshot.

The comparison audit also covers Bench 9/10's exact-model coordinates. The legacy files paired mapped exact values with unmapped base-mesh radii. New exact bundles use mapped `testtomo2` coordinates; the baseline radii correspond to their reference coordinates, while independently checked physical radii match each original map (`h*r` and `h*cos(theta)*r`) including the original exterior taper. This is a coordinate-label correction, not a same-physical-coordinate comparison. No solver setting or comparison threshold changed. For OutputZiheng's outer harmonic vector, the addition-theorem pointwise error upper bound was `9.00762e-13`, below the conservative `GM/R` absolute floor `4.31071e-11`; baseline printed volume/density/radius rounding endpoints were accounted for. The comparison procedure is recorded in `/tmp/gplspec-phase4/compare-final-workflows.py`.

## Performance evidence and tradeoffs

The candidate performance campaign used 30 one-thread processes, three independent processes per case, and repeated batches of at least 0.4 seconds. The separate baseline campaign used 18 processes. Full measurements are in `/tmp/gplspec-phase4/performance-comparison.json`, with raw runs in `/tmp/gplspec-phase4/baseline-performance-runs/` and `/tmp/gplspec-phase4/candidate-performance-runs/`.

| Workload | Baseline | Candidate | Change |
|---|---:|---:|---:|
| File-backed solve | 0.351453 s | 0.357614 s | +1.75% |
| File-backed matched evaluation/export | 0.031612 s | 0.050188 s | +58.8% wall time |
| File-backed matched complete solve-to-output | 0.382681 s | 0.406740 s | +6.3% |
| Candidate complete path with both interface traces | — | 0.451890 s | +18.1% against baseline complete path; rows rise from 28,560 to 55,440 |

File-backed peak RSS remained about 394 MB and the residual remained `6.394135e-10`. Mapped CSV plus manifest grew from 11,291 to 29,388 bytes and throughput rose from 47 to 80 MB/s. File-backed output grew from 2,495,606 to 4,757,339 bytes; throughput rose from 79 to 95 MB/s, while measured export wall time also rose. For the mapped fixture, solve time was 0.01620 s baseline and 0.01494 s candidate; matched complete time was 0.01461 s baseline and 0.01618 s candidate, while the candidate both-trace path was 0.01580 s. These small-workload measurements are timing-sensitive. Report the richer records and additional traces with their actual wall cost; no performance pass is claimed and no speculative optimization was added.

## Remaining limits and checkpoint

The Stage00 scientific gate passes, but its unchanged strict arithmetic comparator still reports `homogeneous/operator_y/109/0` and the historical `.out` byte comparison differs. The retired Stage06 writer-reference harness had no direct scientific correspondence audit for its seven historical products; frozen hashes and strict comparator remain preserved. These limitations are disclosed and were not hidden by changing references or thresholds.

Independent Luna and cumulative Sol final verdicts are PASS for the reviewed candidate. Carried-forward limits remain: the unchanged Stage00 strict arithmetic and historical byte diagnostics differ; the seven retired Stage06 writer products lack direct product-by-product scientific correspondence audit; the protected experimental tree remains unbuilt and unvalidated. Performance is reported without a pass claim: file-backed solve +1.75%, matched evaluation/export +58.8%, matched complete path +6.3%, and mandatory both-trace complete path +18.1% with rows increasing from 28,560 to 55,440. The reviewed 4F implementation is committed locally as `156c9cce7673b09ea7b91d50453404901c4c4fa3`. Stop for the user's human Phase 4 checkpoint; do not start Phase 5, push, merge, tag, or delete branches.


## Final review disposition — human checkpoint open

Cumulative Sol final verdict: **PASS** for the frozen Phase 4 candidate, subject to the user's final human checkpoint. Independent Luna 4F and affected sweep reviews also PASS. Automated Phase 4 work is complete; the reviewed implementation is committed locally as `156c9cce7673b09ea7b91d50453404901c4c4fa3`. The 18-file reviewed closeout identity and patch hash are recorded in `/tmp/gplspec-phase4/4f-closeout-freeze.json`; this status-only closeout delta follows that reviewed freeze. Earlier REQUEST CHANGES reports and responses remain in the ledger.

Carried-forward limits and performance costs are unchanged: Stage00 strict arithmetic and historical `.out` byte diagnostics differ; seven retired Stage06 products lack direct product-by-product scientific correspondence audit; experimental remains unbuilt/unvalidated. File-backed solve was +1.75%, matched evaluation/export +58.8%, matched complete path +6.3%, and mandatory both-trace complete path +18.1% versus baseline, with rows increasing from 28,560 to 55,440. No performance pass, bitwise equivalence, or experimental validation is claimed. Stop for human review before Phase 5; no push, merge, or tag changes are authorized.
