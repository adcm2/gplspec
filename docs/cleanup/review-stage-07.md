# Independent Luna high review — stage 07

**Verdict: PASS WITH ONE NON-BLOCKING DOCUMENTATION FINDING.** No code, scope, compatibility, provenance, or scientific-claim blocker was found. The finding is a wording ambiguity in the Eigen version sentence in `getting-started.md`.

## Reviewed identity and scope

- Base HEAD: `b813285fc34b13e9d4868860ae0651eb405457f4` (`cleanup/07-final`).
- Frozen candidate: `/tmp/gplspec-stage07-review-candidate`.
- Manifest: 262 paths; SHA-256 `d5cbb774aa38ac03e662b34545fbefcec22ea20e40b745040fbd42adc19793e4`.
- Patch: `/tmp/gplspec-stage07-review.patch`; SHA-256 `ea439b565042c3a12ddc72a282d27954a3ad28dfc1c87bb955d975f1712877b2`.
- Recomputed every candidate file against the supplied manifest: 262 entries, zero mismatches. Recomputed the same entries against the live source tree: zero mismatches. The live repository remains at the stated base HEAD with the candidate's expected documentation modifications and new handoff file.
- The only changed paths are the seven authorized documentation/status files: `AGENTS.md`, `README.md`, `getting-started.md`, `docs/cleanup/campaign.md`, `docs/cleanup/deferred-issues.md`, `implementation_status.md`, and `docs/cleanup/stage-07-handoff.md`. No production, test, build, frozen-reference, or experimental files changed.

## Finding

- **Low, non-blocking — `getting-started.md`, dependency paragraph:** “The tested system package versions were Eigen 3.4.0, FFTW 3.3.10, and NetCDF 4.9.2” groups Eigen with the system packages, while the next sentence says Eigen is fetched from the pinned archive when there is no parent `Eigen3::Eigen`. Both routes are valid, but the wording can imply that the stage-07 build used the system package or that the consumer must install Eigen. The stage-07 handoff says its fresh build used the previously fetched pinned Eigen source tree, while the historical stage-00 manifest says that original baseline used the system Eigen package. Clarify which validation stack the listed package versions describe, for example: “The stage-00 baseline used Eigen 3.4.0, FFTW 3.3.10, and NetCDF 4.9.2; Eigen can also be fetched from the pinned archive when no parent project supplies `Eigen3::Eigen`.”

## Review evidence

- The public consumer setup is consistent with the inspected CMake: GPLSpec calls `find_package(FFTW)` while configuring, so the parent adds the repository `cmake/` module path before `add_subdirectory`; `gplspec` carries the dependency interface targets; the parent sets C++23; and a parent `Eigen3::Eigen` is reused. The documented clone/build executable path and the example's `./work/Bench1` output requirement match the CMake and example source.
- The handoff distinguishes stage-07 fresh-build evidence from reused dependency and earlier archive-fetch evidence. It explicitly states that stage-07 used pinned dependency source copies from the accepted stage-06 cache, reused the pinned Eigen source tree, and that a new direct archive download attempt stalled at zero bytes. It does not present that stalled fetch as a successful stage-07 fetch.
- The reported consumer checks, stage-00 comparison counts, application outputs, and numerical coverage boundaries are precise and do not imply broad scientific validation. Deferred Hermitian observations, angle/density conventions, and excluded/unvalidated experimental sources are candidly retained.
- Read the persisted `/tmp/gplspec-stage07-build/Testing/Temporary/LastTest.log`: it records all 11/11 CTests passing, including seven stage-06 output files matching the original-source baseline. The handoff also records the successful all-target build, stage-00 runner, link executable, exact consumer configurations, and byte-identical three-file `clean_bench_1` comparison. No redundant build or test was run for this documentation-only review.

## Narrow re-review of the Sol candidate

**Resolution verdict: PASS; the prior non-blocking finding is resolved.** Rechecked the frozen `/tmp/gplspec-stage07-sol-candidate` snapshot against manifest SHA-256 `576f7f0c2681a817941f79d1087707c33e74bf51741a5eec9d5ae5aa69a8a86d`: 262 entries, zero mismatches. Compared with the original reviewed snapshot, only `getting-started.md` and `implementation_status.md` differ.

The revised dependency paragraph now distinguishes the stage-00 system-package record from stage-07's use of Eigen 3.4.0 from the pinned source archive, and says that an embedding parent may supply `Eigen3::Eigen`. This resolves the ambiguity without changing consumer requirements or validation claims. The implementation-status addition accurately records the review finding, this limited clarification, and that no test rerun was needed. No new issue was found in this two-file delta; no repository files were edited and no tests were run.
