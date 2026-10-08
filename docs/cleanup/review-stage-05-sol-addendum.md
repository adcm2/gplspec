# Independent Sol stage-05 review addendum

**Verdict: PASS WITH NON-BLOCKING FINDINGS** for full reviewed commit `f0fa6920ccc9b53b1d000365fd6ff2f64b5ed700` on `cleanup/base`.

The previously accepted code candidate was `20a5ef194a22c3d8e71ae0900a1b0f9ce57f5289`; its documentation tip was `3b76dad4c0d16692d6f6670951298d78de77229c`. The post-candidate history adds documentation/review records, the stage-05 annotated tag at `83e1a9fb3021f40e6384278432f502baba40bcde`, and commit `f0fa6920ccc9b53b1d000365fd6ff2f64b5ed700`. The only post-tag production build changes are `CMakeLists.txt` and a comment in `cmake/eigen-system/CMakeLists.txt`. No `gplspec/`, `tests/`, or `examples/` production or test source changed between the reviewed candidate and current HEAD. Six new paths are confined to `experimental/ellipticity/`; documentation/status files also changed. Commit metadata and reflog can establish the Git event, not the human or tool that created it.

Reviewer invocation: independent agent configured with `model=gpt-6-sol`, `reasoning_effort=high`. Source review and focused checks used an isolated archive; no production edits were made during review.

## Git evidence

The intervening commits, in order, are `3b76dad4c0d16692d6f6670951298d78de77229c`, `9cc987d6bc46ae82762e0170e77c8ae3e06c5019`, `83e1a9fb3021f40e6384278432f502baba40bcde`, and `f0fa6920ccc9b53b1d000365fd6ff2f64b5ed700`. The first three add only documentation. The last has parent `83e1a9fb3021f40e6384278432f502baba40bcde` and changes nine files (893 insertions, eight deletions), combining the Eigen correction, six experimental files and status records. Its commit message describes the experiments but omits the Eigen correction.

Both HEAD and `cleanup/base` reflogs record `commit: Add experimental ellipticity comparison tools and documentation` at `2026-10-01 16:40:51 +0100`. Author and committer metadata both name Alex Myhill (`adcm2@cam.ac.uk`). These fields do not establish who or what executed the commit. The working tree was clean at review start. No attribution beyond this Git evidence is made.

## Findings

- The prior mandatory `find_package(Eigen3 3.4.0 EXACT CONFIG REQUIRED)` was removed. Absent a parent `Eigen3::Eigen` target, CMake declares the official Eigen 3.4.0 archive with SHA-256 `8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72`, populates headers without configuring Eigen's CMake project, and creates a global imported interface target for those headers. The downloaded archive hash and header version macros independently check as 3.4.0. `CMP0135 NEW` affects archive extraction timestamps, not library arithmetic.
- With a parent-provided `Eigen3::Eigen`, the local shim populates instead and does not replace the target. A prior parent fixture preserved its include-directory property and compiled a consumer linked to `gplspec`. It explicitly supplied C++23 and GPLSpec's FFTW module path; those are material conditions for this embedded test. The changed CMake still declares Eigen before GSHTrans, preventing GSHTrans's transitive declaration from replacing the selected source. GPLSpec's interface continues to propagate Eigen, GaussQuad, GSHTrans, FFTWpp, PlanetaryModel and TomographyModels.
- There is no hidden numerical or C++ public-interface change in the combined commit: `git diff` on `gplspec/`, `tests/`, and `examples/` is empty relative to the reviewed code. The added experiment is outside every production target and contains historical includes; its README explicitly says it is uncompiled, unrun and numerically unvalidated. All four historical file SHA-256 values match `provenance.json`. The source repository commit in that manifest was not independently rechecked in this addendum.
- Earlier build results predate the Git commit but the changed CMake file's filesystem timestamp predates those builds, and the compiled include paths select `_deps/eigen3-src`, not `/usr/include/eigen3`. To remove uncertainty about commit identity, I extracted `f0fa692` with `git archive` into `/tmp/gplspec-sol-addendum-f0fa692`, configured it offline against build-owned copies of the seven recorded pinned Git dependencies plus Eigen 3.4.0, then built and ran focused numerical and linkage targets. The seven copied dependency repository HEADs match the seven CMake `GIT_TAG` pins. The compiler flags point to the archived source and build-owned Eigen headers.

## Exact-commit validation

- Isolated offline CMake configure: PASS, CMake 3.28.3, GCC 13.3.0; FFTW 3.3.10 and NetCDF 4.9.2 in cache. The only configure diagnostic was the already known non-failing CMP0148/FindPythonInterp developer warning.
- Built `stage00_reference`, `stage05_operator_reference`, `stage05_file_perturbation_reference`, and the two-translation-unit `stage02_header_link`: PASS. Linked executable ran successfully.
- `tests/run_stage00.sh`: three 20,700-record comparisons each reported zero changed components; representative `MatrixSolution.out` byte-identical: PASS.
- `tests/run_stage05_1d_reference.sh` and `tests/run_stage05_file_perturbation_reference.sh`: byte-exact CSV comparisons against frozen original-source references: PASS (1,178 and 3,819 records respectively).
- The earlier fetched-Eigen full example build and nine CTests remain usable corroboration for identical precommit source and dependency configuration. The direct isolated build above is the exact-commit check; it did not rerun all examples or every CTest because production/test sources have not changed.

## Non-blocking findings and checkpoint

The previous findings remain: finite fixtures and degree-2 perturbation do not establish all-input or cross-platform equivalence; the original-source 1D/3D sesquilinear differences remain characterized rather than corrected; FFTW and NetCDF versions are recorded but not enforced by CMake. Parent-target reuse may select a parent Eigen version other than 3.4.0; that is intentional reuse, and the present numeric checks cover 3.4.0. Embedded consumers may need to supply the existing FFTW module path and C++23 setting, as the passing fixture did. The `implementation_status.md` lines saying this correction is uncommitted/in progress are stale at current HEAD and should be superseded in the checkpoint record. A one-line trailing whitespace defect in the preserved experimental plotting script is inherited from the historical copy and has no build effect.

No blocking correction is required. The existing annotated `gplspec-cleanup-stage05` tag peels to `83e1a9fb3021f40e6384278432f502baba40bcde`, not reviewed HEAD `f0fa6920ccc9b53b1d000365fd6ff2f64b5ed700`; do not move it. Mark the automated addendum complete for `f0fa692`, retain human review open, and stop at stage 05.

## Coordinator closure record

The automated stage-05 review is complete for the full commit above. No blocking fixes were requested, so no Luna correction/review cycle was needed. Stale present-tense status wording has been clarified in the uncommitted documentation update accompanying this addendum. These documentation edits are subsequent to the reviewed commit and do not change its code. The existing annotated tag (object `36565137c541ea6cbfddf77d3c724e1dde717ce9`) was left unchanged at its earlier target. Human review remains open. No commits, pushes, stages 06–07 or GSHTrans upgrades were performed during this closure.
