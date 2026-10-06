# PR 07 correction pass handoff

This correction pass follows the read-only full-PR review at `/tmp/gplspec-pr7-review/sol-review.md`, which compared base `4fd6488a2b8aa989e92fd04f91ad28bb816ed359` with head `a9409aed3902bf1dd4262ff7cdb62287fd4ea2e1`. Existing review records were left unchanged. The approved scope was limited to the rotated boundary convention and its tests/documentation, the Phobos tutorial, Stage-07 current-state records, MathJax ordering, density normalization, and a PR description. This correction pass made no commit, push, tag, merge, GitHub review/comment, or live PR-description edit.

## Changes made

- Retained the established inner-side trace at duplicated radii in `ReferentialOutputRotated`, `PhysicalOutputRotated`, and `ModelDensityOutputRotated`. Source comments explain the preceding element's upper node, intentional density plotting through the surface, and unchanged first/final endpoints.
- Corrected `ModelDensityOutputRotated` to multiply by `DensityNorm()`. Referential and physical potential output continue to use `PotentialNorm()`.
- Added `tests/pr07_rotated_output.cpp`, which builds a two-layer jump with a vacuum exterior, a nonunit radial mapping/Jacobian, and `DensityNorm()=8` versus `PotentialNorm()=4`. It checks internal traces and endpoints in all rotated writers, continuous-potential output, density scale, and `rho_ref = J rho_phys`.
- Added a build target that extracts and compiles the actual C++ fence from `_tutorials/tutorial5_phobos.md`. Its angle values and output paths are initialized, and its calls use declared vectors.
- Moved the MathJax configuration before the asynchronous loader.
- Corrected current Stage-07 campaign/handoff status to distinguish the historical review checkpoint from the completed close-out. The current local and remote `cleanup/base` resolve to `a9409aed3902bf1dd4262ff7cdb62287fd4ea2e1`; annotated tag `gplspec-clean-baseline-v1` points there.
- Added [pr-07-description.md](pr-07-description.md), which separates behavior changes predating `4ef3a66` from the cleanup and correction pass and states the evidence limits.

## Validation

- Focused `pr07_rotated_output` build and CTest: passed, 1/1.
- `tutorial5_phobos_compile`: passed; the generated source came from the Markdown fence.
- Fresh validation build: `/tmp/gplspec-pr7-validation-build`, configured with examples and the baseline harness enabled. Configure used the pinned dependency sources and GNU 13.3.0, CMake 3.28.3, Python 3.12.3, OpenMP 4.5, FFTW 3.3.10, and NetCDF 4.9.2.
- Fresh configure: `cmake -S . -B /tmp/gplspec-pr7-validation-build -DMY_PROJECT_BUILD_EXAMPLES=ON -DGPLSPEC_BUILD_BASELINE_HARNESS=ON`.
- Full build: `cmake --build /tmp/gplspec-pr7-validation-build -j2` passed all targets, including 15 examples, four public-header smoke objects, `stage02_header_link`, `pr07_rotated_output`, and `tutorial5_phobos_compile`. The configure transcript is `/tmp/gplspec-pr7-reconfigure.log`; the final build transcript is `/tmp/gplspec-pr7-final-build.log`.
- Complete CTest: `ctest --test-dir /tmp/gplspec-pr7-validation-build --output-on-failure` passed 12/12. The persisted result is `/tmp/gplspec-pr7-validation-build/Testing/Temporary/LastTest.log`.
- Stage-00 runner: `tests/run_stage00.sh /tmp/gplspec-pr7-validation-build/bin/stage00_reference` passed all three 20,700-record comparisons with zero changed components and maximum absolute difference 0; `MatrixSolution.out` was byte-identical to the frozen baseline. Its output was streamed during execution and is summarized in `implementation_status.md`; no separate runner log was saved.
- Public-header link: `/tmp/gplspec-pr7-validation-build/bin/stage02_header_link` passed. `git diff --check` passed and `git diff --exit-code -- tests/reference` returned clean. No reference files or tolerances were changed.

## Independent review and human gate

Luna v2 returned **PASS** after resolving one non-blocking documentation finding. Sol returned **PASS** with no findings. Both reviewed the frozen 268-file candidate-v2 manifest `075f5cfe1ee5f7064cf9e35042b07c729b6495c732c1a9167d909bf045919194`; their verbatim reports are [review-pr-07-luna.md](review-pr-07-luna.md) and [review-pr-07-sol.md](review-pr-07-sol.md). The post-review delta consists only of these two review reports and updates to this handoff and `implementation_status.md`; it is separate from the frozen candidate reviewed by both agents. No source, test, build, or reference file changed, and no validation was rerun. Awaiting human diff approval; no integration or GitHub action is authorized by these review reports.

## Evidence boundary

The retained boundary convention predates the cleanup baseline. The stage-06 saved output fixture has unit normalization scales and a smooth model, so it cannot detect the density norm correction or establish compatibility with `main` at discontinuous interfaces. `clean_bench_2` resolution, the Phobos solve tolerance, and `_vec_da` zero-initialization are earlier behavior changes, not cleanup equivalence claims. `OutputZiheng` output is only structurally inspected. The experiments under `experimental/` remain untouched and unvalidated.
