# Stage 07 handoff

**Implementation status: Stage 07 was committed, integrated, tagged, and pushed.** Historical stage-07 close-out commit: `a9409aed3902bf1dd4262ff7cdb62287fd4ea2e1`; annotated tag `gplspec-clean-baseline-v1` points here. Luna high: PASS after resolving one low documentation finding. Cumulative Sol high: PASS WITH NON-BLOCKING FINDINGS (inherited limitations only). The branch inventories below are preserved as capture-time records.

## Starting point and exact scope

- Branch: `cleanup/07-final`.
- Starting HEAD: `b813285fc34b13e9d4868860ae0651eb405457f4`.
- Before creating `cleanup/07-final`, the coordinator fast-forwarded `cleanup/base` from `3045834c55c2906c1d6f60b660feb4e922618cae` to the accepted stage-06 commit `b813285fc34b13e9d4868860ae0651eb405457f4` under the user's authorization. `cleanup/07-final` was then created at that updated base.
- Changed files: `AGENTS.md`, `README.md`, `getting-started.md`, `docs/cleanup/campaign.md`, `docs/cleanup/deferred-issues.md`, `implementation_status.md`, and this handoff; final review records add `docs/cleanup/review-stage-07.md` and `docs/cleanup/review-stage-07-sol.md`.
- The change reconciles the repository's contributor guidance, public build/consumer instructions, current campaign status, preserved limitations, and validation evidence. It changes no production implementation, tests, build files, frozen reference, or content under `experimental/`.
- There was no demonstrably redundant production code in scope for removal. Stage 06's `ForwardAndExpandRealScalarPair` helper is retained as human-approved.

## Build and dependency provenance

Fresh candidate build: `/tmp/gplspec-stage07-build`, configured out of source with examples and the baseline harness enabled. All source dependencies came from the already populated stage-06 cache, were copied into the new build's own `_deps` tree, and were checked at their pinned revisions before building. This avoided external source overrides that the ODR patch guard correctly rejects. Dependency sources were not upgraded.

The candidate configure can be reproduced from the populated stage-06 cache with:

```sh
mkdir -p /tmp/gplspec-stage07-build/_deps
for x in eigen3 numericconcepts fftwpp gaussquad interpolation planetarymodel gshtrans tomographymodels; do
  cp -a "/tmp/gplspec-stage06-build/_deps/${x}-src" "/tmp/gplspec-stage07-build/_deps/${x}-src"
done
cmake -S . -B /tmp/gplspec-stage07-build \
  -DMY_PROJECT_BUILD_EXAMPLES=ON \
  -DGPLSPEC_BUILD_BASELINE_HARNESS=ON \
  -DFETCHCONTENT_SOURCE_DIR_EIGEN3=/tmp/gplspec-stage07-build/_deps/eigen3-src \
  -DFETCHCONTENT_SOURCE_DIR_NUMERICCONCEPTS=/tmp/gplspec-stage07-build/_deps/numericconcepts-src \
  -DFETCHCONTENT_SOURCE_DIR_FFTWPP=/tmp/gplspec-stage07-build/_deps/fftwpp-src \
  -DFETCHCONTENT_SOURCE_DIR_GAUSSQUAD=/tmp/gplspec-stage07-build/_deps/gaussquad-src \
  -DFETCHCONTENT_SOURCE_DIR_INTERPOLATION=/tmp/gplspec-stage07-build/_deps/interpolation-src \
  -DFETCHCONTENT_SOURCE_DIR_PLANETARYMODEL=/tmp/gplspec-stage07-build/_deps/planetarymodel-src \
  -DFETCHCONTENT_SOURCE_DIR_GSHTRANS=/tmp/gplspec-stage07-build/_deps/gshtrans-src \
  -DFETCHCONTENT_SOURCE_DIR_TOMOGRAPHYMODELS=/tmp/gplspec-stage07-build/_deps/tomographymodels-src
cmake --build /tmp/gplspec-stage07-build -j2
ctest --test-dir /tmp/gplspec-stage07-build --output-on-failure
/tmp/gplspec-stage07-build/bin/stage02_header_link
tests/run_stage00.sh /tmp/gplspec-stage07-build/bin/stage00_reference
```

The build and stage-00 command transcripts were streamed by the execution tool and not saved as separate shell logs. CTest detail is persisted at `/tmp/gplspec-stage07-build/Testing/Temporary/LastTest.log`. The stage-00 runner removes its temporary outputs on exit. The parent-Eigen consumer compile flags are at `/tmp/gplspec-stage07-consumer-parent-eigen-real/build/CMakeFiles/my_consumer.dir/flags.make`; consumer executables and CMake caches are under `/tmp/gplspec-stage07-consumer-{parent-eigen-real,fetched-eigen}/build`. The app stdout logs and output files are under `/tmp/gplspec-stage06-bench/run` and `/tmp/gplspec-stage07-bench/run`; the output hashes are recorded below.

- Compiler: `/usr/bin/c++`, GNU 13.3.0 (`13.3.0-6ubuntu2~24.04.1`); CMake 3.28.3; C++23; OpenMP 4.5; Python 3.12.3.
- System packages: FFTW 3.3.10 and NetCDF 4.9.2. These versions are recorded test provenance and are not enforced by the project's CMake configuration.
- Eigen: 3.4.0. The configured build used the previously fetched pinned source tree. The upstream archive pin is SHA-256 `8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72`.
- Git source pins: NumericConcepts `888126b44a979fda5dcc9d0b11d80cb8b3f61b52`; FFTWpp `06cc1fb04c4398407839e637cb5414b179b0e9c3`; GaussQuad `fbe37c7eef93695317dbf9c48715110bf64c800f`; Interpolation `1557aad571ceac1899c2d3795696aa4047f41d8a`; PlanetaryModel `a801b9e9c2205c3b32daa745fd523d596a9d2ede`; GSHTrans `f0a0e24a3579ceb7c99f8782328cae14bef195de`; TomographyModels `97869aeea3f3d87901705db59b5de975960bd9e1`.

An additional attempt to exercise the Eigen archive URL directly during this stage stalled with a zero-byte download in the restricted environment and was stopped. The previous stage-05 Sol addendum records successful SHA-pinned Eigen fetching, including version 3.4.0; the stage-07 no-parent consumer check below builds using that same pinned source tree. The consumer variant with a parent-provided `Eigen3::Eigen` target was independently reconfigured without an Eigen source override and compiled using the parent's `/usr/include/eigen3` include path.

## Validation results

- `cmake --build /tmp/gplspec-stage07-build -j2`: passed all configured targets, including all examples, four public-header smoke targets, the two-translation-unit link target, focused regression executables, and the baseline harness.
- `ctest --test-dir /tmp/gplspec-stage07-build --output-on-failure`: all 11 tests passed, including stage-05 operator/perturbation references and stage-06 output comparisons.
- `/tmp/gplspec-stage07-build/bin/stage02_header_link`: passed.
- `tests/run_stage00.sh /tmp/gplspec-stage07-build/bin/stage00_reference`: passed. Both 20,700-record runs matched the frozen original-source baseline with zero changed components, the two candidate runs matched each other, and both representative `MatrixSolution.out` files were byte-identical to the frozen output. No reference or tolerance was regenerated or changed.
- The exact minimal-consumer CMake example in `getting-started.md` built and ran in two configurations: with a parent-supplied `Eigen3::Eigen` and without one. The parent-target compile line retained `/usr/include/eigen3`; the no-parent configuration used the locked Eigen 3.4.0 source tree. Both also exercised the required FFTW module path and C++23 setting.
- `clean_bench_1` from `/tmp/gplspec-stage06-build/bin/clean_bench_1` and the fresh `/tmp/gplspec-stage07-build/bin/clean_bench_1` were run in separate disposable working directories, each with `work/Bench1` created. The accepted stage-06 build cache points to this repository, uses `/usr/bin/c++`, and had examples plus the baseline harness enabled. Stage 07 changed no production/test/build source from the accepted commit. All three output files matched byte-for-byte: `ExactSolution.out` SHA-256 `605931f6ec08ce9682da77c53245ed7710461ed764ae15d817d12c1afb071bd9`; `IntegralSolution.out` `e417783dc150672a33b7248198563ef8fff753103268d9334ffc92588c98cb91`; `MatrixSolution.out` `48925cac486b020ed2d16feb01693582fc1c937f5fbc9a16bfead033230784c1`.
- The pre-edit SHA-256 inventory of tracked production, test, build-configuration, and example files was saved at `/tmp/gplspec-stage07-code-before.sha256`; the final candidate was checked against it. Frozen reference files were unchanged. The tracked `experimental/` tree was unchanged and was not built or run.

## Preserved limitations

The deferred Hermitian-form discrepancies and their original-source values are in [deferred-issues.md](deferred-issues.md). The historical stage-07 candidate retained caller-specific rotation-angle rules and the existing density behavior; the later density-scale correction is documented in [pr-07-corrections.md](pr-07-corrections.md). Evidence covers finite deterministic fixtures on the toolchain above; it does not prove behavior for arbitrary models, degrees, or platforms. System FFTW/NetCDF versions are recorded rather than enforced. Experimental ellipticity sources remain preserved, excluded from builds, and numerically unvalidated. No GSHTrans modernization or scientific correction is included.

## Branch inventory captured before close-out (historical)

At the review checkpoint, the coordinator's ancestry audit found local branches `cleanup/00-baseline`, `cleanup/01-hygiene`, `cleanup/02-headers`, `cleanup/03-utilities`, `cleanup/04-models`, `cleanup/05-operators`, and `cleanup/06-output` to be ancestors of `cleanup/base` at stage-06 commit `b813285fc34b13e9d4868860ae0651eb405457f4`. `cleanup/07-final` was attached to an active worktree with candidate changes. After the user authorized close-out, the coordinator integrated stage 07 and removed the eligible cleanup branches following ancestry checks. Within the cleanup-stage branch inventory, the retained branch is `cleanup/base` locally and `origin/cleanup/base` remotely; unrelated branches are outside this inventory.

At that review checkpoint, the live remote heads audit found only `origin/cleanup/06-output` at `b813285fc34b13e9d4868860ae0651eb405457f4` and `origin/cleanup/base` at `3045834c55c2906c1d6f60b660feb4e922618cae`. At that time the stage-06 remote branch was proposed for deletion after final integration and another ancestry check; `main`, `develop`, unrelated `gh-pages`, and unrelated remote archives were to be retained. The baseline tag had not yet been created at that capture point. Current close-out status is recorded above.

## Independent review and final human gate (historical review snapshot)

Independent [Luna high review](review-stage-07.md) passed after a low Eigen-provenance wording ambiguity was clarified. The original 262-file manifest was `d5cbb774aa38ac03e662b34545fbefcec22ea20e40b745040fbd42adc19793e4`; Luna then independently checked the two-file documentation correction. The revised 262-file snapshot reviewed by both Luna and [cumulative Sol high](review-stage-07-sol.md) has manifest SHA-256 `576f7f0c2681a817941f79d1087707c33e74bf51741a5eec9d5ae5aa69a8a86d`. Sol's verdict is **PASS WITH NON-BLOCKING FINDINGS**, limited to the inherited numerical, coverage, and dependency limitations above. No blocking finding remains.

After those fixed-snapshot reviews, the coordinator added the two review reports and updated only this handoff, `campaign.md`, and `implementation_status.md` to record the verdicts and final human gate. This final documentation-only delta is distinct from the implementation snapshot identified above. The complete final working-tree manifest and patch, including the new review files, are at `build/cleanup-review/stage07-final-manifest.json` and `build/cleanup-review/stage07-working-tree.patch`; the accepted implementation-review manifest is retained at `build/cleanup-review/stage07-sol-reviewed-manifest.json`. These are local ignored review artifacts, not source changes. The narrow Sol follow-up passed for that final record delta; it does not cover the later user-requested documentation corrections.

**Historical capture note:** the branch inventory and pre-close-out statements above describe the repository at the stage-07 review checkpoint. The reviewed work was later committed as `a9409aed3902bf1dd4262ff7cdb62287fd4ea2e1`, pushed to `cleanup/base`, and tagged with annotated tag `gplspec-clean-baseline-v1` (tag object `2591d902f6e5414703d9d2a4435302cb93e8f4a3`). At stage-07 close-out, local and remote `cleanup/base` resolved to that commit. After the user-authorized close-out and ancestry checks, the listed cleanup branches were deleted; unrelated remote archive and `cleanup-2026-10-01/*` branches were retained. GSHTrans modernization remains excluded.
