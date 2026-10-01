# GPLSpec cleanup campaign

This campaign extracts existing GPLSpec/GSHTrans functionality while preserving numerical behavior and supported public interfaces. Numerical bug fixes, new GSHTrans capabilities, solver/discretization changes, optimizations, batching, real-field compression, and output refactoring are excluded. The user has resumed stages 03–05, followed by a cumulative independent Sol high review and a human-approved stop. Stages 06–07, GSHTrans modernization, performance work, and pushes remain outside scope.

## Stage status

| Stage | Branch | Status | Acceptance |
|---|---|---|---|
| 00 Baseline | `cleanup/00-baseline` | Accepted after independent Luna review of `dece12e59f86de59978e6cccfe684e22b4fcb06e` | Pinned reproducible build, immutable original-source reference, repeatable five-fixture regression and output check |
| 01 Hygiene | `cleanup/01-hygiene` | Accepted after independent Luna PASS at `3f7ccf46625d5b09d13d966c58573d437deb5c99` | Documented deletion/interface decisions; examples and baseline regressions pass; Luna approval |
| 02 Headers | `cleanup/02-headers` | Accepted and integrated: reviewed candidate `3337718bdd0a5a5746cb833b1b5371ea392c215a` | Independent Luna high review passed; public include/link checks, examples, fail-closed patch guard, and frozen baseline regressions pass |
| 03 Utilities | `cleanup/03-utilities` | Accepted and integrated: reviewed candidate `f429202bf5318adbd46995aea75bf6f3140838cf` | Seven equivalent derivative builders included after review correction; focused comparisons and frozen regressions pass |
| 04 Models | `cleanup/04-models` | Accepted and integrated at `e352efcf09bcf0d0fbbbded62f1b218736e2fe79` | Eight storage initializers covered; four focused old/new checks, all examples/public headers, and frozen intermediate-array/solution regressions pass |
| 05 Operators | `cleanup/05-operators` | Accepted and integrated: code `20a5ef194a22c3d8e71ae0900a1b0f9ce57f5289`; Luna PASS, Sol PASS WITH NON-BLOCKING FINDINGS | Direct operator/source/perturbation comparisons and baseline regressions; human checkpoint remains after integration |

Stages begin from the latest accepted `cleanup/base`. Accepted stage branches are fast-forwarded only into `cleanup/base`. Keep `main` and `develop` unchanged; do not push. Each stage requires a fixed candidate and separate independent Luna high review before integration into `cleanup/base`. After stage 05 passes, request one cumulative Sol high review. Once stages 00–05 are integrated and coordinator tag/bookkeeping is complete, stop for the human-approved checkpoint; automated reviews do not constitute human approval. No work on stages 06–07 or GSHTrans modernization is authorized.

## Stage-00 comparison contract

The frozen reference is source commit `4ef3a66c62d64408d99989dd51c3ccbdc46d0b0c`, using the dependency revisions in `baseline-manifest.md`. It covers homogeneous and layered spheres, analytic lateral tomography values, a smooth aspherical mapping, a mapping perturbation, full model intermediate arrays, deterministic complex matrix-free operator products, force vectors, solutions, and `MatrixSolution.out`.

Repeated reference runs were byte-identical on the pinned compiler/dependency stack. Cross-stage checks use `rtol=1e-13` and `atol=1e-15`; no stage may relax these tolerances. The representative output file is compared byte-for-byte.

## Resumed campaign boundary

Stage 02 passed independent review and was integrated. The user has since authorized stages 03–05, with separate Luna review/integration per stage, followed by cumulative Sol high review and a human checkpoint. Stage 03 passed independent review at `f429202bf5318adbd46995aea75bf6f3140838cf` and was integrated. Stage 04 passed its independent Luna review at `689d6086738a9c7c740826945d26ec1f191a90b7` and was integrated at `e352efcf09bcf0d0fbbbded62f1b218736e2fe79`. Stage 05 was implemented on `cleanup/05-operators` from that accepted base. Fixed candidate `20a5ef194a22c3d8e71ae0900a1b0f9ce57f5289` passed independent Luna high review and cumulative Sol high review (non-blocking findings only); both reviews have no required code corrections. Sol evidence is in `review-stage-05-sol.md`. Stage 05 and its review records were fast-forward integrated through `9cc987d6bc46ae82762e0170e77c8ae3e06c5019`; final integration records are included in annotated tag `gplspec-cleanup-stage05`. Stages 00–05 are complete at the automated gates. Work is stopped awaiting human review and approval. No later stages are authorized.
