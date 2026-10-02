# GPLSpec cleanup campaign

This campaign preserves numerical behavior and supported public interfaces. The user has authorized stage 06 output and diagnostics consolidation after approving its bounded proposal. Stage 06 is implemented, independently reviewed and human-approved. It is not integrated into cleanup/base. Numerical fixes, solver/discretization changes, dependency upgrades, performance work, stage 07 and GSHTrans modernization remain excluded. The user explicitly authorized the stage-06 commit and non-force push of cleanup/06-output. Integration into cleanup/base and tag changes remain unauthorized.

## Stage status

| Stage | Branch | Status | Acceptance |
|---|---|---|---|
| 00 Baseline | `cleanup/00-baseline` | Accepted after independent Luna review of `dece12e59f86de59978e6cccfe684e22b4fcb06e` | Pinned reproducible build, immutable original-source reference, repeatable five-fixture regression and output check |
| 01 Hygiene | `cleanup/01-hygiene` | Accepted after independent Luna PASS at `3f7ccf46625d5b09d13d966c58573d437deb5c99` | Documented deletion/interface decisions; examples and baseline regressions pass; Luna approval |
| 02 Headers | `cleanup/02-headers` | Accepted and integrated: reviewed candidate `3337718bdd0a5a5746cb833b1b5371ea392c215a` | Independent Luna high review passed; public include/link checks, examples, fail-closed patch guard, and frozen baseline regressions pass |
| 03 Utilities | `cleanup/03-utilities` | Accepted and integrated: reviewed candidate `f429202bf5318adbd46995aea75bf6f3140838cf` | Seven equivalent derivative builders included after review correction; focused comparisons and frozen regressions pass |
| 04 Models | `cleanup/04-models` | Accepted and integrated at `e352efcf09bcf0d0fbbbded62f1b218736e2fe79` | Eight storage initializers covered; four focused old/new checks, all examples/public headers, and frozen intermediate-array/solution regressions pass |
| 05 Operators | `cleanup/05-operators` | Accepted and integrated: code `20a5ef194a22c3d8e71ae0900a1b0f9ce57f5289`; Luna PASS, Sol PASS WITH NON-BLOCKING FINDINGS | Direct operator/source/perturbation comparisons and baseline regressions; human checkpoint remains after integration |

## Historical stages 00–05 workflow

Stages begin from the latest accepted `cleanup/base`. Accepted stage branches are fast-forwarded only into `cleanup/base`. Keep `main` and `develop` unchanged; do not push. Each stage requires a fixed candidate and separate independent Luna high review before integration into `cleanup/base`. After stage 05 passes, request one cumulative Sol high review. Once stages 00–05 are integrated and coordinator tag/bookkeeping is complete, stop for the human-approved checkpoint; automated reviews do not constitute human approval. No work on stages 06–07 or GSHTrans modernization is authorized.

## Stage-00 comparison contract

The frozen reference is source commit `4ef3a66c62d64408d99989dd51c3ccbdc46d0b0c`, using the dependency revisions in `baseline-manifest.md`. It covers homogeneous and layered spheres, analytic lateral tomography values, a smooth aspherical mapping, a mapping perturbation, full model intermediate arrays, deterministic complex matrix-free operator products, force vectors, solutions, and `MatrixSolution.out`.

Repeated reference runs were byte-identical on the pinned compiler/dependency stack. Cross-stage checks use `rtol=1e-13` and `atol=1e-15`; no stage may relax these tolerances. The representative output file is compared byte-for-byte.

## Historical stage-05 campaign boundary

Stage 02 passed independent review and was integrated. The user has since authorized stages 03–05, with separate Luna review/integration per stage, followed by cumulative Sol high review and a human checkpoint. Stage 03 passed independent review at `f429202bf5318adbd46995aea75bf6f3140838cf` and was integrated. Stage 04 passed its independent Luna review at `689d6086738a9c7c740826945d26ec1f191a90b7` and was integrated at `e352efcf09bcf0d0fbbbded62f1b218736e2fe79`. Stage 05 was implemented on `cleanup/05-operators` from that accepted base. Fixed candidate `20a5ef194a22c3d8e71ae0900a1b0f9ce57f5289` passed independent Luna high review and cumulative Sol high review (non-blocking findings only); both reviews have no required code corrections. Sol evidence is in `review-stage-05-sol.md`. Stage 05 and its review records were fast-forward integrated through `9cc987d6bc46ae82762e0170e77c8ae3e06c5019`; final integration records are included in annotated tag `gplspec-cleanup-stage05`. Stages 00–05 are complete at the automated gates. Work is stopped awaiting human review and approval. No later stages are authorized.

## Historical stage-05 addendum closure

Independent Sol high review accepted full commit `f0fa6920ccc9b53b1d000365fd6ff2f64b5ed700` with **PASS WITH NON-BLOCKING FINDINGS**. The addendum covers all changes since previous Sol candidate `20a5ef194a22c3d8e71ae0900a1b0f9ce57f5289`, including restored Eigen fetching and isolated experimental ellipticity files. Exact-commit numerical/link checks passed; parent-target reuse and retained full-build evidence were checked. See [the review addendum](review-stage-05-sol-addendum.md) for validation, remaining limitations and Git evidence.

The automated stage-05 gate is complete for this commit. **Human review remains OPEN.** The existing annotated `gplspec-cleanup-stage05` tag remains at `83e1a9fb3021f40e6384278432f502baba40bcde`; it does not identify the newly reviewed HEAD and was not moved. This closure adds only uncommitted documentation. Stop at stage 05; stages 06–07 and GSHTrans modernization remain unauthorized.

## Current stage-06 checkpoint — 2026-10-02

**Stage 06 human-approved.** The reviewed diff, including the unchanged paired-transform helper, is accepted for commit and push on cleanup/06-output only. Not integrated into cleanup/base.

- Branch: `cleanup/06-output`; starting HEAD: `3045834c55c2906c1d6f60b660feb4e922618cae`. At start, `cleanup/base` matched local `origin/cleanup/base`; no reset or fetch was performed. The user's stage-06 request supersedes the historical stop boundaries above.
- Implementer: Luna medium. Independent reviewer: separate Luna high, **PASS**, no blocking code findings. See [review](review-stage-06.md) and [human-review handoff](stage-06-handoff.md).
- Validation: all 11 CTests, public umbrella/link checks, three example builds, frozen stage00 regressions and seven baseline-derived output comparisons passed. Existing references were not regenerated. `experimental/` remains unchanged and excluded from production.
- Sole intended console change: `Error:` becomes `Requested tolerance:` with the same printed tolerance. Public interfaces and output-file formatting/values are preserved in the tested cases.
- The existing stage-05 tag remains at `83e1a9fb3021f40e6384278432f502baba40bcde`; no tags were changed. The user has approved committing and pushing cleanup/06-output; do not integrate into cleanup/base, change tags, start stage 07 or upgrade GSHTrans.

### Post-review clarification and finalization authorization

Human approval explicitly retains `ForwardAndExpandRealScalarPair` unchanged. After independent review, Luna supplied a test-only clarification: the assertion diagnostic now says `angle fixture does not distinguish original slice and output conventions`, with a comment identifying its fixture-check purpose and the separate production-output baseline test. The assertion expression, production code, references and checksums remain unchanged. This small correction is distinct from the historical independently reviewed snapshot; the existing full validation continues to apply to unchanged production code. The final helper test result is recorded in `implementation_status.md`. Commit/push authorization applies only to `cleanup/06-output`; stage 06 is human-approved, not integrated.
