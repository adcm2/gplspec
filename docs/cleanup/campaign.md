# GPLSpec cleanup campaign

The original cleanup stages consolidate duplicated implementation while preserving public interfaces, mathematical operations, solver behavior, coefficient conventions, and output formats. A later PR 07 correction pass is documented separately in [pr-07-corrections.md](pr-07-corrections.md); its density-normalization correction is outside the earlier cleanup-equivalence claims. No dependency upgrade, performance change, or GSHTrans modernization is part of the cleanup campaign.

## Stage status

| Stage | Branch / reviewed candidate | Status |
|---|---|---|
| 00 Baseline | `cleanup/00-baseline`, `dece12e59f86de59978e6cccfe684e22b4fcb06e` | Accepted; immutable original-source references and regression runner established. |
| 01 Hygiene | `cleanup/01-hygiene`, `3f7ccf46625d5b09d13d966c58573d437deb5c99` | Accepted after independent Luna review. |
| 02 Headers | `cleanup/02-headers`, `3337718bdd0a5a5746cb833b1b5371ea392c215a` | Accepted and integrated. |
| 03 Utilities | `cleanup/03-utilities`, `f429202bf5318adbd46995aea75bf6f3140838cf` | Accepted and integrated. |
| 04 Models | `cleanup/04-models`, `e352efcf09bcf0d0fbbbded62f1b218736e2fe79` | Accepted and integrated after independent review. |
| 05 Operators | `cleanup/05-operators`, `20a5ef194a22c3d8e71ae0900a1b0f9ce57f5289` | Accepted and integrated; Luna PASS and Sol PASS WITH NON-BLOCKING FINDINGS. Later documentation/build-configuration closure is recorded in the Sol addendum. |
| 06 Output | `cleanup/06-output`, `b813285fc34b13e9d4868860ae0651eb405457f4` | Human-approved and integrated into `cleanup/base`; paired transform helper retained unchanged. |
| 07 Final | `cleanup/base`, `a9409aed3902bf1dd4262ff7cdb62287fd4ea2e1` | Accepted and integrated; annotated tag `gplspec-clean-baseline-v1` points here. Luna PASS; Sol PASS WITH NON-BLOCKING FINDINGS. |

Earlier stage reports remain the detailed record of their individual diffs and reviews. The stage-07 close-out later fast-forwarded `cleanup/base` to `a9409aed3902bf1dd4262ff7cdb62287fd4ea2e1` and created annotated tag `gplspec-clean-baseline-v1` at that commit. Current local and remote `cleanup/base` resolve to the same commit; `main` and `develop` were not changed.

## Behavior and regression contract

The original-source baseline is commit `4ef3a66c62d64408d99989dd51c3ccbdc46d0b0c`, with toolchain and dependency provenance in [baseline-manifest.md](baseline-manifest.md). Stage-00 covers five finite deterministic model/perturbation fixtures, intermediate arrays, forces, matrix-free operator results, solutions, and a representative output file. Repeated reference runs were byte-identical on the recorded stack. Cross-stage numerical comparisons use `rtol=1e-13` and `atol=1e-15`; the representative file is compared byte-for-byte. Frozen references and tolerances are not regenerated or relaxed.

The full cleanup suite includes four public-header smoke translation units, a two-translation-unit link executable, stage-03/04/05/06 focused checks, original-source operator and output comparisons, and the stage-00 frozen runner. The exact current commands and fresh candidate results are recorded in [stage-07-handoff.md](stage-07-handoff.md) after validation.

## Accepted stage 06

Stage 06 extracted common output-angle and Wigner-matrix work, plus a paired scalar transform helper, and clarified one diagnostic label from `Error:` to `Requested tolerance:` while printing the same tolerance. It preserves `ForwardAndExpandRealScalarPair`, caller allocation and density selection, and caller-specific slice-rotation angle behavior. Its separate Luna review passed. The user approved the final stage-06 diff and authorized its commit and push on `cleanup/06-output`; the accepted commit is now the stage-07 base. See [review-stage-06.md](review-stage-06.md) and [stage-06-handoff.md](stage-06-handoff.md).

## Stage 07 boundary

Stage 07 reconciles user-facing build instructions and the accumulated campaign status. It does not modify production code, tests, build configuration, frozen references, or `experimental/`. The documentation records preserved Hermitian observations, angle and density conventions, finite fixture coverage, system dependency versions that are not enforced by CMake, and the setup required for embedded CMake consumers. The experiments under `experimental/`, including ellipticity calculations, remain untouched, excluded from builds, and numerically unvalidated.

The exact candidate files, validation provenance, branch inventory, independent review verdicts, and deferred limitations are summarized in [stage-07-handoff.md](stage-07-handoff.md). The user approved the final diff and authorized the coordinator to commit stage 07, fast-forward `cleanup/base`, create the annotated baseline tag, push and verify, then clean up only the explicitly listed branches after ancestry checks. The commit, tag, and push were completed; local and remote `cleanup/base` agree. Historical review inventories retain their capture-time state. Unrelated remote archive and `cleanup-2026-10-01/*` branches are retained.
