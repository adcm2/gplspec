# GPLSpec cleanup campaign instructions

## Goal and scope

Build the simplest readable, maintainable implementation of GPLSpec’s scientific capabilities while preserving correctness, reproducibility and performance. Phases 3 and 4 are accepted. The user has authorized Phase 5 stages 5A–5E to retire remaining migration scaffolding and validate the completed reconstruction, under the scope and gates recorded in the current user request and current handoff.

Treat the existing implementation as a scientific and regression starting point, not a permanent architecture. Replace APIs, types, files and implementations when this simplifies the code. Preserve useful mathematical explanations. Avoid unnecessary abstractions, compatibility layers and duplicate implementations.

Follow the latest explicit user authorization and the active phase plan in `docs/core-cleanup/`. Historical instructions do not override newer authorization. Phases 2 and 3 are accepted. Earlier Phase 3 continuation proceeded from its reviewed partial checkpoint; the later user authorization explicitly approved the Phase 3 corrections and all Phase 4 stages 4A–4F. Phase 4 uses baseline `a8b1ca99fdadd9541a4227129c8fce872c98c1aa`; the Phase 3 baseline `15582edab8112dbb63e0101be5b26f5135298644` remains historical provenance.

## Execution and review

- Use one active writer: Luna implements; a separate independent Luna reviews each frozen stage.
- Resolve blocking findings, commit accepted stages locally, and proceed automatically within the authorized phase after each independent Luna stage review accepts the checkpoint. No separate permission is needed for those stage commits.
- Keep every accepted stage working and independently reviewable.
- At the Phase 5 boundary, obtain cumulative Sol-high review, resolve blockers, then stop for human review with the complete cumulative patch and validation evidence.
- Use additional Sol reviews only for concrete mathematical or architectural uncertainties.
- Do not silently substitute models if the requested model is unavailable.
- Do not push, merge, change tags, delete branches or rewrite history without explicit authorization.
- Preserve unrelated user changes and do not stage them implicitly.

## Scientific and simplification requirements

- Treat existing GPLSpec as the accepted scientific and algorithmic reference, not an architecture to preserve. Reconstruct its supported capabilities and fundamental algorithms through the simplest clear, maintainable and efficient implementation. Baseline agreement alone does not prove scientific correctness; use independent analytical checks and reproducible workflows.
- Preserve clean_bench and Phobos paper-verification inputs, parameters and meaningful comparisons. Their source code and APIs may be rewritten.
- Before replacing a component, state its mathematical contract and establish relevant regression and independent checks. Define physical comparisons with units, scale-aware absolute and relative tolerances, and a near-zero absolute floor before evaluating the replacement; do not tune criteria to obtain a pass. The clean_bench_3 and clean_bench_6 residual observations accepted by the user are workflow-specific and do not authorize solver-algorithm changes.
- Validate replacements, migrate callers, then remove superseded implementations. Temporary adapters must have an explicit removal stage.
- Unknown hypothetical consumers do not require indefinite compatibility. API changes within the approved migration scope may proceed.
- Retain existing names when they suit the simplest design; renaming or deleting every old header is not a goal.
- Assess historical prototypes for unique scientific value. Preserve useful independent checks and Git provenance without maintaining unnecessary archive copies.
- Do not regenerate frozen references, loosen tolerances, silently fix numerical discrepancies or upgrade dependencies as part of cleanup.
- Preserve datasets and user outputs unless an explicit preservation/retrieval decision authorizes movement or deletion.
- Leave `experimental/`, including `experimental/ellipticity/`, unchanged, unbuilt and unvalidated unless separately authorized.

## Validation

Use focused checks during stages and a fresh out-of-source Release build for final phase validation:

```sh
cmake -S . -B /tmp/gplspec-phaseN-validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DMY_PROJECT_BUILD_EXAMPLES=ON \
  -DGPLSPEC_BUILD_BASELINE_HARNESS=ON
cmake --build /tmp/gplspec-phaseN-validation -j2
ctest --test-dir /tmp/gplspec-phaseN-validation --output-on-failure
tests/run_stage00.sh /tmp/gplspec-phaseN-validation/bin/stage00_reference
/tmp/gplspec-phaseN-validation/bin/stage02_header_link
```

Replace `phaseN` with the active phase and use a new build directory.

- Establish transform indexing, spin, normalization and grid checks before rewriting dependent numerical code.
- Capture baseline scientific and performance results before replacement; preserve the baseline source at its immutable Git commit and use an isolated source archive/build when a detached worktree is unavailable.
- Check operator actions and independently recomputed solver residuals. Retain genuinely independent analytic, integral and perturbation checks.
- At numerical-phase boundaries, run the identified paper-verification cases at their recorded settings and compare scientific outputs.
- Compare performance using matched Release builds, fixed inputs and thread settings, repeated substantial workloads, memory measurements, iterations and residuals. Investigate material regressions; do not invent thresholds from noisy timings.
- Run applications in isolated working directories to protect saved outputs.
- Record relevant compiler/dependency versions and source/input/reference identities. Distinguish freshly executed validation from historical evidence and unrun workflows.

## Records and provenance

Keep concise stage summaries, one authoritative phase status and one consolidated review record. Report actual simplification and remaining limitations without duplicating narratives across documents.

Historical evidence and dependency information are in:
- `docs/cleanup/campaign.md`
- `docs/cleanup/deferred-issues.md`
- `docs/cleanup/baseline-manifest.md`

The Stage 00 baseline source and frozen references remain retrievable unchanged at the Phase 3 baseline commit. During an explicitly authorized API or ownership reconstruction, the current harness may migrate with the implementation while preserving benchmark problems, physical inputs, numerical settings, fundamental equations and scientifically meaningful comparisons. Intermediate arrays that exist only because of retired internal APIs may be removed when an explicit old-to-new coverage map identifies their frozen record counts and replacement checks. Keep frozen references byte-for-byte unchanged and preserve the exact comparator, its recorded `rtol=1e-13`, `atol=1e-15`, and baseline source as diagnostic evidence. The explicit `--retire-mapping-intermediates` option may omit only `mapping_perturbation/{dxi,dxilm,df}` rows. A strict comparator failure is not itself a scientific acceptance failure: judge retained quantities using predeclared, unit-appropriate physical criteria, and report exact arithmetic differences separately. The scientific gate must reject unaccounted keys, missing retained quantities, nonfinite values, and out-of-bound physical results. Do not regenerate references or loosen comparator tolerances.

## Current phase status

Phase 3 is accepted at implementation candidate `8418063561161f89e9b15f0c7040d2207855b255`, based on `15582edab8112dbb63e0101be5b26f5135298644`. Phase 4, including the reviewed seven-file correction, has now been accepted by the user; retain its prior review verdicts and limitations as historical evidence in `docs/core-cleanup/phase-4-handoff.md`. The user's current request authorizes Phase 5 stages 5A–5E. The sole writer is Luna; a separate independent Luna reviews each frozen stage; the coordinator records/commits accepted checkpoints locally and advances within the authorized phase. Phase 5A must produce a concrete keep/remove/replace inventory before implementation. Preserve `experimental/`, frozen references, independent checks, datasets and user outputs. At the Phase 5 boundary, complete Sol-high review and stop for human review. Do not push, merge, tag, delete branches, or start further modernization. The user reports that the Phase 4 correction was committed and pushed; the checkout is currently at `7cda0caa723b40bfed865fe1d8951a7f98ee653d` on `refactor/core-cleanup`. Do not claim independent verification of the remote push. The coordinator's local Phase 4 acceptance-record commit will become the Phase 5 baseline and must be recorded after it exists.
