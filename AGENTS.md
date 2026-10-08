# GPLSpec cleanup campaign instructions

## Goal and scope

Build the simplest readable, maintainable implementation of GPLSpec’s scientific capabilities while preserving correctness, reproducibility and performance.

Treat the existing implementation as a scientific and regression starting point, not a permanent architecture. Replace APIs, types, files and implementations when this simplifies the code. Preserve useful mathematical explanations. Avoid unnecessary abstractions, compatibility layers and duplicate implementations.

Follow the latest explicit user authorization and the active phase plan in `docs/core-cleanup/`. Historical instructions do not override newer authorization. Phase 2 is accepted, and Phase 3 is authorized under the user's pasted request. The coordinator records the exact Phase 3 baseline before implementation begins.

## Execution and review

- Use one active writer: Luna implements; a separate independent Luna reviews each frozen stage.
- Resolve blocking findings, commit accepted stages locally, and proceed automatically within the authorized phase. No separate permission is needed for those stage commits.
- Keep every accepted stage working and independently reviewable.
- At each phase boundary, obtain cumulative Sol-high review, resolve blockers, then stop for human review with the complete cumulative patch and validation evidence.
- Use additional Sol reviews only for concrete mathematical or architectural uncertainties.
- Do not silently substitute models if the requested model is unavailable.
- Do not push, merge, change tags, delete branches or rewrite history without explicit authorization.
- Preserve unrelated user changes and do not stage them implicitly.

## Scientific and simplification requirements

- Preserve supported mathematics, numerical accuracy, performance and scientific workflows. Baseline agreement alone does not prove scientific correctness; document known discrepancies.
- Preserve clean_bench and Phobos paper-verification inputs, parameters and meaningful comparisons. Their source code and APIs may be rewritten.
- Before replacing a component, state its mathematical contract and establish relevant regression and independent checks.
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
- Capture baseline scientific and performance results before replacement; preserve the baseline in an isolated checkout.
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

The Stage 00 baseline source and frozen references remain retrievable unchanged at the Phase 3 baseline commit. The current tracked Stage 00 harness may receive only an API migration required by an explicitly authorized Phase 3 public-API retirement; preserve its fixtures, mathematics, requested tolerances, output records and reference comparator exactly, and compare it against the immutable baseline checkout. Do not use API migration to alter scientific behavior, regenerate references or loosen tolerances. Record newly discovered scientific issues separately and pause dependent changes for human review.
