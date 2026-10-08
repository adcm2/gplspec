# Phase 1 independent review record

This record contains reviews of the current Phase 1 campaign only. Earlier `docs/cleanup/review-stage-*` files concern a prior campaign.

## Stage 1A — initial Luna review

**Verdict: BLOCKED.**

The reviewer found that the first inventory did not state the actual radial mesh invariant: each input-model layer boundary is explicitly an element node, and the mesh may also contain additional element boundaries within each layer. Evidence cited by the review is `Radial_Tools::node_info` around lines 134–155 and 179–181. The review also requested concise examples mapping named public routines/classes to tracked caller paths.

**Implementation response:** `docs/core-cleanup/phase-1.md` now records the invariant and its source path, and lists representative callers for `Density3D` constructors, `TestTools::EarthModel`, `TestTools::PhobosRead`, and `Gravity_Tools` solvers. The AGENTS stage gate is now phrased for every stage rather than only Stage 1A.

**Disposition: resolved pending independent re-review.** This response does not claim PASS or Stage 1A acceptance. No production code, tests, numerical references, datasets, outputs, or `experimental/` files were changed.

**Independent Luna re-review: PASS.** The reviewer checked the explicit mesh invariant and named caller mappings against source; the blocker is resolved with no remaining findings. `git diff --check` passed. The coordinator accepted Stage 1A and recorded local commit `31e75c18f4b52ad1855ae425249a8c2b264ac5a9`.


## Stage 1B — initial Luna review

**Verdict: BLOCKED.**

The reviewer found that Stage 1C's solver plan required “meaningful nonzero iterations” without a measurable minimum. The review requires one representative solve with at least two iterations and an independently computed residual plus solver-reported error, without loosening tolerance to force the result. It also asked that the proposed transform seam not imply implementing an adapter in Stage 1C.

**Implementation response:** `docs/core-cleanup/phase-1.md` now requires at least two iterations at the existing requested tolerance and records both independent residuals and solver-reported error. It clarifies that Stage 1C checks parity at the existing GSHTrans call/index boundary only; facade implementation and modernization remain later work.

**Disposition: resolved pending independent re-review.** No production code, tests, references, dependencies or numerical settings were changed.

**Stage 1B independent Luna re-review: PASS.** The reviewer verified the minimum-two-iteration criterion, independent residual contract, and existing-call transform parity scope. No remaining findings. The coordinator acceptance update follows the reviewed documentation snapshot.
