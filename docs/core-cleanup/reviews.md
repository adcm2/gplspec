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

**Cumulative evidence correction:** The preceding implementation-response wording records the Stage 1B plan at the time, not a completed direct transform-parity test. Stage 1C provides indirect frozen whole-workflow behavior coverage; it does not directly call GSHTrans forward/inverse methods or assert coefficient range/indexing, spin, normalization, or sign/conjugation. Add and review that focused direct parity gate before introducing any Phase 3 transform adapter. No Phase 1 numerical test is added for this correction.


## Stage 1C — initial Luna review

**Verdict: BLOCKED.**

The reviewer found that the standalone public-umbrella compile checks covered `Core`, `GeneralModels`, `SimpleModels`, and `All`, but omitted the supported `Test` and `Timer` umbrellas.

**Implementation response:** Added `tests/stage02_test_smoke.cpp` and `tests/stage02_timer_smoke.cpp` and included both umbrellas in the standalone smoke-target loop in `tests/CMakeLists.txt`. Both compile as separate translation units in the fresh Release build; no header or production implementation repair was needed. Added both sources and all build/CI inputs to the Stage 1C content manifest, which compares candidate hashes with the exact Stage 1C parent and records inputs/references separately. Updated campaign and implementation status accordingly.

**Disposition: resolved pending independent re-review.** The prior full build, 14/14 CTests and frozen Stage 00 comparison remain valid because the response only adds two compile-only object targets; both new targets were force-recompiled successfully. No CTest count or numerical result changed. No PASS or Stage 1C acceptance is claimed yet.

**Stage 1C independent Luna re-review: PASS.** All six public umbrellas have standalone coverage; follow-up compilation and link evidence resolves the blocker. The reviewer independently verified all 96 base/candidate manifest hashes with zero mismatches. No remaining findings; production, inputs, references, dependencies and experiments remain unchanged. Coordinator acceptance record follows the frozen reviewed candidate.


## Stage 1D — initial Luna review

**Verdict: BLOCKED.**

The initial build-time evidence timed only the `phase1_performance` target after cleaning the 70 configured targets. The review requested timing the full all-target clean Release build. It also requested explicit methodology limits: the first public solver call occurs after setup and six matrix-vector actions; full public-call timings include output capture/copy overhead; independent residual diagnostics run once per process outside the timed calls rather than for every timed solve.

**Implementation response:** Repeated the clean-first build for all 70 configured targets with `-j2` and controlled thread variables. It succeeded in 3:41.47 wall time with peak RSS 980,376 KiB; the complete build and `/usr/bin/time -v` logs are saved as `stage1d-clean-full-build.log` and `stage1d-clean-full-time.txt`. Updated the evidence record to state that this excludes configure/dependency fetch, and to distinguish the first public call from a cold end-to-end solve and describe capture overhead and per-process residual diagnostics. The timing harness and runtime batch were not changed or rerun.

**Disposition: resolved pending independent re-review.** No production algorithms, tolerance, references, inputs, or dependency pins changed. Stage 1E remains gated on Stage 1D acceptance.

**Stage 1D rereview: substantive PASS.** The sole documentation finding was a stale sentence in `phase-1.md` that still described target-only build timing. It now reports the measured clean-first all-target Release build (70 steps, 3:41.47, 980,376 KiB peak RSS); no rebuild was needed.

**Stage 1D final independent Luna re-review: PASS.** Full clean build evidence and the aligned summary resolve all findings. No further runtime or numerical rerun was needed. Coordinator acceptance follows the reviewed candidate.


## Stage 1E — Phase 1 decision package

The package in `docs/core-cleanup/phase-1.md` consolidates the accepted Stage 1A–1D inventory, architecture proposal, correctness baseline and bounded performance evidence, with explicit limitations, user decisions and later-phase gates. Checkpoint commits are recorded there: base `ac62330e1a43d06f2d4d2c64492817fd7cc6c740`; 1A `31e75c18f4b52ad1855ae425249a8c2b264ac5a9`; 1B `b68bb219e94eac50ecc145c9eb6388ac5fb61f57`; 1C `22306e257c64240f4b0536abdf1ea49f0e2cb326`; 1D `a491da647efa2baa827224ea53d21a616701adf9`.

**Stage 1E independent Luna review: PASS** on the prior frozen package snapshot. **Cumulative Sol-high review: BLOCKED** on two documentation/evidence reconciliations; no production regression was found and no runtime rerun was requested.

1. The tracked PREM500/S40RTS fixture constructs a model and checks a nonzero complex operator action; it does not solve that fixture. Earlier zero-iteration wording had no saved provenance. The independent 3-iteration/residual evidence belongs to the separate mapped model.
2. Direct GSHTrans transform/index parity was promised in the plan and described as checked in the Stage 1B implementation record, but no direct parity test exists. Stage 1C provides indirect frozen whole-workflow coverage only. Make direct checks of range/indexing, spin, normalization, sign/conjugation and grid conventions a required gate before any Phase 3 transform adapter.

The review also requested precision that all 30 timed public calls report 3 iterations while independent residuals were recomputed in five per-process diagnostics; the clean-build 980,376 KiB figure is GNU `/usr/bin/time -v` maximum RSS, not aggregate concurrent `-j2` worker memory; and the 96-row content manifest is a Stage 1C snapshot, not a final Phase 1 manifest. These clarifications are included in the evidence summary. The accepted prior Luna snapshot and this documentation-only response are distinct; affected Luna review and Sol-high re-review remain pending. No Phase 1 PASS is claimed. Phase 2 remains unauthorized pending explicit user approval.

### Final affected reviews and phase acceptance

**Independent Luna affected re-review: PASS, no remaining findings. Cumulative Sol-high affected re-review: PASS, no remaining blocking findings.** Sol verified both evidence corrections, the iteration/residual distinction, and the non-aggregate build-memory label. The cumulative source, test, CI, raw-data and scope assessment remains applicable; no numerical rerun was required. Hosted CI remains unrun, direct transform parity is a pre-adapter gate, and other scientific/performance limitations remain explicit in the decision package.

The fixed Sol-reviewed candidate was HEAD `a491da647efa2baa827224ea53d21a616701adf9` plus these four documentation hashes:

| Path | SHA-256 |
|---|---|
| `docs/core-cleanup/phase-1.md` | `a84cf5867e7ea448c9f200bc57801900562c0117f0e7c537e8436b90ef40c2ca` |
| `docs/core-cleanup/reviews.md` | `fd544e3b36fd451304d03310001669fb2337997c02bec85edf27a0af9a97c4a8` |
| `implementation_status.md` | `c6fbd3640aec68e6a928960c6fcc4ecfcdb611abe4ccb532ba4e4ad082bbdda6` |
| `docs/core-cleanup/evidence/README.md` | `55ad25ed0a9dfdaceb60a45bd020d8f1b6e899fd814ab507f522ceea4a458b87` |

The coordinator verified those hashes before updating only the campaign status, this review record and `implementation_status.md` to record final acceptance. This small post-review record delta is separate from the reviewed snapshot and does not alter tests, source, inputs or numerical evidence. Stage 1E is accepted for local commit; the Phase 1 human-review checkpoint remains open. No Phase 2 work, push, merge or tag change is authorized.
