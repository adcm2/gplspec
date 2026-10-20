# Phase 2 review record

## Stage 2A — scoped map and decisions

**Verdict: PASS, no blocking findings.** Independent Luna review of the frozen Stage 2A diff:

> Stage 2A review: PASS, no blocking findings. I reviewed the frozen diff against the Phase 2 authorization, workflow protections, and phase boundaries. git diff --check passes, and all 132 protected-file hashes match. The plan’s build and documentation candidates are grounded in the repository. The nonblocking follow-up for Stages 2C/2D is to make the per-file documentation migration and Phase 3 retirement decisions concrete, as the plan already requires. No changes made.

Disposition: the nonblocking 2C/2D follow-up is retained as a requirement in the Phase 2 plan. Stage 2B may proceed. The review made no repository changes.


## Stage 2B — build setup simplification

**Verdict: PASS, no blocking findings.** Independent Luna review accepted the build changes and validation. The reviewer noted one nonblocking hash-accounting clarification: the Phase 2A manifest has 132 entries, three of which are the root/examples/tests CMake files intentionally changed in Stage 2B; the remaining 129 protected files match their starting hashes. Saved figures are also unchanged.

Disposition: recorded the count distinction in the Phase 2 plan, implementation status, and validation evidence. Stage 2C may proceed.


## Stage 2C — initial independent Luna review

**Verdict: BLOCKED, three factual documentation groups.** The reviewer found that the `clean_bench_3` equivalent-sphere density was inverted; `OutputZiheng` also reads `modeldata/prem.200.no` for a separate PREM integral product; and the Phobos plot table included inactive or overwritten alternatives as if they were active file reads for `PhobosSliceDensity.py`, `PhobosSlicePlot.py`, `PhobosSlicePlot2.py`, and `PowerSpectrum.py`.

**Response:** corrected the density to `1/(1.2^3)`; documented the separate PREM input and `OneDimensionalPlanetFromFile(path, 5, 2, 0.001, 1.2)` call; changed the plot consumer rows to list only active final file paths, including the final `path6` value in `PhobosSlicePlot2.py`. This is documentation-only. Source/build/plot consumers are unchanged. The candidate is awaiting independent Luna re-review.


## Stage 2C — independent Luna re-review

**Verdict: PASS, no blocking findings.** All three previously reported factual groups were corrected and source-verified: the `clean_bench_3` equivalent-sphere density, `OutputZiheng`'s separate `prem.200.no` input and constructor call, and the active final file paths read by the four reviewed Phobos plot scripts. No new blockers were reported. The changes remain documentation-only; Stage 2D may proceed.


## Stage 2D — initial independent Luna review

**Verdict: BLOCKED on three recommendation/staging issues.** (1) The handoff mandated a GPLSpec-owned harmonic-transform backend without evidence that this is simpler or better; (2) it scheduled `testremove/` deletion in 3A before the proposed 3C solver exists for comparison; and (3) perturbation migration was assigned to 3D while its operations were assigned removal to 3E.

**Response:** retain the pinned GSHTrans backend and make only a call-ownership consolidation the selected 3D path; a new backend is an unselected option needing evidence and fresh human approval. Preserve all six `testremove/` programs through 3C comparison and allow removal only after that comparison and any unique check is promoted. Keep existing perturbation mathematics during 3C/3D, then migrate and remove it with its temporary internal solver adapter in 3E. No source/API implementation was changed. The corrected documentation candidate awaits independent re-review.


## Stage 2D — independent Luna re-review

**Verdict: PASS.** The corrected handoff retains pinned GSHTrans, delays `testremove/` decisions until solver comparison in 3C, and aligns perturbation retirement to 3E. No blocking findings remained. The cumulative Sol-high review that followed identified two additional handoff gaps; those are recorded below.

## Stage 2D — initial cumulative Sol-high review

**Verdict: BLOCKED on two Phase 3 handoff gaps.** The solver/RHS proposal omitted explicit `FindForce` 1-D/3-D mathematical and layout contracts, source-unit/quadrature/index/origin/interface behavior, the replacement owner, named callers/checks, and removal points. `FindBoundaryPerturbationForce` also lacked a 3E contract, callers, and zero-map/finite-difference checks. Separately, the 3C solver retirement gate compared parity and residuals but omitted controlled performance evidence before deleting the old path. The reviewer required matched Release builds on the same host, toolchain, dependency pins, and thread settings; representative mapped and file-tomography fixtures; repeated raw setup/source/operator/solve timing; peak RSS; iteration counts; independent residuals; and investigation of material regressions without invented thresholds given Stage 1D noise.

**Response:** added a source/RHS contract row describing both `FindForce` overloads and their model-normalized density/G, quadrature, harmonic layout, origin/interface and zero-slot conventions; assigned 3C solver ownership and 3E boundary-perturbation ownership; named tracked callers and current checks; and required independent monopole/quadrature, low-degree harmonic, zero-map, and mapped-source finite-difference checks before removal. Added a pre-retirement 3C same-host matched Release comparison with repeated raw setup/source/operator/solve timings, peak RSS, iterations, and residuals; require investigation of regressions, and assign no hard performance threshold. The 2D public/private wording now states plainly that the intended supported surface is public and low-level solver/implementation headers are intended to be private. No source/API implementation changed; corrected handoff awaits renewed independent reviews.


## Stage 2D — affected independent Luna re-review

**Verdict: PASS on the force/RHS contracts and 3C performance gate; one wording clarification requested.** Source checks confirmed the newly specified 1-D/3-D `FindForce` contracts, `FindBoundaryPerturbationForce` ownership/checks, and matched Release pre-retirement comparison. The remaining wording issue was that the solver row called the replacement operator/preconditioner private while saying to preserve their current public status.

**Response:** clarified that the replacement matrix-free operator and preconditioner become private implementation details after tracked callers migrate, with no forwarding public API. This record-only note and wording correction await cumulative Sol re-review.

## Final cumulative Sol-high re-review and acceptance

**Verdict: PASS, no remaining blocking findings.** Sol verified the explicit 1-D/3-D source/RHS and boundary-perturbation contracts, callers, independent checks and removal stages, the matched Release performance gate before solver retirement, and the clarified public/private boundary. The Stage 2B build/regression and Stage 2C focused evidence remain applicable; no numerical rerun was needed for the documentation-only corrections. All independent Luna stage reviews have passed.

The reviewed candidate was HEAD `d4d250a2c4a5861fbffad35ca23966a3f50658f8` plus 21 changed paths (14 modified, three deleted, four added), including the starting rewrite clarification. The complete manifest `/tmp/gplspec-phase2-sol-candidate-manifest-rereview.txt` has SHA-256 `d6e9eef79b8d62cfb75582a9da4096513d4fa8a9462bf8f2c70fc2e9908d7b34`; the original Sol report and re-review addendum remain at `/tmp/gplspec-phase2-sol-review.md` and `/tmp/gplspec-phase2-sol-review-addendum.md`. The coordinator verified every manifest content hash before updating only this record, `phase-2.md` and `implementation_status.md` with the completed verdict. That post-review delta records acceptance, without changing reviewed implementation or evidence.

Remaining limits: declared CMake 3.20 was tested with 3.28.3 only; Jekyll/site rendering and hosted CI were not run; full paper workflows were mapped but not executed; the prior scientific and performance limitations remain. The Phase 2 diff is uncommitted and ready for human review. Phase 3 is not authorized.


## Post-review formatting note

After cumulative Sol-high review accepted Phase 2, trailing spaces/tabs were removed from six lines in the new `evidence/phase2b-validation.log`. No logged commands, outputs, review findings, numerical/source content, or validation results changed. This is documentation formatting only.
