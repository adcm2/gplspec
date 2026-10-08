# Phase 3 review record

## Pre-baseline Phase 2 acceptance

The Phase 2 candidate passed its independent Luna stage reviews and cumulative Sol-high re-review after the handoff corrections. The final Phase 2 review record and its findings are in `phase-2-reviews.md`. The coordinator matched 20 manifest paths and separately reviewed the post-handoff `AGENTS.md` revision, preserving its authorization/status edits in the baseline. The accepted Phase 2 tree was committed as `15582edab8112dbb63e0101be5b26f5135298644` on `refactor/core-cleanup`. This records the pre-baseline Luna PASS; it is not a review of Phase 3 changes.

## Stage 3A.1

Accepted candidate: commit `bf844159ffd82383ee1c5abfae1374a785f46a2a`, based on `15582edab8112dbb63e0101be5b26f5135298644`. Scope is `docs/core-cleanup/phase-3.md`, `tests/stage3a_harmonic_contract.cpp`, `tests/CMakeLists.txt`, `tests/external_consumer/main.cpp`, the Stage 3A evidence directory, and the matching `implementation_status.md` entry. No production numerical implementation changed. Targeted Release validation passed: final CTest transcript `evidence/phase3a/candidate-ctest.log` (2/2, 22.43 s; SHA-256 `b3551a1032ef4e26ee43c30bc58d660b04f37bfdba6d9b9ccc0b3084ac1d3f4c`) covers the harmonic contract test and external consumer with both translation units using GPLSpec. Independent Luna review: **PASS**, no blocking science or implementation findings. Stage 3A.1 is accepted and committed; overall Stage 3A remains incomplete.


## Stage 3A.2.1 — file-backed source and solve

Candidate: test-only expansion of `tests/stage1c_file_tomography_solver.cpp` plus the evidence summary/log at `evidence/phase3a/file-backed-solve-quadrature.md` and matching status/plan updates. Preserves the earlier file-backed complex-action and synthetic mapped solve checks. Fresh Release build passed; focused CTest passed 1/1. The actual PREM500/S40RTS solve reported 0 CG iterations from the preconditioned guess, with independently recomputed relative residual `5.06196e-12`; direct analytic `l=0,1,2; m=0` source quadrature differed from `FindForce` by at most `3.55271e-15`. No production numerical implementation changed.

**Initial independent Luna review: BLOCKED.** The candidate's record claimed the prior deterministic complex file-backed action regression was retained, but that input/action/nonzero-finite assertion was absent from the candidate source. **Response:** restored the exact original complex vector formula and action check (`Complex(0.125*(i+1), -0.0625*(i+2))`); rebuilt and reran the direct check and focused CTest. Final source SHA-256 is `dbfa373d8a59ee2e592a23325bc5cbca50110c5c983a5bcbae6cfae792927a41`; CTest passes 1/1 in 0.27 s. Corrected transcript hash is recorded in `file-backed-solve-quadrature.md`. Independent Luna re-review: **PASS**, no blocking findings. Stage 3A.2.1 is accepted and awaits coordinator commit; overall Stage 3A remains incomplete.
