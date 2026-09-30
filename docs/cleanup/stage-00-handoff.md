# Stage 00 handoff — reproducible baseline

- Starting source commit: `4ef3a66c62d64408d99989dd51c3ccbdc46d0b0c` (`cleanup/base`).
- Candidate: the fixed head of `cleanup/00-baseline`; the coordinator will supply its commit SHA with the independent review request.
- Changes: SHA-pinned CMake dependency closure; exact system Eigen 3.4.0 selection; the isolated interface linkage repair recorded in `implementation_status.md`; opt-in deterministic five-fixture baseline harness and frozen raw numerical/output artifacts.
- Validation: clean SHA-only out-of-source configure/build and `tests/run_stage00.sh`; 20,700 numerical records compared twice to the original-source reference and to each other, zero changed components; representative `MatrixSolution.out` byte-identical. Full commands and provenance are in `baseline-manifest.md`.
- Numerical differences: none. The reference and candidate captures are byte-identical on the recorded stack.
- Deferred observations: zero-iteration solver counts are recorded in `deferred-issues.md`; no suspected bug was silently corrected.
- Independent Luna review: pending.
