# Stage 06 human diff review

**Stage 06 human-approved.** The user approved the reviewed diff, including `ForwardAndExpandRealScalarPair` unchanged, and explicitly authorized committing and pushing `cleanup/06-output`. Stage 06 is not integrated into `cleanup/base`.

## Starting state

Repository relocated to `/home/adcm2/Documents/Research/Projects/gplspec`. Starting HEAD is `3045834c55c2906c1d6f60b660feb4e922618cae`; the starting worktree was clean and `cleanup/base` matched local `origin/cleanup/base`. This commit contains only stage-05 closure documentation beyond independently reviewed `f0fa6920ccc9b53b1d000365fd6ff2f64b5ed700`. At the original diff handoff, no reset, commit, integration, tag or push had been performed; later commit/push authorization is recorded below. The stage-05 tag remains at its earlier `83e1a9fb3021f40e6384278432f502baba40bcde` target.

## Reading order and scope

1. `gplspec/src/GeneralModels/Density_Model_Return.h`: shared output angle derivation for three rotated writers; shared Wigner matrix construction for those plus `RotateSliceToEquator`; paired scalar mapping/density transformation and coefficient expansion. The slice method retains its different angle sign and one-sided epsilon checks. Coefficient accumulation, normalization, allocation, physical/referential selection and file writers remain caller-owned.
2. `gplspec/src/GeneralModels/Gravity_Tools.h`: only `Error:` to `Requested tolerance:`. Same tolerance value, output location and solver behavior. This is the sole intentional console change.
3. `tests/stage06_output_helpers.cpp`: exact original-expression comparisons and transform ordering. `tests/stage06_output_reference.cpp`, `tests/run_stage06_output_reference.py`, and `tests/reference/stage06-output-baseline.sha256`: public-method output probe and seven-file baseline checksum check. `tests/CMakeLists.txt` registers both new checks.
4. `review-stage-06.md`, `campaign.md`, and root `implementation_status.md`: review, status and provenance. This handoff is a new documentation file.

## Validation and provenance

Implementation used Luna medium; a separate Luna high reviewer returned **PASS**, no blocking code findings. The reviewer independently verified all 259 files against snapshot manifest SHA-256 `820785a30df90634512baecd6a7268853cc9b5759beb53097e6045b3cb7af951`, reran all 11 CTests, and compared all seven files from separate original/candidate runs byte-for-byte.

The original output reference uses an archive of `3045834c55c2906c1d6f60b660feb4e922618cae` at `/tmp/gplspec-stage06-baseline/source`, built under `/tmp/gplspec-stage06-baseline/build`. Both builds use identical probe source SHA-256 `843e8f70eec0f34c3d86da7d5064ca3ff4375979d6faa39b4960302efafa0425`. The checksum manifest is generated only from that baseline. Existing frozen references are untouched.

```sh
cmake --build /tmp/gplspec-stage06-build -j2
ctest --test-dir /tmp/gplspec-stage06-build --output-on-failure
bash tests/run_stage00.sh /tmp/gplspec-stage06-build/bin/stage00_reference
cmake --build /tmp/gplspec-stage06-build --target clean_bench_1 phobos_gravity phobos_heterogeneous -j2
```

All checks passed: 11 CTests; public umbrella/header/link checks; three examples; three 20,700-record frozen comparisons with zero differences and byte-identical representative output. The seven new file comparisons cover rotated referential/physical potential, both model-density modes, slice coefficients serialized at 17 digits, and unrotated referential/physical potential on a five-element mapped model. Helper tests exercise four angle pairs and degrees 0, 1, 3 and 5. Solver stdout has four approved label substitutions; iteration counts, tolerance values and remaining text match after allowing the two elapsed-time measurements to vary.

Coverage remains finite and tied to the recorded pinned stack. The pre-existing angle/sign and density-normalization quirks are preserved, not repaired. The experimental ellipticity sources are unchanged and remain excluded from production builds. Review requested only documentation reconciliation, now applied; no implementation correction was needed.

## Historical pre-commit inspection commands

```sh
cd /home/adcm2/Documents/Research/Projects/gplspec
git status --short
git diff --stat HEAD
git diff HEAD -- gplspec tests docs/cleanup implementation_status.md
git ls-files --others --exclude-standard
# Complete patch, including new files:
less build/cleanup-review/stage06-working-tree.patch
```

`git diff HEAD` excludes untracked additions; the complete patch includes them. To inspect a new file directly, use `git diff --no-index /dev/null PATH` (exit status 1 means differences were displayed). No files need staging to review this candidate.

## Human approval and small post-review correction — 2026-10-02

The user approved the final stage-06 diff and the paired-transform helper unchanged, and authorized a commit plus non-force push of `cleanup/06-output` with no additional approval. Historical independent review in `review-stage-06.md` remains unchanged and describes its original snapshot.

Luna medium supplied the assertion clarification, applied by the coordinator after the agent's shell failed to start. Only the message and explanatory comment in `tests/stage06_output_helpers.cpp` changed: `angle fixture does not distinguish original slice and output conventions`. The expression is unchanged. The comment explains the fixture distinction and points to the separate baseline test of production `RotateSliceToEquator` output. This correction occurred after independent review and is not attributed to that reviewed snapshot.

The existing helper target is rebuilt and its focused CTest rerun; the final result is recorded in `implementation_status.md`. Earlier full validation remains applicable to unchanged production code. No reference outputs or checksums were regenerated. The human-approved changes may be committed and pushed only on `cleanup/06-output`; no integration into cleanup/base, tag changes, stage 07 or GSHTrans modernization is authorized.
