# Phase 3 continuation handoff — 2026-10-09

The reviewed, committed checkpoint is `0912fd62a2ecf1ebab6ab0c2b4fbc05ac9b14c45` on `refactor/core-cleanup`, based on exact Phase 3 baseline `15582edab8112dbb63e0101be5b26f5135298644`. Independent Luna stage reviews passed after their recorded fixes. Sol-high cumulatively **ACCEPTED this partial checkpoint**, subject to the record corrections documented in the [review ledger](phase-3-reviews.md). The exact read-only Sol report is archived at [`evidence/phase3-final/sol-high-cumulative-review.md`](evidence/phase3-final/sol-high-cumulative-review.md), SHA-256 `386cf721e9f03e0bbc47a389415bb660e73180bd1346eecd17a977c71a45203f`. The user has since superseded the earlier stop instruction and authorized Phase 3 continuation through 3F.

This acceptance covers the reviewed Stage 3A evidence/contracts and bounded Stage 3B typed-input, radial-owner, and Phobos parsed-input migration at that commit. It does not accept full Phase 3 completion, paper reproduction, or performance preservation. 3A and 3B retain the bounded outstanding items in the [active plan](phase-3.md); 3C–3F are now authorized and not started at this checkpoint. No Phase 4 work was performed.

The accepted workflow-specific clean_bench_3 residual is approximately `2.08e-12`; the accepted clean_bench_6 range is `4.43e-10`–`7.76e-10`. These observations do not block reconstruction and do not authorize solver-algorithm changes. The captured performance runs are baseline observations, not a replacement-versus-baseline comparison. Some paper-to-example labels and the interactive Phobos `lsmall` remain unmapped; this uncertainty must be documented without blocking reproducible supported-workflow reconstruction. Physical comparison criteria are fixed in the active plan before replacement results are evaluated.

Captured Phobos gravity and heterogeneous workflows retain strict componentwise output mismatches, reaching 74,231 and 216,724 components respectively. Their cause is not established by baseline variability, so no blanket parity claim is made. Focused parser contract checks remain for malformed trailing PREM header/count tokens, Phobos triangular row order, and coefficient finiteness; no open-ended malformed-input campaign is planned.

The full Release build, 16/16 CTest, Stage 00 checks, and captured workflows in `evidence/phase3b2/` predate the final double-only `RadialEarthModel` scalar constraint and OutputZiheng timer-placement correction. The exact committed source was checked by the focused three-target rebuild and 2/2 CTest after those changes. Do not describe the earlier broad results as a full post-correction rerun. The 3B.2 evidence manifest is a preserved reviewed snapshot (SHA-256 `dd6ddbd541fb1b237efc6069f716545b7c88936945dc00096d5502a224362007`); its earlier pre-review hash `cae285945554f8fd113499edeb3ce12cfe4a7730e76fa99218b37540f78be7d7` is retained as historical provenance. The final record-only amendment hashes are listed below.

For scope and provenance, see the [Phase 3 plan](phase-3.md), [review ledger](phase-3-reviews.md), [3B.2 evidence README](evidence/phase3b2/README.md), and [cumulative Sol review](evidence/phase3-final/sol-high-cumulative-review.md). The continuation status and current first-stage candidate are recorded in [`../../implementation_status.md`](../../implementation_status.md).


## Historical partial-checkpoint amendment hashes

These documentation hashes identify the 2026-10-08 closeout amendment only; they are historical and do not identify the current Phase 3 continuation documents, change the reviewed source commit, or alter archived workflow artifacts.

- `docs/core-cleanup/phase-3.md`: `536d4b8526b2dad67bb246dc1a8775a77f6c8cdbdcd5f91a6fb6ab306d6b5175`
- `docs/core-cleanup/phase-3-reviews.md`: `cd283edb5e5640a02b40233c592cecc1749a1f91cdfae9ea939621b19cb2f0d0`
- `implementation_status.md`: `188be615aeb19face6efa0cb2e768458e7921d86a7accd2c06892c6bc0562d7c`
- `docs/core-cleanup/evidence/phase3b2/README.md`: `3697f56ce70c200264344a13bf1e1a8a7c6f8c7dffd2442196d9b2950cbb730a`
- `docs/core-cleanup/evidence/phase3-final/sol-high-cumulative-review.md`: `386cf721e9f03e0bbc47a389415bb660e73180bd1346eecd17a977c71a45203f` (identical to source report)
