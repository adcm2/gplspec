# Phase 3 partial handoff — stop after 3B

The reviewed, committed checkpoint is `0912fd62a2ecf1ebab6ab0c2b4fbc05ac9b14c45` on `refactor/core-cleanup`, based on the exact Phase 3 baseline `15582edab8112dbb63e0101be5b26f5135298644`. Independent Luna stage reviews passed after their recorded fixes. Sol-high cumulatively **ACCEPTED this partial checkpoint**, subject to the record corrections documented in the [review ledger](phase-3-reviews.md). The exact read-only Sol report is archived at [`evidence/phase3-final/sol-high-cumulative-review.md`](evidence/phase3-final/sol-high-cumulative-review.md), SHA-256 `386cf721e9f03e0bbc47a389415bb660e73180bd1346eecd17a977c71a45203f`.

This acceptance covers the reviewed Stage 3A evidence/contracts and the bounded Stage 3B typed-input, radial-owner, and Phobos parsed-input migration at that commit. It does not accept all of 3A or 3B, solver accuracy, full Phase 3 completion, paper reproduction, or performance preservation. Stages 3A and 3B remain incomplete. The user directed the team to stop after 3B; stages 3C–3F are not started, solver-dependent work remains paused, and no Phase 4 work was performed.

The solver accuracy hold remains material. The exact mapped clean_bench_3 run requests `1e-12` but has independent relative residual `2.079844310e-12`. The exact file-backed clean_bench_6 run requests `1e-14` but has independent relative residuals around `4.43e-10` to `7.76e-10`. No algorithm or tolerance was changed. The captured performance runs are baseline observations, not a replacement-versus-baseline gate. Paper-to-example settings and the interactive Phobos `lsmall` remain unmapped.

Captured Phobos gravity and heterogeneous workflows retain strict componentwise output mismatches, reaching 74,231 and 216,724 components respectively. Their cause is not established by baseline variability, so no blanket parity claim is made. Known parser inputs and exact comparisons pass, but exhaustive malformed-input validation is incomplete: optional PREM header/count parsing may stop at a nonnumeric trailing token, and the typed Phobos reader does not validate triangular row order or coefficient finiteness.

The full Release build, 16/16 CTest, Stage 00 checks, and captured workflows in `evidence/phase3b2/` predate the final double-only `RadialEarthModel` scalar constraint and OutputZiheng timer-placement correction. The exact committed source was checked by the focused three-target rebuild and 2/2 CTest after those changes. Do not describe the earlier broad results as a full post-correction rerun. The 3B.2 evidence manifest is a preserved reviewed snapshot (SHA-256 `dd6ddbd541fb1b237efc6069f716545b7c88936945dc00096d5502a224362007`); its earlier pre-review hash `cae285945554f8fd113499edeb3ce12cfe4a7730e76fa99218b37540f78be7d7` is retained as historical provenance. The final record-only amendment hashes are listed below.

For scope and provenance, see the [Phase 3 plan](phase-3.md), [review ledger](phase-3-reviews.md), [3B.2 evidence README](evidence/phase3b2/README.md), and [cumulative Sol review](evidence/phase3-final/sol-high-cumulative-review.md). The final acceptance and stop status is also recorded in [`../../implementation_status.md`](../../implementation_status.md).


## Final record-only amendment hashes

These documentation hashes identify the closeout amendment; they do not change the reviewed source commit or archived workflow artifacts.

- `docs/core-cleanup/phase-3.md`: `536d4b8526b2dad67bb246dc1a8775a77f6c8cdbdcd5f91a6fb6ab306d6b5175`
- `docs/core-cleanup/phase-3-reviews.md`: `cd283edb5e5640a02b40233c592cecc1749a1f91cdfae9ea939621b19cb2f0d0`
- `implementation_status.md`: `188be615aeb19face6efa0cb2e768458e7921d86a7accd2c06892c6bc0562d7c`
- `docs/core-cleanup/evidence/phase3b2/README.md`: `3697f56ce70c200264344a13bf1e1a8a7c6f8c7dffd2442196d9b2950cbb730a`
- `docs/core-cleanup/evidence/phase3-final/sol-high-cumulative-review.md`: `386cf721e9f03e0bbc47a389415bb660e73180bd1346eecd17a977c71a45203f` (identical to source report)
