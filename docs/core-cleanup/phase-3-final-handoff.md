# GPLSpec Phase 3 final handoff

**Phase 3 is complete and accepted by independent Luna and cumulative Sol-high review. The human review checkpoint is OPEN. No Phase 4 work has begun or is authorized.**

## Baseline and reviewed candidate

- Scientific baseline: `15582edab8112dbb63e0101be5b26f5135298644`.
- Full Phase 3 implementation candidate reviewed by Sol-high: `8418063561161f89e9b15f0c7040d2207855b255`.
- Independent Luna review of 3F: **PASS**, no scientific, numerical, or API blockers.
- Cumulative Sol-high review of baseline..candidate: **PASS**, no blocking findings or required fixes. The verbatim report is [archived here](evidence/phase3f/sol-high-cumulative-review.md). It accepts Phase 3 for this human checkpoint; it does not authorize Phase 4 or establish every paper figure or experimental compatibility.
- This handoff and Sol report are a later documentation-only close-out. `8418063` identifies the reviewed implementation candidate; it is not a claim that it remains repository `HEAD` after the close-out commit.
- Human review and the user's decision whether to authorize any next phase remain outstanding.

## Reconstruction delivered

GPLSpec’s scientific capabilities and fundamental solver algorithms were retained while its internal organization and APIs were simplified. `RadialMesh` is the canonical radial geometry owner, including shared global endpoints and element-local material traces. Harmonic transformations use one borrowed const grid accessor and cleared operation-local workspaces. Source assembly and the spherical/mapped weak operator actions are explicit mathematical operations; Eigen product adapters are thin. The solve result owns a flat complex solution and ordered per-solve diagnostics. Named potential, mapping, classical perturbation, and sensitivity operations return this result, and `ToHarmonicField` checks status, vector length, and finiteness at output boundaries.

The existing CG, BiCGSTAB, preconditioner, initial-guess, ordering, and caller-tolerance behavior remain. Redundant node/material arrays and grid aliases, old potential/sensitivity wrappers, temporary perturbation getters, inactive MappingTools helpers, and duplicate integral/solver prototypes were retired after active callers migrated. The old `Spherical_Integrator.h` algorithms were removed; only its `Radial_Node`/`All_Node` forwarding utilities remain for the named protected experimental include. That experimental tree remains unchanged, unbuilt, and unvalidated. clean_bench and Phobos problem definitions, parameters, and scientific output meaning were preserved; current serialization was retained where practical.

## Validation and performance

- Fresh disconnected Release configure and complete build passed; all targets built. CTest passed **25/25**; standalone header-link and external two-TU consumer checks passed.
- The scientific-reference gate passed Stage 00 twice (18,648/20,700 records retained), Stage 05 1-D (1,178/1,178), and Stage 05 file perturbation (1,766/3,818; only 2,052 `dxi`/`dxilm`/`df` intermediate rows retired). Its six positive/negative contract checks passed. Frozen references were not regenerated, and all 121 protected paths remained unchanged. Exact/bytewise diagnostics remain unchanged and can report small differences independently of scientific acceptance.
- Ten original-settings workflows exited 0 and produced all 40 expected output files. **39/40** meet the declared scientific criteria. The Phobos physical-radius coordinate difference passes the existing scale-aware geometry formula. The sole comparison exception is clean_bench_10 classical `perturbed/MatrixSolution.out`, maximum difference `461813.485791 m^2/s^2`; parity fails and it is not counted as equivalent. A bounded baseline-only grid-view ownership probe attributed the source discrepancy to a temporary borrowed-view lifetime in the old baseline, corrected by the accepted const-reference grid owner.
- Matched public API medians (three one-thread processes, median of process medians): mapped 13.893254 ms baseline and 14.026254 ms candidate (+0.96%); file-backed 333.042666 ms and 342.357401 ms (+2.80%). Median RSS changed from 16,784 to 16,792 KiB mapped and 416,476 to 393,060 KiB file-backed. The baseline timer included old nested output conversion; the candidate timer covers `SolvePotential` and its independent residual but stops before explicit `ToHarmonicField`. Candidate residuals were `2.69037115e-12` mapped and `6.39413501e-10` file-backed. Existing solver algorithms were not altered.

Detailed commands, manifests, raw performance samples, diagnostics, and workflow comparisons are in [Phase 3F evidence](evidence/phase3f/README.md) and its linked evidence files. The protected manifest confirms all 121 protected paths are unchanged. The source-freeze record notes one post-validation explanation-only comment update; the executable Python AST is unchanged.

## Numerical corrections and limits

The radial homogeneous-sphere exterior formula now uses the mathematically required `R^3/r` scaling; the previous intermediate-radius branch used `R/r`. The confirmed mass normalization constant was corrected by about `1.43e-11` relative. The classical boundary-source prefactor uses full-precision pi; the previous literal differed by about `2.86e-11` relatively. These are bounded, independently checked formula/normalization corrections, not solver changes. The original baseline’s high-degree clean_bench_10 source leakage is separately attributed and is not hidden by tolerance changes.

The user-accepted clean_bench_3 and clean_bench_6 baseline residual observations remain workflow-specific, not global accuracy thresholds. Candidate residuals are reported above. The supported workflow settings are documented, but ambiguous paper-to-program mappings and interactive Phobos `lsmall` figure interpretation remain unverified. Interactive clean_bench_5/7/8 and the high-degree interactive Phobos perturbation workflow were not among the ten final runs; no paper claim is made for them.

## Proposed Phase 4 scope — not authorized

If the user later authorizes Phase 4, a focused proposal is output/result and plotting migration: define clear named quantities, units, coordinates, and interface-trace semantics; migrate producers and output-processing/plotting consumers together; keep reproducible benchmark and Phobos capabilities and documentation; and retire obsolete output paths only after the scientific comparisons pass. Serialization bytes and filenames are not acceptance requirements under the current scientific reconstruction principle. This is a proposal only; Phase 4 is not authorized. Do not start it before the human review checkpoint is resolved.

## Complete patch for human review

After the documentation-only close-out is committed, review the complete Phase 3 patch from the exact scientific baseline with:

```sh
git diff --no-ext-diff --find-renames 15582edab8112dbb63e0101be5b26f5135298644 HEAD
git diff --stat 15582edab8112dbb63e0101be5b26f5135298644 HEAD
```

The final `HEAD` is the close-out commit, while `8418063561161f89e9b15f0c7040d2207855b255` remains the full implementation candidate cumulatively reviewed by Sol-high.
