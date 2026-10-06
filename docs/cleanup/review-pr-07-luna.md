# PR 07 correction pass — Luna high independent review

**Verdict: PASS WITH ONE NON-BLOCKING DOCUMENTATION FINDING.** I found no blocking source or regression-test issue in the fixed correction candidate. Update the current validation status in `docs/cleanup/pr-07-corrections.md` before the next independent Sol review, then freeze and identify a revised candidate.

## Fixed identity and review scope

- Base: `a9409aed3902bf1dd4262ff7cdb62287fd4ea2e1` (`cleanup/base`). The candidate is the complete proposed correction patch in `/tmp/gplspec-pr7-candidate-v1.patch`; its SHA-256 is `e98b9c02af2e1096f01b93a37d4f403ceb565832c6c9865673b2988a06ff037b`.
- Fixed snapshot: `/tmp/gplspec-pr7-candidate-v1`, 268 entries. Its manifest SHA-256 is `1eb1989075ef914f069c39534bfb6b72074d054a1c45feff6845f89a904909d0` (as supplied in the handoff; actual digest verified). All 268 manifest entries match the snapshot bytes. The existing CMake validation build points at `/home/adcm2/Documents/Research/Projects/gplspec`, rather than the `/tmp` snapshot; I compared every manifest entry at that configured source path against the fixed manifest and found zero missing or mismatched files. The candidate exercised by that build therefore matches this reviewed snapshot for all 268 manifested files.
- This was a read-only review. I did not modify the repository, run tests, regenerate references, or mutate Git state.
- Read the repository instructions, the correction handoff and PR description, the prior full-PR Sol review, and the focused regression and production changes. Examined the persisted full-build log and CTest log, snapshot manifest, and Git close-out refs.

## Findings

### [Non-blocking; please correct before Sol review] Correction handoff still reports validation as pending

`docs/cleanup/pr-07-corrections.md:20` says the full build, complete CTest, stage-00 runner, public-header link, and final diff/reference checks “are still in progress.” The same candidate’s `implementation_status.md:309` records them as complete: fresh all-target build, 12/12 CTests, three 20,700-record comparisons with zero changes and byte-identical `MatrixSolution.out`, public-header link, tutorial-fence compilation, clean `git diff --check`, and no reference changes. This stale line makes the handoff internally inconsistent. Amend it to state the completed results and preserve the evidence limits, then freeze a revised snapshot before the cumulative Sol review. It does not undermine the source/test findings below.

## Source and test review

- `ModelDensityOutputRotated` now uses `DensityNorm()` for density serialization. `ReferentialOutputRotated` and `PhysicalOutputRotated` continue to use `PotentialNorm()`. The focused fixture creates distinct scales (8 and 4), so the former incorrect potential-scale multiplication would fail its density expectations.
- The regression constructs a two-layer density jump and exercises a nonunit radial Jacobian. It checks discontinuous markers in all rotated potential writers, inner-side trace selection, first and final endpoints, continuous potential rows, referential and physical density scaling, `rho_ref = J rho_phys`, a nonzero body row immediately inside the surface, and a zero vacuum endpoint. The nonzero-jump and nonunit-J checks prevent those portions of the fixture from passing vacuously.
- The physical-density assertion uses the model’s sampled reference density and Jacobian as its expected value; it verifies output selection/scaling and the documented pointwise relation for this fixture. It does not establish correctness of the full solver or physical density for arbitrary models. The supplied continuous potential array tests writer continuity only; it is not solver-produced evidence.
- The tutorial build target extracts the actual C++ fence from `_tutorials/tutorial5_phobos.md`; the captured generated source has initialized coordinates, declared angle vectors and solution, and created output directories. The persisted build log shows `tutorial5_phobos_compile` built from that extracted source.
- MathJax configuration precedes the asynchronous loader. The relevant source ordering is correct; no site render is claimed.
- The PR description properly distinguishes the pre-cleanup behavior and example changes from cleanup and the new density normalization correction. It does not claim full-PR equivalence to `main` or numerical validation of `OutputZiheng`/`experimental/`.
- Stage-07 close-out claims were checked against Git: `HEAD`, local `cleanup/base`, `origin/cleanup/base`, and the peeled annotated `gplspec-clean-baseline-v1` tag all resolve to `a9409aed3902bf1dd4262ff7cdb62287fd4ea2e1`.

## Existing validation evidence inspected

- `/tmp/gplspec-pr7-final-build.log` reaches 100% and includes all example targets, four public-header smoke objects, `stage02_header_link`, `pr07_rotated_output`, and `tutorial5_phobos_compile`.
- `/tmp/gplspec-pr7-validation-build/Testing/Temporary/LastTest.log` records 12/12 passing CTests, including the new rotated-output regression. The stage-06 fixture still matches its seven frozen output files byte-for-byte.
- The worktree’s `implementation_status.md:309` records the stage-00 runner, reference preservation, `git diff --check`, and link checks. The console transcript for stage-00 is not persisted in the inspected artifacts; I rely on that explicit recorded result rather than claiming independent reproduction.

The stage-06 smooth, unit-scale fixture does not detect the density normalization correction or establish main-to-head behavior at discontinuous interfaces. No frozen reference or tolerance change appears in the proposed correction record. These are evidence boundaries, not blockers for this correction candidate.

## Narrow v2 re-review — documentation finding resolved

**Verdict: PASS.** The sole finding above is resolved in the revised fixed snapshot. `docs/cleanup/pr-07-corrections.md` now records the completed configure, full build, 12/12 CTests, stage-00 result and its transcript limitation, public-header link, and reference/diff checks. This is consistent with the already-recorded completion in `implementation_status.md`. The status entry also preserves that v2 is the new frozen candidate for cumulative Sol review.

- Exact v2 identity: `/tmp/gplspec-pr7-candidate-v2`, 268 files; manifest SHA-256 `075f5cfe1ee5f7064cf9e35042b07c729b6495c732c1a9167d909bf045919194`.
- I recomputed the manifest digest and checked all 268 snapshot file hashes: zero missing entries and zero mismatches. Comparing v1/v2 manifests shows only `docs/cleanup/pr-07-corrections.md` and `implementation_status.md` changed.
- No repository source or tests were edited, and no build or tests were rerun, as requested. The previous source/test review and its evidence remain applicable.
