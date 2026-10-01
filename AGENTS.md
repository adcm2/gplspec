# GPLSpec cleanup campaign instructions

The user has resumed the campaign through stage 05 and the independent Sol review. Implement stages 03–05 sequentially with a separate Luna review and integration after each stage, then obtain the consolidated Sol high review and stop for the human-approved checkpoint. Stages 06–07, GSHTrans modernization, performance work, and pushes remain unauthorized. Preserve the mathematical operations, public interfaces, output formats, and provenance described in `docs/cleanup/campaign.md`.

Before each stage, start from the latest accepted `cleanup/base` commit and use its named `cleanup/NN-*` branch. Keep `main` and `develop` unchanged. Do not push. A stage advances only after its fixed candidate commit passes its independent Luna review. Do not start a stage before the preceding stage receives independent Luna approval and is integrated into cleanup/base. After stage 05 passes Luna review, obtain the independent Sol high review, integrate only if accepted, then stop for the human checkpoint. Do not begin stages 06–07, modernize GSHTrans, or start performance work.

Keep one active implementation writer. Review the candidate in an isolated read-only worktree. Record each code/build change and validation result in `implementation_status.md`. Never regenerate the frozen stage-00 baseline from a later stage.

Implementation uses gpt-6-luna medium; each fixed stage candidate receives an independent gpt-6-luna high review. The cumulative stage-00–05 checkpoint receives gpt-6-sol high after stage 05. Preserve one active writer and keep reviewers isolated and read-only.

Stage-00 baseline command:

```sh
cmake -S . -B /tmp/gplspec-cleanup-baseline \
  -DMY_PROJECT_BUILD_EXAMPLES=OFF \
  -DGPLSPEC_BUILD_BASELINE_HARNESS=ON
cmake --build /tmp/gplspec-cleanup-baseline --target stage00_reference -j2
tests/run_stage00.sh /tmp/gplspec-cleanup-baseline/bin/stage00_reference
```

## Current checkpoint

Stages 00–05 are integrated and recorded by annotated tag `gplspec-cleanup-stage05`. Independent Luna stage reviews and consolidated Sol review are complete (Sol: PASS WITH NON-BLOCKING FINDINGS). Work is stopped awaiting the requested human review and approval. Do not begin any later phase without explicit user authorization.

Automated review addendum completed for `f0fa6920ccc9b53b1d000365fd6ff2f64b5ed700`: Sol high **PASS WITH NON-BLOCKING FINDINGS**; see `docs/cleanup/review-stage-05-sol-addendum.md`. Human review remains OPEN. The existing stage-05 tag stays at the earlier `83e1a9fb3021f40e6384278432f502baba40bcde`; do not move it. Do not commit unless explicitly requested and confirmed with the user. No later stage is authorized.
