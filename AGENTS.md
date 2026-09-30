# GPLSpec cleanup campaign instructions

This campaign is limited to behaviour-preserving cleanup stages 00–05. Preserve the mathematical operations, public interfaces, output formats, and provenance described in `docs/cleanup/campaign.md`.

Before each stage, start from the latest accepted `cleanup/base` commit and use its named `cleanup/NN-*` branch. Keep `main` and `develop` unchanged. Do not push. A stage advances only after its fixed candidate commit passes its independent Luna review. Do not start stage 01 before stage 00 review acceptance; stop after the cumulative stage-05 Sol checkpoint. Stages 06–07 and GSHTrans modernization are out of scope.

Keep one active implementation writer. Review the candidate in an isolated read-only worktree. Record each code/build change and validation result in `implementation_status.md`. Never regenerate the frozen stage-00 baseline from a later stage.

Stage-00 baseline command:

```sh
cmake -S . -B /tmp/gplspec-cleanup-baseline \
  -DMY_PROJECT_BUILD_EXAMPLES=OFF \
  -DGPLSPEC_BUILD_BASELINE_HARNESS=ON
cmake --build /tmp/gplspec-cleanup-baseline --target stage00_reference -j2
tests/run_stage00.sh /tmp/gplspec-cleanup-baseline/bin/stage00_reference
```
