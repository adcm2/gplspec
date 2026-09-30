# GPLSpec cleanup campaign instructions

The currently authorized campaign boundary is behaviour-preserving cleanup through stage 02. The user has explicitly instructed the work to stop after stage 02; stages 03–07 and GSHTrans modernization are not authorized in this run. Preserve the mathematical operations, public interfaces, output formats, and provenance described in `docs/cleanup/campaign.md`.

Before each stage, start from the latest accepted `cleanup/base` commit and use its named `cleanup/NN-*` branch. Keep `main` and `develop` unchanged. Do not push. A stage advances only after its fixed candidate commit passes its independent Luna review. Do not start a stage before the preceding stage receives independent Luna approval and is integrated into cleanup/base. Stop after stage 02 is accepted and integrated. Do not begin stages 03–07, request the stage-05 Sol checkpoint, modernize GSHTrans, or start performance work.

Keep one active implementation writer. Review the candidate in an isolated read-only worktree. Record each code/build change and validation result in `implementation_status.md`. Never regenerate the frozen stage-00 baseline from a later stage.

Stage-00 baseline command:

```sh
cmake -S . -B /tmp/gplspec-cleanup-baseline \
  -DMY_PROJECT_BUILD_EXAMPLES=OFF \
  -DGPLSPEC_BUILD_BASELINE_HARNESS=ON
cmake --build /tmp/gplspec-cleanup-baseline --target stage00_reference -j2
tests/run_stage00.sh /tmp/gplspec-cleanup-baseline/bin/stage00_reference
```
