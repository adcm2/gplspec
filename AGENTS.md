# GPLSpec cleanup campaign instructions

## Scope and review gates

- Follow the latest explicit user authorization. Older status notes and stop boundaries are historical records and do not override a newer user instruction.
- Keep changes within the currently approved phase and exact file scope. After each logical code, build, or documentation change, record what changed and its validation state in `implementation_status.md`.
- Use one active writer. Freeze the candidate before independent review; reviewers inspect a fixed snapshot without editing it.
- Preserve public interfaces, mathematical operations, solver behavior, coefficient conventions, and output formats unless a later approved proposal says otherwise.
- Do not regenerate frozen numerical references, loosen tolerances, or claim broader scientific validation than the fixtures support.
- Preserve every user-added file under `experimental/`. Do not move, rewrite, automatically build, or claim numerical validation of those experiments unless explicitly authorized.
- Keep stage changes uncommitted and do not integrate, tag, push, or delete branches unless the user explicitly authorizes that close-out step.

## Cleanup provenance

`docs/cleanup/campaign.md`, `docs/cleanup/deferred-issues.md`, and `docs/cleanup/baseline-manifest.md` contain the campaign record, preserved observations, and dependency pins. Stage 00's source commit and its checked-in references are immutable. The tested system versions of Eigen, FFTW, and NetCDF are documented; only dependency revisions explicitly pinned by the build are enforced. Stage-07 evidence and review gates are in [docs/cleanup/stage-07-handoff.md](docs/cleanup/stage-07-handoff.md).

## Validation workflow

Use a fresh out-of-source build for a new final-candidate validation. The standard command enables all examples and the baseline harness:

```sh
cmake -S . -B /tmp/gplspec-build \
  -DMY_PROJECT_BUILD_EXAMPLES=ON \
  -DGPLSPEC_BUILD_BASELINE_HARNESS=ON
cmake --build /tmp/gplspec-build -j2
ctest --test-dir /tmp/gplspec-build --output-on-failure
tests/run_stage00.sh /tmp/gplspec-build/bin/stage00_reference
/tmp/gplspec-build/bin/stage02_header_link
```

Run representative applications from isolated working directories because they write scientific output files relative to the current directory. For changes that touch code, tests, build files, or validation provenance, record compiler and relevant dependency versions and compare the code/reference checksums before and after validation.
