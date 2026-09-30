# Stage 02 handoff — header organization

- Starting accepted base: `62e5e315301de4b198edf7f6d4724fd47f364d56` (`cleanup/base`).
- Branch: `cleanup/02-headers`. Dependency-only commit: `49006f882550c79c8cd6a08739a2b547507a6124`; fixed candidate reviewed: `944181c77a5a02003144b4945b492ba60993c56e`.
- Public umbrella paths are unchanged: `gplspec/Core`, `gplspec/SimpleModels`, `gplspec/GeneralModels`, and `gplspec/All`. Their roles and implementation-fragment boundaries are in `stage-02-decisions.md`. Matrix wrappers now include model declarations they directly name, removing include-order assumptions.
- GPLSpec header ODR repairs add `inline` only to affected non-template class declarations and header free functions; function bodies, member layout, ownership, and mathematical operations are unchanged.
- The first full `<gplspec/All>` two-TU link exposed seven external duplicate definitions in the exact pinned dependency revisions: FFTWpp `CleanUp` (`FFTWpp/src/Core.h`); `ExportWisdom`, `ImportWisdom`, `ForgetWisdom` (`FFTWpp/src/Wisdom.h`); TomographyModels `Tomography::GetValueAt`, `ReverseLatitude` (`TomographyModels/src/Tomography.hpp`), and `ShellExec` (`TomographyModels/src/ShellExec.hpp`). Local patch files add only `inline`; CMake applies them after dependency population in per-build FetchContent copies and checks all seven markers. No dependency pin or shared cache was changed.
- Validation: `cmake -S . -B /tmp/gplspec-stage02-build -DMY_PROJECT_BUILD_EXAMPLES=ON -DGPLSPEC_BUILD_BASELINE_HARNESS=ON`; `cmake --build /tmp/gplspec-stage02-build -j2`; `tests/run_stage00.sh /tmp/gplspec-stage02-build/bin/stage00_reference`; `/tmp/gplspec-stage02-build/bin/stage02_header_link`; `git diff --check`. Four umbrella smoke translation units compiled; the `All` two-TU consumer linked and ran; every example target built (clean_bench_1–10, phobos_gravity, phobos_gravity_perturbation, phobos_heterogeneous, wignertest, OutputZiheng). All three numerical comparisons passed with zero differences across 20,700 records; representative output was byte-identical.
- Numerical differences: none observed. Deferred issue: dependencies still contain their upstream definitions; the reproducible stage-local CMake patch repairs only the private build copies.
- Implementer model: `gpt-6-luna` medium. Independent Luna review: pending.
- Stop boundary: after stage 02 review and integration, stop. Stages 03–07, the stage-05 Sol checkpoint, and GSHTrans modernization are not authorized by the user's latest instruction.

## User-requested review stop

Stage 02 is committed but NOT accepted or integrated. The user requested cutting the independent review short; it must resume before acceptance. See `docs/cleanup/review-stage-02.md`. `cleanup/base` remains at accepted stage 01 (`62e5e315301de4b198edf7f6d4724fd47f364d56`). No later stage was started.
