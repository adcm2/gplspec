# Stage 01 handoff — repository hygiene

- Starting accepted base: `a73600de54652e4ac246f90f1de1c6f6eabf14a0` (`cleanup/base`).
- Candidate branch: `cleanup/01-hygiene`; fixed candidate commit is provided with this handoff.
- Changes: removed stale tracked CMake cache/build metadata; narrowed ignored generated files; removed five isolated prototype/stub files and two commented internal umbrella includes; removed temporary rotation-angle/check-marker output and solver matrix-dimension probes. File classifications and interface/reference retention decisions are in `stage-01-decisions.md`; every change is logged in `implementation_status.md`.
- Preserved: public/active headers, all independent spherical/aspherical potential and spectral-element numerical references, equations, unit conventions, solver iteration/error reports, and frozen baseline artifacts.
- Validation: `cmake --build /tmp/gplspec-stage00-repro2 --target stage00_reference clean_bench_1 phobos_gravity -j2`; `tests/run_stage00.sh /tmp/gplspec-stage00-repro2/bin/stage00_reference`; `git diff --check`. The harness passed three comparisons of 20,700 records with zero differences at `rtol=1e-13`, `atol=1e-15`; representative output matched byte-for-byte; both examples built and linked.
- Numerical changes: none observed. Stage 02 has not started.
- Independent Luna review: pending.
