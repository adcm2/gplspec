# GPLSpec cleanup campaign

This campaign extracts existing GPLSpec/GSHTrans functionality while preserving numerical behavior and supported public interfaces. Numerical bug fixes, new GSHTrans capabilities, solver/discretization changes, optimizations, batching, real-field compression, and output refactoring are excluded. The user has updated the stopping boundary: complete and review stage 02, integrate it, then stop; stages 03–07 and the stage-05 Sol checkpoint are not authorized in this run.

## Stage status

| Stage | Branch | Status | Acceptance |
|---|---|---|---|
| 00 Baseline | `cleanup/00-baseline` | Accepted after independent Luna review of `dece12e59f86de59978e6cccfe684e22b4fcb06e` | Pinned reproducible build, immutable original-source reference, repeatable five-fixture regression and output check |
| 01 Hygiene | `cleanup/01-hygiene` | Accepted after independent Luna PASS at `3f7ccf46625d5b09d13d966c58573d437deb5c99` | Documented deletion/interface decisions; examples and baseline regressions pass; Luna approval |
| 02 Headers | `cleanup/02-headers` | Candidate `944181c77a5a02003144b4945b492ba60993c56e` complete; awaiting independent Luna review | Public include/link checks and examples pass; baseline regressions pass; Luna approval |
| 03 Utilities | `cleanup/03-utilities` | Not authorized under current stop boundary | Focused old/new helper comparisons and baseline regressions pass; Luna approval |
| 04 Models | `cleanup/04-models` | Not authorized under current stop boundary | Intermediate arrays and solutions match frozen baseline; public constructors remain compatible; Luna approval |
| 05 Operators | `cleanup/05-operators` | Not authorized under current stop boundary | Direct operator/source/perturbation comparisons and baseline regressions pass; Luna approval, then cumulative Sol review |

Stages begin from the latest accepted `cleanup/base`. Accepted stage branches are fast-forwarded only into `cleanup/base`. Keep `main` and `develop` unchanged, do not push, and stop after stage 05. Each authorized stage requires its own fixed candidate commit and separate high-reasoning Luna review. No work beyond stage 02 is permitted under the current user instruction.

## Stage-00 comparison contract

The frozen reference is source commit `4ef3a66c62d64408d99989dd51c3ccbdc46d0b0c`, using the dependency revisions in `baseline-manifest.md`. It covers homogeneous and layered spheres, analytic lateral tomography values, a smooth aspherical mapping, a mapping perturbation, full model intermediate arrays, deterministic complex matrix-free operator products, force vectors, solutions, and `MatrixSolution.out`.

Repeated reference runs were byte-identical on the pinned compiler/dependency stack. Cross-stage checks use `rtol=1e-13` and `atol=1e-15`; no stage may relax these tolerances. The representative output file is compared byte-for-byte.

## User-requested review stop

Stage 02 is committed but NOT accepted or integrated. The user requested cutting the independent review short; it must resume before acceptance. See `docs/cleanup/review-stage-02.md`. `cleanup/base` remains at accepted stage 01 (`62e5e315301de4b198edf7f6d4724fd47f364d56`). No later stage was started.
