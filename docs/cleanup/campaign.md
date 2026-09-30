# GPLSpec cleanup campaign

The campaign extracts existing GPLSpec/GSHTrans functionality while preserving numerical behavior and supported public interfaces. Numerical bug fixes, new GSHTrans capabilities, solver/discretization changes, optimizations, batching, real-field compression, and output refactoring are excluded.

## Stage status

| Stage | Branch | Status | Acceptance |
|---|---|---|---|
| 00 Baseline | `cleanup/00-baseline` | Accepted after independent Luna review of `dece12e59f86de59978e6cccfe684e22b4fcb06e` | Pinned reproducible build, immutable original-source reference, repeatable five-fixture regression and output check |
| 01 Hygiene | `cleanup/01-hygiene` | Not started | Documented deletion/interface decisions; examples and baseline regressions pass; Luna approval |
| 02 Headers | `cleanup/02-headers` | Not started | Public include/link checks and examples pass; baseline regressions pass; Luna approval |
| 03 Utilities | `cleanup/03-utilities` | Not started | Focused old/new helper comparisons and baseline regressions pass; Luna approval |
| 04 Models | `cleanup/04-models` | Not started | Intermediate arrays and solutions match frozen baseline; public constructors remain compatible; Luna approval |
| 05 Operators | `cleanup/05-operators` | Not started | Direct operator/source/perturbation comparisons and baseline regressions pass; Luna approval, then cumulative Sol review |

Stages begin from the latest accepted `cleanup/base`. Accepted stage branches are fast-forwarded only into `cleanup/base`. Keep `main` and `develop` unchanged, do not push, and stop after stage 05. Each stage requires its own fixed candidate commit and separate high-reasoning Luna review. Sol is reserved for the cumulative post-stage-05 checkpoint unless the specified concrete blocker occurs.

## Stage-00 comparison contract

The frozen reference is source commit `4ef3a66c62d64408d99989dd51c3ccbdc46d0b0c`, using the dependency revisions in `baseline-manifest.md`. It covers homogeneous and layered spheres, analytic lateral tomography values, a smooth aspherical mapping, a mapping perturbation, full model intermediate arrays, deterministic complex matrix-free operator products, force vectors, solutions, and `MatrixSolution.out`.

Repeated reference runs were byte-identical on the pinned compiler/dependency stack. Cross-stage checks use `rtol=1e-13` and `atol=1e-15`; no stage may relax these tolerances. The representative output file is compared byte-for-byte.
