# Phase 3E.3 focused evidence

Candidate source is the working tree based on commit `47b503b`; the independent clean_bench baseline is commit `15582edab8112dbb63e0101be5b26f5135298644`. Frozen references were not modified. The protected `experimental/` tree was not changed or built.

The focused test is `stage3e_classical_perturbation`. Build and run it with:

```sh
cmake -S . -B /tmp/gplspec-phase3-resume/candidate-build
cmake --build /tmp/gplspec-phase3-resume/candidate-build --target stage3e_classical_perturbation -j1
ctest --test-dir /tmp/gplspec-phase3-resume/candidate-build -R '^stage3e_classical_perturbation$' --output-on-failure
```

The test checks an analytic normalized `Y10` source coefficient over a two-density layered sphere with separate one-sided interface traces, manufactured `r^2 Y20` advection with layer-specific displacement and exterior taper, exact zero-map source/advection, and full-mapped versus classical field tangents for the clean_bench_10 spherical problem. The mapped fixture uses `p=5`, `lmax=2`, outer radius `1.5`, maximum element steps `0.1` and `0.05`, solver tolerance `1e-12`, and central-difference amplitudes `1e-3` and `5e-4`. It compares at the same reference nodes, uses positive spherical volume weights and full-harmonic Parseval RMS in `m^2/s^2`, checks referential source invariance, and applies the predeclared individual tangent and noise bounds in `phase-3.md`.

The focused Release CTest passed 1/1. At mesh step `0.1`, the classical and full-mapped tangent RMS values are `7.14099e7` and `7.14100e7 m^2/s^2`; the `1e-3` tangent difference is `24.5314` against a propagated solve-noise floor of `24.5305`. At `5e-4`, the difference is `7.2756` against floor `49.0609`, so the amplitude ratio `3.37173` is descriptive and noise-limited. At fixed amplitude `1e-3`, the mesh-step `0.05` difference is `24.0585`, below its floor `24.5305`; the mesh-refinement comparison is therefore noise-limited. Both amplitudes pass the fixed individual error bound; no solver algorithm or workflow setting changed.

Clean_bench 9 and 10 retain their original settings (`p=5`, `maxstep=0.1`, `outer radius=1.5`, `h=0.001`, and `lmax=20`/`30`). Reproduce in fresh output directories so the programs' relative output paths stay isolated:

```sh
mkdir -p /tmp/bench9/work/Bench9/{exact,perturbed,unperturbed}
(cd /tmp/bench9 && printf '0.001\n' | /tmp/gplspec-phase3-resume/candidate-build/bin/clean_bench_9)
mkdir -p /tmp/bench10/work/Bench10/{exact,perturbed,unperturbed,perturbed2}
(cd /tmp/bench10 && printf '0.001\n' | /tmp/gplspec-phase3-resume/candidate-build/bin/clean_bench_10)
```

Both candidate workflows exited `0` with empty stderr. Candidate outputs and hashes are currently captured under `/tmp/gplspec-phase3-resume/candidate-classical-workflows/bench9` and `.../bench10`; exact baseline outputs and their manifest are under `/tmp/gplspec-phase3-resume/baseline-classical-workflows`. The candidate output digests are:

| Workflow | Output | SHA-256 |
|---|---|---|
| clean_bench_9 | exact/MatrixSolution.out | `3815485232acd8a79ca094494e7cd75f2908d572a0239e0d631d05618b56f5b0` |
| clean_bench_9 | perturbed/MatrixSolution.out | `7383dabc57ed71b113a21c0a9d9ca90692810388e4288b226e06afa652ff6553` |
| clean_bench_9 | unperturbed/MatrixSolution.out | `2ee7d98a34ff2ed41bbfa9cb80aa40baa598f240d156cb124af1482f9ff66a6d` |
| clean_bench_10 | exact/MatrixSolution.out | `c93399b1b381abf0f128730aee8ed8a37bedc6ebe8d816ee4e1cf25a716eae70` |
| clean_bench_10 | perturbed/MatrixSolution.out | `452a319608addbe85559138381d81083029d6edbde53afd31efa7a24323a377b` |
| clean_bench_10 | unperturbed/MatrixSolution.out | `ab24b6df090001bcd066d5259ec07b28def3f2bc9d9d49059cde09e491fd14e1` |
| clean_bench_10 | perturbed2/MatrixSolution.out | `19b78d1fc8732f73a31443856a98b29d6d6d2507bad65c1c271a8fd6394ffc77` |

These semicolon-delimited outputs contain repeated `(radius, theta, phi, real potential, imaginary potential)` tuples, with radii in meters and potential in `m^2/s^2`. The matched baseline/candidate comparison found six of seven files within the fixed pointwise physical-potential floor `1.7471e-4 m^2/s^2`; `clean_bench_10/perturbed/MatrixSolution.out` differs by up to `4.6181e5 m^2/s^2`. That classical-output discrepancy fails the declared criterion and is under bounded source/advection attribution; no workflow-equivalence pass is claimed here. Candidate digests establish captured-file identity only.

## Bounded attribution of the clean_bench_10 classical-output difference

The predeclared candidate comparison still records a real failure: baseline commit `15582ed` and the current workflow differ by up to `461813.485791 m^2/s^2` in `clean_bench_10/perturbed/MatrixSolution.out`, above the fixed pointwise floor. This is not a 3E.3 arithmetic regression: the exact 3E.2 parent (`47b503b`) already differs from baseline by the same amount, while the 3E.3 candidate differs from that parent by at most `3.33786e-6 m^2/s^2`. The original baseline-vs-candidate failure remains in the record and is not relabeled as a tolerance pass.

A bounded solver-component probe localized the difference to the baseline boundary-source vector. At `p=5`, `lmax=30`, `maxstep=0.1`, outer radius `1.5`, and `h=0.001`, baseline source expected to be degree one also contains spurious degree norms: `l=0: 1.10228e-4`, `l=2: 2.36047e-4`, and through `l=30: 4.43354e-5`; the intended `l=1` norm is `0.519373`. The accepted parent source is pure `l=1` with norm `0.519184` (other degrees below `1e-12`). Base RHS norms are identical (`60.2221`). Baseline and parent solves both report Eigen `Success`; the base solves take one iteration, boundary corrections zero iterations from their preconditioned guesses, and independently recomputed residuals are respectively `3.07703e-11` / `2.31576e-14` and `3.3836e-11` / `2.31993e-14`. Solver convergence therefore does not explain the artifact.

To isolate the suspected geometry lifetime, the baseline source was copied to a scratch include tree and only `Density3D::GSH_Grid()` was changed from returning the grid by value to returning `const Grid&`. Rebuilding and rerunning the same baseline probe removes the leakage: only `l=1` is printed, with source norm `0.519184`; the base and boundary residuals are `3.07703e-11` and `2.31582e-14`. This matches the accepted parent to roundoff. In the tested GNU C++ 13.3.0 build, the old by-value accessor made `for (auto theta : model.GSH_Grid().CoLatitudes())` consume borrowed views backed by a temporary grid. The accessor ownership correction was already made and reviewed in 3D (`5ee1ed0`); the scratch test is attribution evidence, not a new production patch or a claim about every standard-library implementation.

Probe source, compact logs, two-line scratch patch, and reproduction command are under `attribution/`. Large component CSVs and scratch build artifacts remain in `/tmp/gplspec-phase3-resume/` and are not needed to reproduce the summary. This attribution explains the historical baseline discrepancy; it does not establish full workflow equivalence or replace the predeclared physical gate.

The earlier Stage 05 status statement that `MappingTools::dxitodf` remained unchanged is historical. The current 3E.3 caller inventory found no active tracked consumer, so the helper has been retired; frozen references and historical reports remain untouched.
