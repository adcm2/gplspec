# Stage 3C.3 validation evidence

The paired baseline is the clean Git archive at `15582edab8112dbb63e0101be5b26f5135298644`; the candidate is the tree based on `472c146` with the Stage 3C.3 edits. Both use the same `phase3a_performance.cpp` source (SHA-256 `1781a778ff5f665c070884934bd17c4199691bba06f6fe92e88dc252d2ba62e2`), the same model inputs, Release flags, pinned dependencies, one thread, and 0.25-second measurement batches. Baseline and candidate manifests verify the input/reference hashes and exact process outcomes. The copied raw candidate measurements and baseline three-run summary/manifests are alongside this file.

Build and run the focused checks:

```sh
cmake --build /tmp/gplspec-phase3-resume/candidate-build --target \
  stage1c_file_tomography_solver stage3c_source_discretization \
  stage3c_operator_action phase3a_performance -j1
ctest --test-dir /tmp/gplspec-phase3-resume/candidate-build --output-on-failure \
  -R '^(stage1c_file_tomography_solver|stage3c_source_discretization|stage3c_operator_action)$'
```

The focused CTests passed 3/3. The solver fixture reported zero-RHS success with zero residual, retained the PREM 1-D BiCGSTAB and mapped 3-D CG guess/update checks, and compared the public result with its independent residual calculation. The performance executable can be run as `phase3a_performance mapped 0.25` and `phase3a_performance file_backed 0.25` from the project root, repeated three times per scenario with `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and `OMP_DYNAMIC=FALSE`.

The timing harness manually constructs the same Eigen solver lifecycle on both source snapshots. It does not time the new `PreparedSolver` public result API or its extra independent-residual matrix-vector product; that cost remains for the bounded 3F performance check. Three-run medians (baseline → candidate) were:

| Scenario | End-to-end solve batch | Peak RSS median | Independent residual |
|---|---:|---:|---:|
| Mapped | 17.413 → 15.786 ms | 17,232 → 17,224 KiB | `2.07984e-12` → `2.69037e-12` |
| File-backed | 436.140 → 420.600 ms | 573,656 → 564,340 KiB | `5.48440e-10` → `6.39414e-10` |

All baseline and candidate performance processes returned historical exit code 3 only because the independent residual exceeded the requested tolerance. This is recorded as residual evidence, not a blanket accuracy gate. The candidate retains successful Eigen status and finite outputs.

Because the benchmark checksum is only a solution norm, a paired full-vector capture was also made with `solution-probe.cpp`, a copy of the same performance harness with the first end-to-end complex solution, `PotentialNorm`, physical `Mass()`, and physical outer radius written to CSV. It was compiled separately against the exact baseline archive and candidate tree, then run for both fixture scenarios. Probe source SHA-256: `f4647282f9bdab667020c6c1d3e69a8d34d7f45efa03b710b5d03e511290d510`. The four vector CSVs remain in `/tmp/gplspec-phase3-resume/` with hashes:

To recreate the probe binaries, compile this source once with `source_root` set to the baseline archive and once to the candidate tree, using that snapshot's dependency build paths:

```sh
c++ -O3 -DNDEBUG -std=gnu++23 -fopenmp \
  -I"$source_root" -I"$deps/gaussquad-src" \
  -I"$deps/numericconcepts-src/include" -I"$deps/gshtrans-src" \
  -I"$deps/fftwpp-src" -I"$deps/planetarymodel-src" \
  -I"$deps/interpolation-src" -I"$deps/tomographymodels-src" \
  -isystem "$deps/eigen3-src" solution-probe.cpp -o solution-probe \
  -fopenmp /usr/lib/x86_64-linux-gnu/libfftw3.so \
  /usr/lib/x86_64-linux-gnu/libfftw3f.so \
  /usr/lib/x86_64-linux-gnu/libfftw3l.so -lnetcdf
./solution-probe mapped 0.02 /tmp/solution-mapped.csv
./solution-probe file_backed 0.02 /tmp/solution-file-backed.csv
```

Run from the repository root so the fixture input paths resolve. Each probe exits 3 at the existing requested-residual gate after writing the solution; the vector output and independent physical comparison are still valid. The exact original source for the shared harness is `tests/phase3a_performance.cpp`, SHA-256 `1781a778ff5f665c070884934bd17c4199691bba06f6fe92e88dc252d2ba62e2`.

| File | SHA-256 |
|---|---|
| `solution-baseline-mapped.csv` | `555c03de31a2907c6cacf996b1dea6c846243840c9aff3063447e2cdd758cbfc` |
| `solution-candidate-mapped.csv` | `072856c9891077b8d1a3d81df783e272fa2bef340c0a7921c86ceb23cbd1483e` |
| `solution-baseline-file_backed.csv` | `3f1508aeb4ffc0ede25cad787d406f7657cb8a48695e0905d5768e4aae235c01` |
| `solution-candidate-file_backed.csv` | `e52e67f47e79c30cbc4412be5e8c1be7dd3409c478b18251dd36f6a5f87ba8b9` |

At each radial node, group the full complex harmonic coefficient differences into `N=(lmax+1)^2` entries. The addition theorem bounds the physical pointwise potential error by `(lmax+1)/sqrt(4*pi) * PotentialNorm * max_node(||delta_a||_2)`. A pointwise bound also bounds every positive-quadrature weighted RMS. Under the predeclared potential criterion `rtol=1e-8`, `atol=1e-12*S_q`, `S_q=max(G*M/R,RMS(|q_ref|))`, the absolute floor is at least `1e-12*G*M/R`. Using the physical `Mass()` and outer radius in each exact fixture gives:

| Scenario | Mass (kg) | Outer radius (m) | `GM/R` (m²/s²) | Conservative pointwise bound (m²/s²) | Lower-bound absolute floor (m²/s²) |
|---|---:|---:|---:|---:|---:|
| Mapped, `lmax=2` | `2.501545510262682e25` | `8,919,400` | `1.871882520288863e8` | `1.356e-6` | `1.872e-4` |
| File-backed, `lmax=20` | `5.973168921908029e24` | `7,645,200` | `5.214609205154102e7` | `7.119e-6` | `5.215e-5` |

Each physical upper bound is below even the conservative lower bound on the fixed absolute tolerance. `RMS(|q_ref|)` can only increase the actual `S_q`; no assumption that model mass normalization equals physical fixture mass is used. The independent residual changed slightly from baseline and is reported above; the full physical field comparison passes the already-declared criterion without changing the solver or tolerance.

Candidate measurement CSV SHA-256: `ad3b43cfc924418b9542fed1aec370e0428ec6dad43401386279a5926aea0f41`; candidate manifest SHA-256: `f92d59efb0f7aec3873a0d4b53fbb6a7c3945f04396a63c2b366aba6791391ae`. Baseline summary SHA-256: `1a9c91feff4248524a14af30d6aa2630e758034c98ee17f4239ec61258405d13`; baseline manifest SHA-256: `f12f2b4cb77e30cfebf6237d651a47ae22b870e4a4dc0d71cc3ff027ee5dfb21`. The candidate performance binary SHA-256 is `960cb7e898dd494a628ab25fc3c792675033ca4ccaf8b4ead13a82a7e38ab11f`.
