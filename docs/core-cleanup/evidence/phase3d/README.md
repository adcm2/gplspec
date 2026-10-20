# Phase 3D harmonic workspace and action timing

The Stage 3D independent review passed. Focused Release CTests passed 5/5: `stage1c_file_tomography_solver`, `stage3a_harmonic_contract`, `stage3b_model_input`, `stage3c_source_discretization`, and `stage3c_operator_action`. The action test interleaves two inputs and repeats the first, checking operation-local buffers do not retain transformed values. `git diff --check` is clean. No frozen reference or protected path changed.

The action timings use the existing `phase3a_performance` mapped and file-backed fixtures, one thread, and 0.25-second batches. The before values are the three raw operator-action measurements in `../phase3c3/candidate-measurements.csv`; that capture is immediately before workspace reuse (3C.3), while 3C.4 changed geometry ownership but did not change action formula/workspace. This is a bounded before/after comparison across those accepted snapshots, not a claim that 3D alone caused the full difference. The same harness source is `tests/phase3a_performance.cpp` (SHA-256 `a15305d3f38ddff20a91ff1827ad57e938762afe5f702097dba70fa5ae1e0797`). The post-change Release binary SHA-256 is `2d729252a6044baa57b164fbec0ace595e74b8f9ee43f2c5c9a9b2e39ec0d2c1`.

| Fixture | Before median | After run medians | After median | Change |
|---|---:|---|---:|---:|
| Mapped action | 0.893441 ms | 0.905499, 0.874380, 0.815001 ms | 0.874380 ms | -2.1% |
| File-backed action | 34.5646 ms | 34.0890, 33.7631, 33.5433 ms | 33.7631 ms | -2.3% |

The matched setup uses `OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_DYNAMIC=FALSE`, running `phase3a_performance mapped 0.25` and `phase3a_performance file_backed 0.25` three times each from the project root. The raw post-change CSVs are stored beside this README. All six commands returned exit 3 at the harness's historical independent requested-residual gate. The mapped residual was `2.69037e-12`; the file-backed residual was `6.39414e-10`. These observations are reported as workflow evidence and did not cause a solver change. The harness `operator_action` metric isolates the repeated matrix-vector action; its solver and end-to-end metrics are not used for this comparison.

Raw CSV SHA-256:

| File | SHA-256 |
|---|---|
| `mapped-1.csv` | `16a4309a163a46dfa10d43d6b3a63fa88f330dfa39c7961f73364b4fba29e795` |
| `mapped-2.csv` | `3c41abe131b0cad61463772d8dd4776d8ccc72986f8bd6f1ab5eb67acd10f86c` |
| `mapped-3.csv` | `d13d40dbef08657bace6074c8bf17268a1148485eaf5870b3a5592cffb9c1170` |
| `file-backed-1.csv` | `f4ec58347d82926e7c6fb1c9f1ee8bef3eec51fc824871db0a2d438fd81d0865` |
| `file-backed-2.csv` | `b29b28f0712cb20483d260e3afc8f3442c34b8a16e6ce9a63372acdd42feaa26` |
| `file-backed-3.csv` | `5b4be9464ae9de07b3163cbaa7ea3ae9258d0c54b5a6c5e89abb17e9e083104a` |

Review correction: an initial timing note reported a mapped median of 0.8447 ms by taking the median of a subset of the run medians. The summary above correctly uses all three raw runs, yielding 0.874380 ms. No timing was rerun; the CSV captures remain unchanged.
