# Phase 3C.1 source/discretization evidence

Candidate base: accepted checkpoint `7e0b2780419edae52b606f800d00dd9e2b4df41d`; all frozen reference files retain their historical content. SHA-256: Stage 00 `ef20bd7e3016661c60903c290c74c599151fdfc0417443bcf81c4215401f105e`; Stage 05 1-D `a9e89b9ac1250f3f01b695ec5be0fcb79990a85da42e6eb80ae001f5f88a43d6`; Stage 05 file perturbation `08cf04f1243fa128645490c4eba13eed5f7b5fceff76b032f2f5d8141efe3263`. Candidate source SHA-256: `Gravity_Tools.h` `18db0ecbc18ddea6def1a3764b970f3f1117613edfa63e7c8753eb01411c9e21`; `Radial_Tools.h` `c4e6956d7b0436971d49820e22d5c19455712da1a9deb79c8a30df6a22b9cf92`; `Earth_General_Models_1D.h` `4dc0e28e60abbc6986c256581671f94b9a1aad4e04f30b6a099e4f56d722d21a`; `tests/CMakeLists.txt` `286f506f003ee86a13d9cc6e66d371d1ad1c87209895ae5712361b89b684292c`; test `stage3c_source_discretization.cpp` `ac2d0486dd4b701d45b868e3deca1b481487603ee507ee5b83563deed5a5a611`. The source/discretization contract and physical comparison criteria are in `docs/core-cleanup/phase-3.md`.

The focused Release build produced `stage3c_source_discretization` and the existing `stage1c_file_tomography_solver`. Both passed CTest (2/2). The tests cover canonical global radial indexing and shared element endpoints, one-sided interface density samples, direct monopole quadrature, low-degree 3-D source coefficients, and source zero slots.

The unchanged Stage 00 harness retained all 20,700 record keys. The componentwise exact harness reports roundoff differences, with 17, 58, 17, 189, and 198 changed components in homogeneous, lateral-density, layered, mapping-perturbation, and smooth-mapping cases respectively. The potential solutions contain `lmax=2` complex orthonormal harmonic coefficients. At each radial node, the addition theorem gives the physical pointwise bound

`|delta Phi| <= (lmax+1)/sqrt(4*pi) * ||delta a||_2 * PotentialNorm`.

Using the maximum per-radius coefficient-vector norm for each case and `PotentialNorm=(6371000/3600)^2`, the respective bounds are `1.51e-7`, `4.52e-7`, `2.11e-6`, `6.03e-7`, and `6.03e-7 m^2/s^2`. Each is below the predeclared absolute field floor `1e-12*(GM/R)=6.26e-5 m^2/s^2`; the same pointwise bound also bounds positive quadrature-weighted RMS. The frozen reference and strict comparator remain unchanged.

Stage 05 1-D has 1,178 numeric records. Twenty-nine source scalar components changed across 29 source rows, maximum absolute difference `3.55e-15` and maximum scaled error `6.96e-16` using the predeclared `1e-13` source/algebra floor. Stage 05 file-perturbation has 3,818 numeric records plus its recognized `Finding a` preamble. The summarizer reports seven source scalar components changed across seven source rows (maximum absolute difference `1.42e-14`, scaled error `1.92e-16`) and 104 scalar components changed across 87 complex solution coefficients. The maximum per-radius coefficient-vector difference is `1.71e-13`, giving the addition-theorem pointwise potential bound `4.52e-7 m^2/s^2`, below the predeclared absolute floor `6.26e-5 m^2/s^2`. The unchanged strict comparators report arithmetic differences as diagnostics.

Candidate captures used for the summaries are under `/tmp/gplspec-stage3c1-stage00-2lmkkucb/actual.csv`, `/tmp/gplspec-stage3c1-stage05-1d.csv`, and `/tmp/gplspec-stage3c1-stage05-file.csv`. With the configured Release tree `/tmp/gplspec-phase3-resume/candidate-build`, reproduce the build and focused test using:

```sh
cmake --build /tmp/gplspec-phase3-resume/candidate-build --target stage3c_source_discretization stage1c_file_tomography_solver stage00_reference stage05_operator_reference stage05_file_perturbation_reference -j1
ctest --test-dir /tmp/gplspec-phase3-resume/candidate-build --output-on-failure -R 'stage3c_source_discretization|stage1c_file_tomography_solver'
```

Capture comparisons in isolated temporary outputs:

```sh
/tmp/gplspec-phase3-resume/candidate-build/bin/stage00_reference /tmp/phase3c1/out /tmp/phase3c1/stage00.csv
/tmp/gplspec-phase3-resume/candidate-build/bin/stage05_operator_reference > /tmp/phase3c1/stage05-1d.csv
/tmp/gplspec-phase3-resume/candidate-build/bin/stage05_file_perturbation_reference > /tmp/phase3c1/stage05-file.csv
python3 docs/core-cleanup/evidence/phase3c1/summarize.py tests/reference/stage00-baseline.csv /tmp/phase3c1/stage00.csv
python3 docs/core-cleanup/evidence/phase3c1/summarize.py tests/reference/stage05-1d-original.csv /tmp/phase3c1/stage05-1d.csv
python3 docs/core-cleanup/evidence/phase3c1/summarize.py tests/reference/stage05-file-perturbation-original.csv /tmp/phase3c1/stage05-file.csv
```

Frozen Stage 00/Stage 05 references were not edited. No full workflow or performance campaign was repeated for this source-only checkpoint.
