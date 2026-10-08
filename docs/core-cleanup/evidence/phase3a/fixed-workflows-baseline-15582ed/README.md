# Fixed-workflow baseline capture

The raw manifest and small process logs were copied from `/tmp/gplspec-phase3-fixed-baseline/`. The manifest SHA-256 is `bef92a2b0a7f676fea6ae749d43810949ece5c81705fb7f93b1cf15fb09cb641`. It records baseline commit `15582edab8112dbb63e0101be5b26f5135298644`, one OpenMP thread, process return codes, wall times, and each generated output's relative path, byte count, and SHA-256. Raw `stdout.log`, `stderr.log`, and available `time.txt` files are preserved in the sibling directories.

All eight runs exited 0: `clean_bench_1` through `_4`, `clean_bench_6`, `phobos_gravity`, `phobos_heterogeneous`, and `OutputZiheng`. They produced 33 files. `clean_bench_1..4` each took below 0.12 s; these short runs are fixed-workflow checks, not controlled timings. `clean_bench_6` took 2.45 s and wrote five files (about 5.9 MB). The two Phobos runs wrote 7 and 9 files (about 129 MB and 182 MB) in 102.7 s and 82.6 s. `OutputZiheng` wrote three files (about 16.2 MB) in 74.0 s. See `manifest.json` for exact output identities.

Raw process-log provenance (all `stderr.log` files are empty, SHA-256 `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`):

| Program | `stdout.log` SHA-256 | `time.txt` SHA-256 |
|---|---|---|
| `clean_bench_1` | `47ca07b81e4e2781dd1ef9299ba4601145f310cc161eed05f62c9afdbcce65ff` | not captured |
| `clean_bench_2` | `884361c3143f75f9dc3e2e3396e529fb5820b0b5dd90fcb34a6f7fbf7f9cd4d1` | not captured |
| `clean_bench_3` | `95225d21c051fd65d9a3ced6b52e3477cb94b56634037a33df0f8d81ee0fc6d5` | not captured |
| `clean_bench_4` | `5d1e057c474cc0aa1216b91b88c5b50aa1aab6b57d68e4204db69851c0bcaeb9` | not captured |
| `clean_bench_6` | `2d6d83dbba96412fa5a6f4af81a9f35b374a3aaa6b27f6ff7c8e2d29dd3554e1` | `56acca4a12bcf07338a2c37d4588bbc6e99efbc45e28390271d7afc5376868b8` |
| `phobos_gravity` | `ad6bf4b228fa7ed92bd1fe9edd6c063a06d57d6e2340fb42e1bfac085ff290cc` | `61e85ff65355ba6e01fabf9657c32dc4508db3628f3f8525ec5f861731864bd3` |
| `phobos_heterogeneous` | `c611ffd78afe09fa68287937a15ae0b73558f4d89ab581fb06e4af2bd9b1cf6d` | `8b37eee022bb6ecde5486aa166ddd6cd57f20052912d350e13a0c78d7bbed242` |
| `OutputZiheng` | `fcc7dc959a62f4584f5473767ee19a37b7520e0e0de6577ab82f201ffe37dc80` | `1c46e35559ef192fbba89e3d88a4896e16601b8ad97291de2ca3993a6b17af8f` |

The generated outputs remain isolated under `/tmp/gplspec-phase3-fixed-baseline/<program>/work/`; they were not copied over tracked workflow outputs. No paper mapping, numerical accuracy, or controlled performance claim follows from exit status or text scanning. The Phobos gravity variant produced seven files; no absent physical-rotated file is inferred. The parent report's textual NaN/Inf scan is a smoke scan, not a numerical comparison.

The detached baseline checkout provides source, build, and input provenance; it is not the runtime working directory. The unchanged executables were built from that checkout in `/tmp/gplspec-phase3-baseline-release`. To recreate a lost output tree, run `/tmp/gplspec-phase3-baseline-release/bin/<program>` from a fresh isolated directory with `modeldata` linked to `/tmp/gplspec-phase3-baseline/modeldata`, using `OMP_NUM_THREADS=1`, `OMP_DYNAMIC=FALSE`, `OPENBLAS_NUM_THREADS=1`, and `MKL_NUM_THREADS=1`. Create the workflow's required output directories there first: `work/Bench1`, `work/Bench2`, `work/Bench3`, `work/Bench4`, or `work/Bench6` (also `work/Bench6/lm`); Phobos uses `work/Phobos/Cartesian`, `work/Phobos/Spherical`, and `work/Phobos/Sensitivity`; `OutputZiheng` uses `work/Ziheng`. Compare each regenerated file with the per-file identities in `manifest.json`. The separate interactive `phobos_gravity_perturbation` run remains unmapped because the requested `lsmall` input has not been supplied.
