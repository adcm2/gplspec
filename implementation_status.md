# Implementation status

## Stage 00 — baseline

- Build repair (isolated): the `clean_bench_*` targets included GPLSpec headers that transitively include `GSHTrans/All`, but their link interface did not provide GSHTrans include paths and dependencies. GPLSpec's interface target now links GSHTrans. This changes build propagation only; no numerical source was changed.
- Validation: an out-of-source configure succeeded against the locally available clean dependency stack; `phobos_gravity` compiled before the repair. `clean_bench_1` first failed with `GSHTrans/All: No such file or directory`, then configured, compiled, linked, and ran after the repair; the homogeneous-sphere solve reported 0 iterations and error `1e-06`.
