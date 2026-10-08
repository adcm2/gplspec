# Independent-residual audit artifacts

`gplspec-phase3-sol-residual-audit.md` is Sol's report copied byte-for-byte from the read-only baseline audit (SHA-256 `53ed55013e504c61a3823a892a4b177448aca9e48a6ac97c378af34698ef339e`). It contains exact settings, source interpretation, raw-result ranges, and limits on what the residual/backward-error values establish.

The two probe sources and final raw CSVs are also preserved here:

- File-backed PREM200/S40RTS probe: `gplspec-phase3-residual-probe.cpp`, SHA-256 `bf007eb203607a2455611e5e7dc51b7e506eb019d01dc6586a5cd647c4d17eef`; result `gplspec-phase3-residual-prem200-final.csv`, SHA-256 `b46c0d75dfc073f1f0aab56afe24a4012e63799bde129c86298554f82a6bb2f6`.
- Mapped clean_bench_3 probe: `gplspec-phase3-residual-mapped-probe.cpp`, SHA-256 `1d1fabe9900a45ddddb56a36122abb5ef5266ee73ac8871c4e36f5a821afbd5e`; result `gplspec-phase3-residual-mapped-final.csv`, SHA-256 `146b8176ba4b9be81181dc426c29392060bf2596a5841c106f74ec5d90a049e5`.

Both probes used the detached baseline commit, pinned build dependencies, and fixed one-thread settings. They are diagnostics and raw evidence, not frozen solution references or solver acceptance tests.
