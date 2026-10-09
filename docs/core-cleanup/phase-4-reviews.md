# Phase 4 review ledger

## Authorization and baseline

The user approved the documented Phase 3 corrections and explicitly authorized Phase 4 stages 4A–4F. Exact clean starting point: `refactor/core-cleanup` at `a8b1ca99fdadd9541a4227129c8fce872c98c1aa`. The isolated baseline archive and protected manifest are `/tmp/gplspec-phase4/baseline` and `/tmp/gplspec-phase4/baseline-identity.json`; the manifest covers 97 protected dataset/output/reference/experimental paths. Plotting source `work/*.py` is authorized to migrate, but saved products remain protected. See `phase-4.md` for scope, inventory and planned criteria.

## 4A — inventory and minimal design

Independent Luna review: **PASS**, no blockers. The current proposal separates callback-based bounded field evaluation, `manifest.json` + long-form `samples.csv` serialization, and Python plot readers. It distinguishes numerical element boundaries from material interfaces, labels both modeled sides by increasing reference radius, associates continuous potential traces with a shared degree of freedom, and defines the no-exterior-element contract. The required quantitative schema, workflow migration, predefined comparison criteria, Phobos batched-sampling performance constraint and protected experimental include are inventoried in `phase-4.md`. The design was reviewed as a plan only; no production implementation or protected-file modification is part of 4A. Stage 4A is accepted pending coordinator local commit. Stage 4B is next after that commit.
