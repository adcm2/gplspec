# PR 07 description

This PR combines pre-baseline model/output behavior changes with a later source cleanup. The cleanup comparison starts at `4ef3a66`; it does not represent equivalence to `main` (`4fd6488`).

## Behavior changes before the cleanup baseline

- Compared with `main` (`4fd6488`), rotated output at duplicated internal radii uses the inner-side trace: the preceding element's upper node. The first and final endpoints keep their prior element-side selection. For density plots this includes interior material through the surface, then exposes the vacuum endpoint. This behavior predates the cleanup baseline and is intentionally retained for plotting interior density up to the surface.
- `clean_bench_2` changes mesh `maxstep` from `0.01` to `0.1`; `phobos_gravity_perturbation.cpp` changes requested solve tolerance from `1e-12` to `1e-6`. These change example resolution and solver stopping behavior relative to `main`.
- `MappingPerturbation` now explicitly zero-initializes its `_vec_da` matrices. This is an apparent initialization fix relative to `main`, not an equivalence-preserving extraction.
- `OutputZiheng` and its checked-in radius-output data were added. The file has the expected three-field 359-by-720 grid shape; this structural check does not establish that the executable generated it or that its values are numerically validated.

## Corrections in this pass

- `ModelDensityOutputRotated` now scales density output with `DensityNorm()`. Both `main` (`4fd6488`) and the pre-correction PR implementation used `PotentialNorm()` here; this pass changes that behavior to the density scale. A dedicated regression uses distinct norms and verifies the expected density/Jacobian relation. The frozen stage-06 fixture uses unit scales and therefore could not expose this issue.
- A focused regression covers the discontinuous density trace, body-to-vacuum endpoint, all three rotated writers, first/final endpoints, and continuity of a continuous potential fixture.
- The published Phobos C++ fence now initializes rotation angles, uses declared solution/angle variables, creates its output directories, and compiles directly from the Markdown source.
- MathJax configuration now precedes its asynchronous loader.

## Validation boundaries

The original-source references and tolerances remain unchanged. The stage-06 fixture establishes cleanup equivalence to `4ef3a66` for its smooth model and unit norm scales; it does not establish main-to-head numerical equivalence at density jumps. The newly fixed density normalization intentionally changes density-output scaling where the two norms differ. No frozen numerical references were regenerated or tolerances relaxed. Experimental sources remain untouched and numerically unvalidated. Full all-target, complete CTest, stage-00, and public-header validation results are recorded in [pr-07-corrections.md](pr-07-corrections.md).
