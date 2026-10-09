# GPLSpec reconstruction — Phase 4

## Authorization and baseline

The user authorized all Phase 4 stages 4A–4F after accepting Phase 3, including its documented bench-10 grid-lifetime, homogeneous-sphere exterior, and normalization corrections. Phase 3's handoff remains historical. Work on `refactor/core-cleanup` from exact baseline `a8b1ca99fdadd9541a4227129c8fce872c98c1aa`; the checkout was clean at entry. The root coordinator archived the baseline at `/tmp/gplspec-phase4/baseline` and recorded its identity and 97 protected dataset/output/reference/experimental paths in `/tmp/gplspec-phase4/baseline-identity.json`. The tracked plotting sources in `work/*.py` are authorized to migrate; saved products remain protected. Existing Phase 3 final-build/workflow outputs may be used only as provenance-bounded baseline evidence. Do not modify the Phase 3 build or captures. New builds and workflow outputs use isolated Phase 4 directories.

Luna is the sole writer. Freeze each stage for separate independent Luna review; the coordinator commits accepted stages locally and advances automatically. At the Phase 4 boundary, obtain cumulative Sol-high review, present the cumulative diff/evidence, and stop for the user's human review. No push, merge, tag, branch deletion, or Phase 5 work is authorized.

## 4A — inventory and proposed minimal design

The inventory is from active source, tracked benchmark/plot files, and read-only references to protected experiments. It does not validate an experimental workflow or authorize touching saved outputs.

| Producer/capability | Existing scientific role and products | Active consumers / migration mapping |
|---|---|---|
| `Density3D` output methods in `gplspec/src/GeneralModels/Density_Model_Return.h` | Element radial profiles; referential and mapped physical potential samples; Cartesian coordinates; equatorial/rotated slices; referential and physical density; sensitivity-kernel samples; exact/integral spherical comparison values; normalized `PowerSTD`. | `examples/clean_bench_1..10.cpp`, the three Phobos examples and `OutputZiheng.cpp`; corresponding `tests/stage00_reference.cpp` and `tests/stage06_output_reference.cpp`. Replace writer calls with explicit evaluation requests and a serializer; retain mathematical solve/integral operations. |
| Direct `OutputZiheng.cpp` writers | `PhobosRadius.out` is a latitude/longitude radius surface; `PHOBOS_POTENTIAL_OUTER.out` is normalized complex `(l,m)` potential coefficients at the outer node; `PREM200.NO.INTEGRAL.out` is a normalized monopole radial integral profile. | Dedicated paper-verification consumer/input path; retain all three quantities and provenance, migrate to declared record representations. |
| clean_bench plotting | `CleanBench1Plot.py` through `CleanBench10Plot.py`, plus `CleanBench6Plot1.py`, `CleanBench6Plot_SK.py`, `CleanBench7Plot1.py`. | Bench 1 compares matrix/integral/exact sphere; 2 matrix/integral PREM; 3–5 mapped matrix vs analytic sphere; 6 matrix/integral and selected sensitivity; 7–8 mapped/unmapped cases and a radius sweep; 9 classical perturbation; 10 exact/classical/mapping perturbations. Migrate to bundle reader and preserve each plot's benchmark meaning. Exact model inputs, prompts, and numerical settings remain those catalogued in [`benchmark.md`](../../benchmark.md). |
| Phobos analysis/plotting | `PhobosFullSlice.py`, `PhobosPerturbation.py`, `PhobosSensitivityKernel.py`, `PhobosSensitivityKernel2.py`, `PhobosSliceDensity.py`, `PhobosSlicePlot.py`, `PhobosSlicePlot2.py`, `PhobosSurfacePotential.py`, `PhobosSurfacePotentialRotated.py`, `PhobosSurfaceTopography.py`, and `PowerSpectrum.py`. | Preserve gravity maps/slices, rotated potential, referential/physical density contrasts, low-degree perturbation comparisons, sensitivity, surface/topography and degree-power plots. `PowerSTD` is `sum_m |Phi_lm|^2 / |Phi_00|^2` per degree; its output must remain explicitly normalized power. |
| Protected experimental include chain | Six files under `experimental/`; `ellipticity_test.cpp` includes `ellipticitytools.h`, which includes `Spherical_Integrator.h`; it uses `GravityFunctions`, `Radial_Node`, and `All_Node`. | No protected writer or plotting call was found. Preserve only the existing small mesh-utility compatibility surface if output/API retirement would otherwise break this named include. Leave every experimental file unchanged, unbuilt and unvalidated. |

### Proposed separation and schema

Use one small, stateless field evaluator, a serializer that performs no scientific calculation, and plotting readers that perform no field evaluation. The evaluator exposes requested samples through a callback and keeps only bounded per-node angular workspaces; it must not materialize or copy the complete field merely to write it. Radial element/GSH grid transforms are batched and reused; rotated coefficients are prepared once per requested rotation, not once per sample. Phobos `lmax=100` workloads contain roughly 300,000 angular samples: do not replace batched inverse transforms with per-point harmonic summation or a full basis cache. Arbitrary off-grid angular requests, if needed by a real consumer, use a separately bounded path and are not the default grid/slice route.

Proposed bundle: a compact `manifest.json` plus a long-form `samples.csv`. The manifest records schema version, workflow/input hashes, physical parameters, solver method/tolerance/diagnostics, normalization scales, GSH degree/order/index convention, and a per-quantity definition, units, frame, and representation. A sample row carries quantity/representation, numerical element-boundary ID, optional material-interface ID, adjacent element/material IDs, side, reference and physical coordinates when defined, optional harmonic `(l,m)`, and real/imaginary value columns. Empty coordinates/components are explicitly inapplicable, not implicit zeros. A harmonic coefficient, spatial sample, radial profile sample, and degree-power datum are distinguished by the representation field. A schema reader rejects unknown versions/quantities, missing required fields, non-finite values, and duplicate sample keys.

Interface contract:

- Every element boundary has a numerical-interface ID. A material-interface ID exists only where element-local material IDs differ or at the declared body/vacuum boundary; numerical interfaces do not imply density jumps.
- Emit one row per adjacent element trace, with `side=inner|outer` defined by increasing **reference** radius and stable adjacent element/material IDs. Preserve each element's own density sample. This includes both body and vacuum traces whenever an exterior vacuum element is part of the model.
- Continuous potential coefficients use the shared global degree of freedom. They may be represented on both trace rows for association, but must carry the same value; never create a jump by independently interpolating the two sides.
- Always distinguish reference radius from mapped physical position. At mapped interfaces, physical position may depend on angle. A model with no exterior element has no modeled outer vacuum trace: emit its available inner trace and explicit absence metadata; do not extrapolate or invent an outer sample.
- Plot profiles as element/material segments, retaining both interface points but not joining across a material jump. Spatial interpolation is confined to one element at a time; mapped interfaces use physical coordinates while retaining reference-coordinate identity.

Quantity contract to implement and independently check:

| Quantity | Meaning, units and normalization |
|---|---|
| Potential field | Real spatial gravitational potential in `m^2 s^-2`; normalized harmonic coefficients use the pinned GSHTrans convention and are scaled by the model `PotentialNorm`. Reference and mapped physical coordinates are separately labeled. |
| Density | `density_referential` is mass per reference volume in `kg m^-3`; `density_physical` is mass per mapped physical volume in `kg m^-3`, related at samples by `rho_physical=rho_referential/J`. Keep them separate in manifest and rows. |
| Sensitivity | A selector-specific derivative kernel for the requested potential harmonic with respect to referential source density integrated over reference volume, using the existing `PotentialNorm/MassNorm` scale. Its units are potential per mass (`m^2 s^-2 kg^-1`); record reference integration measure and selector explicitly. It is not itself a potential field. A physical-density sensitivity view is emitted only where the Jacobian conversion has an independent contraction check. The accepted 3E.4 homogeneous density perturbation is the reference contract. |
| Harmonic outer potential | Complex coefficient by `(l,m)` at the named outer node, in physical potential units and with GSH normalization; used by the OutputZiheng comparison. |
| Degree power | Existing `PowerSTD` definition above, dimensionless normalized power by degree and radial sample; not potential amplitude. |
| Radius/topography | Physical surface radius in metres at explicit latitude/longitude in radians; preserve the Phobos TABLEA1 input harmonics and reference/physical coordinate distinction. |

### Predeclared comparisons and validation

Freeze acceptance scales before evaluating 4B replacements. Match samples by quantity, reference location, numerical/material interface and side, never by row order or filenames.

- Coordinates and material/interface identities: exact IDs and trace counts; normalized coordinate `rtol=1e-12`, `atol=1e-14*max(1,|x_reference|)` in the corresponding frame. Physical coordinates use the same formula after scaling to metres.
- Density values: `rtol=1e-12`, `atol=1e-14*rho_norm`; independently check `rho_physical*J=rho_referential` at non-unit Jacobian and both one-sided interface traces. Layer and body/vacuum IDs are exact.
- Potential samples: compare at matched points using `|delta Phi| <= 1e-10*|Phi_ref| + 1e-12*G*M/R`, where `M,R` are the fixture's physical mass and outer-radius scales. For harmonic coefficients, also apply the addition-theorem bound to the full retained complex coefficient vector; do not use component-count mismatches as an acceptance rule.
- Sensitivity: use the accepted 3E.4 independent density-perturbation contraction/finite-difference as the scientific check; representation migration requires full sampled-kernel `||delta K||2 <= 1e-8*||K_ref||2`. For an individual near-zero coefficient use `1e-8*||K_ref||2` as its absolute comparison floor. Do not apply potential tolerances to kernels.
- Power: compare `sqrt(P_l)` relative to the monopole amplitude at `1e-8`; near-zero degrees use an absolute `1e-8` amplitude floor in that same normalized scale. Retain the distinction between power and amplitude.
- Serialization: round-trip all metadata/typed rows, reject malformed/non-finite/duplicate records, preserve one-sided trace IDs and equal shared-potential values. Keep frozen Phase 3 references unchanged; exact bytewise comparisons are diagnostic only.
- Rotation/coordinates: verify identity and known axis rotations against independent Cartesian/spherical coordinate identities; compare values at corresponding transformed points and same trace side. Plot checks include a layered discontinuity and mapped Phobos surface/slice.

These are criteria for output/evaluation reconstruction at corresponding discrete points; they do not impose a new solver-accuracy threshold or change workflow settings. If an established baseline comparison and an independent identity disagree, diagnose the specific quantity and geometry before deciding acceptance; do not widen these bounds.

The Phase 4F performance comparison uses the root-captured baseline cohort and separately times solve-only, evaluation/export (including each side's required flat-to-nested or equivalent conversion), and complete solve-to-output. Keep the legacy cohort of angular grid points at each element's lower endpoint plus the final outer endpoint matched exactly on both sides. Count any additional mandatory interface-trace rows/work separately; do not hide that work in a changed sample cohort or attribute it to an implementation speedup. Include repeated one-thread measurements and peak memory.

### Stage sequence and retirement

| Stage | Scope and retirement point | State |
|---|---|---|
| 4A | Inventory and agree this interface/schema, trace policy, comparison criteria, workflow map, and protected-consumer limit. | Accepted by independent Luna with no blockers. No production implementation is part of this checkpoint. |
| 4B | Implement batched field sampling for radial profiles, angular grid samples, and rotated slices; explicit numerical/material IDs and both modeled one-sided traces, including body/vacuum. Add independent density/Jacobian, shared-potential, coordinate/rotation and surface-only tests. Existing writers remain temporary adapters. | Next after the coordinator commits 4A; not started. |
| 4C | Add manifest/CSV serialization and round-trip tests; migrate public output calls/examples, OutputZiheng, and reference capture producers. Keep the current writers temporarily producing the old products while 4D consumers still require them; the adapter-removal point is 4E. | Not started. |
| 4D | Migrate tracked `work/*.py` plot/postprocessing consumers; segment radial curves, avoid cross-interface interpolation, relabel units/frames, visually inspect representative layered and Phobos plots in temporary directories. Retire old file readers after their last consumer migrates. | Not started. |
| 4E | Consolidate clean_bench/Phobos driver and tutorial workflows, document solve-to-plot reproduction and paper-label limits, remove superseded active output code/comments, and retire old writers/adapters once every plot consumer is migrated. Preserve model inputs, parameters, solver settings, benchmark roles and supported quantities. Retain only necessary protected experimental compatibility. | Not started. |
| 4F | Fresh Release build/CTest/header consumers, required interface/normalization/rotation/round-trip/plot checks, migrated Phase 3 workflows and physical comparisons, representative solve/evaluate/export/full-path performance and memory, cumulative Sol-high review. | Not started. |

Every stage records actual source changes and fresh validation in `implementation_status.md`; accepted review verdicts are added to `phase-4-reviews.md`. Preserve all 97 protected files and frozen references; never write outputs into saved `work/` or dataset directories.
