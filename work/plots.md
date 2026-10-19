# Field-bundle plotting

`plot_fields.py` is the shared consumer for schema-version-1 field bundles. Run
it from the repository root with the test/plot virtual environment's interpreter (see [getting started](../getting-started.md) for setup), for example `build/python-test-env/bin/python work/plot_fields.py`. Each command
accepts one or more bundle directories and writes a PNG with `--output PATH`.
The reader validates the manifest and streams the CSV; plot payloads retain only
the selected numeric columns.

For a producer, build from the repository root but run it through
`run_field_bundle.py` with an isolated `/tmp` working directory. The wrapper
binds each bundle to the named unchanged input; use absolute paths for the
executable and keep generated products out of the repository's saved `work/`
tree. Plot commands below read the corresponding isolated `RUN_DIR/work/...`
paths.

Set these to the isolated runs that produced the bundles you want to inspect;
the PREM and Phobos producer recipes below create the corresponding run dirs.

```sh
BENCH6_RUN=/tmp/gplspec-bench6-bundle-run
PHOBOS_RUN=/tmp/gplspec-phobos-bundle-run
ZIHENG_RUN=/tmp/gplspec-ziheng-bundle-run
```

Examples:

```sh
build/python-test-env/bin/python work/plot_fields.py radial "$BENCH6_RUN/work/Bench6/lm/harmonic_bundle" \
  --quantity potential --representation harmonic_coefficients --degree 0 --order 0
build/python-test-env/bin/python work/plot_fields.py radial "$BENCH6_RUN/work/Bench6/lm/sensitivity_bundle" \
  --quantity sensitivity_kernel --representation harmonic_coefficients --degree 1 --order 0
build/python-test-env/bin/python work/plot_fields.py map "$PHOBOS_RUN/work/Phobos/Spherical/physical_rotated_grid_bundle" \
  --quantity potential --element 9 --side inner
build/python-test-env/bin/python work/plot_fields.py slice "$PHOBOS_RUN/work/Phobos/Spherical/physical_rotated_slice_bundle" \
  --quantity potential
build/python-test-env/bin/python work/plot_fields.py slice "$PHOBOS_RUN/work/Phobos/Spherical/rotated_density_slice_bundle" \
  --quantity density_physical
build/python-test-env/bin/python work/plot_fields.py power "$PHOBOS_RUN/work/Phobos/Spherical/field_bundle" \
  --element 10 --radial-node 0
build/python-test-env/bin/python work/plot_fields.py topography "$ZIHENG_RUN/work/Ziheng/PhobosRadius.bundle"
```

For the Bench 6 harmonic and sensitivity plots, the producer requires both
unchanged PREM and S40RTS inputs. Its recorded settings are `npoly=12`,
`lmax=20`, `maxstep=0.05`, and `ballrad=1.2`:

```sh
cmake --build build --target clean_bench_6
mkdir -p "$BENCH6_RUN/modeldata"
ln -s "$PWD/modeldata/prem.200" "$BENCH6_RUN/modeldata/prem.200"
ln -s "$PWD/modeldata/S40RTS_dvs.nc" "$BENCH6_RUN/modeldata/S40RTS_dvs.nc"
python3 work/run_field_bundle.py --workdir "$BENCH6_RUN" \
  --input modeldata/prem.200 --input modeldata/S40RTS_dvs.nc \
  --bundle work/Bench6/field_bundle \
  --bundle work/Bench6/lm/integral_bundle \
  --bundle work/Bench6/lm/harmonic_bundle \
  --bundle work/Bench6/lm/sensitivity_bundle -- "$PWD/build/bin/clean_bench_6"
```

To reproduce the PREM comparison from unchanged `modeldata/prem.200` using the
recorded workflow settings (`npoly=5`, `lmax=2`, `maxstep=0.1`,
`ballrad=1.2`), configure the producer, link the input into an isolated run
directory, and invoke its absolute executable path:

```sh
cmake -S . -B build -DMY_PROJECT_BUILD_EXAMPLES=ON
cmake --build build --target clean_bench_2
RUN_DIR=/tmp/gplspec-prem-bundle-run
mkdir -p "$RUN_DIR/modeldata"
ln -s "$PWD/modeldata/prem.200" "$RUN_DIR/modeldata/prem.200"
python3 work/run_field_bundle.py --workdir "$RUN_DIR" \
  --input modeldata/prem.200 --bundle work/Bench2/field_bundle \
  --bundle work/Bench2/integral_bundle -- "$PWD/build/bin/clean_bench_2"
build/python-test-env/bin/python work/plot_fields.py radial "$RUN_DIR/work/Bench2/field_bundle" \
  "$RUN_DIR/work/Bench2/integral_bundle" --quantity potential \
  --representation harmonic_coefficients --degree 0 --order 0 \
  --output /tmp/prem-potential.png
```

For Phobos, build `phobos_gravity` and run it in a separate isolated directory
with a symlink to unchanged `modeldata/TABLEA1.DAT` and settings recorded in each manifest
(`density=1860 kg m^-3`, `length_norm=10000 m`, `npoly=5`, `lmax=100`,
`max_step=0.1`, `exterior_radius=1.5`; solver tolerance `1e-12`). The same plot
commands above can then read the isolated `work/Phobos/Spherical/field_bundle`,
`physical_rotated_grid_bundle`, and `physical_rotated_slice_bundle` products:

```sh
cmake --build build --target phobos_gravity
PHOBOS_RUN=/tmp/gplspec-phobos-bundle-run
mkdir -p "$PHOBOS_RUN/modeldata"
ln -s "$PWD/modeldata/TABLEA1.DAT" "$PHOBOS_RUN/modeldata/TABLEA1.DAT"
python3 work/run_field_bundle.py --workdir "$PHOBOS_RUN" \
  --input modeldata/TABLEA1.DAT \
  --bundle work/Phobos/Spherical/field_bundle \
  --bundle work/Phobos/Spherical/physical_rotated_grid_bundle \
  --bundle work/Phobos/Spherical/equatorial_slice_bundle \
  --bundle work/Phobos/Spherical/physical_rotated_slice_bundle \
  --bundle work/Phobos/Spherical/referential_rotated_grid_bundle \
  --bundle work/Phobos/Spherical/referential_rotated_slice_bundle \
  --bundle work/Phobos/Spherical/rotated_density_slice_bundle \
  --bundle work/Phobos/Sensitivity/field_bundle \
  --bundle work/Phobos/Sensitivity/sensitivity_slice_bundle \
  -- "$PWD/build/bin/phobos_gravity"
build/python-test-env/bin/python work/plot_fields.py map \
  "$PHOBOS_RUN/work/Phobos/Spherical/physical_rotated_grid_bundle" \
  --quantity potential --element 9 --side inner \
  --output /tmp/phobos-potential.png
```

Radial plots retain both one-sided endpoint samples and segment lines by
element/material, so they do not draw ramps across a discontinuity. Physical
radius is used when present; harmonic and reference-sphere profiles use
reference radius. Harmonic plots label the selected real or imaginary `(l,m)`
coefficient and reference-radius axis; use `--component real` or
`--component imag` (the default is `real`). Profiles with recorded angles label their sampled ray in radians. Angle-free radial profiles use their declared radial coordinates and do not claim a ray. Physical radius is labeled where present. `--difference` adds an absolute-difference panel only after
element/material, side, reference-radius, and physical-coordinate agreement is
checked (coordinate tolerance `rtol=1e-12`, `atol=1e-14`); manifests must
also have matching coordinate and rotation frames. Omit it for perturbation
geometries that intentionally differ. `auto` can compare sampled angular grids
with finite-angle radial profiles on a shared ray, or harmonic coefficients
with harmonic coefficients. Angle-free radial profiles remain valid for
standalone radial plots at their declared coordinates and do not claim a ray;
differences still require matching bundle frames and radial coordinates.
Different scientific representations should be plotted in separate invocations.

Maps require an explicit element and trace side. They show coordinates in the
bundle's recorded frame: a rotated output frame when rotation angles are
present, otherwise the body-fixed frame. A `surface_radius` slice draws only
the exported surface outline in its output frame; it does not infer an interior
mesh. The reader accepts existing angular-grid-only surface manifests; new
manifests declare both angular-grid and slice forms. Slice interpolation for
field quantities uses only the two boundary rings belonging to each element,
with a center fan for the central
element. No interpolation crosses a material interface. Axes and colorbars
carry SI units from the bundle contract.

## CleanBench 7 zero-radius sweep

This synthetic workflow uses a unit-radius, unit-density `HomogeneousSphere` and
`TomographyZeroModel`; it does not consume PREM or S40RTS files. The producer
retains its original `npoly=5`, `lmax=20`, `maxstep=0.1`, `ballrad=1.5`,
`h` prompt, and 100 diagnostic amplitudes `h=0, 0.002, ..., 0.198`. The output
table columns are `h;mapping-perturbation centre potential;analytic exact
centre potential`, with potential in m² s⁻². Bundle metadata records the
synthetic model and sphere parameters.

Run it in an isolated directory, then plot the saved diagnostic table with the
shared CLI:

```sh
cmake --build build --target clean_bench_7
BENCH7_RUN=/tmp/gplspec-bench7-sweep-run
mkdir -p "$BENCH7_RUN"
printf '0.001\n' | python3 work/run_field_bundle.py --workdir "$BENCH7_RUN" \
  --bundle work/Bench7/exact/field_bundle \
  --bundle work/Bench7/perturbed/field_bundle \
  --bundle work/Bench7/unperturbed/field_bundle -- "$PWD/build/bin/clean_bench_7"
build/python-test-env/bin/python work/plot_fields.py sweep \
  "$BENCH7_RUN/work/Bench7/zeroradius/MatrixSolution.out" \
  --output /tmp/bench7-zero-radius-sweep.png
```

The left panel keeps the zero-amplitude sample and labels the analytic and
mapping-perturbation curves correctly. The right panel plots
`abs(exact-approximate)/max(abs(approximate))*100` and the original quadratic
guide normalized at the fourth sample (zero-based row index 3). Zero `h` and
zero-error points are explicitly omitted from the logarithmic panel. The CLI
requires a finite, three-column table, strictly increasing nonnegative `h`, a
nonzero approximate-potential scale, and a nonzero error at its positive guide
anchor.

## Current plot selection reference

Use these selectors with the isolated bundle paths and producer commands above. The supported build settings and output products are recorded in [scientific workflows](../benchmark.md).

| Workflow products | Plot command and selectors |
|---|---|
| CleanBench 1 and 2 field/integral bundles | `radial` with `--quantity potential --representation harmonic_coefficients --degree 0 --order 0`; include the Bench 1 exact bundle for the analytic overlay. |
| CleanBench 3–5 mapped field and exact bundles | `radial --quantity potential`; `auto` can pair finite-angle `radial_profile` and `angular_grid` data on a shared ray. Use `--difference` only when ray and coordinates agree. Angle-free radial profiles plot at their declared radii without ray labels. |
| CleanBench 6 harmonic and sensitivity bundles | Potential: select `(L,M)` with `radial --quantity potential --representation harmonic_coefficients --degree L --order M`. Sensitivity: use `--quantity sensitivity_kernel` with the same selectors and repeat with `--component real` and `--component imag`. |
| CleanBench 7–10 exact, perturbed, and unperturbed field bundles | `radial --quantity potential` overlays the available scalar radial profiles. Leave `--difference` off when mapping geometries or coordinates differ. CleanBench 7's three-column zero-radius diagnostic is plotted separately with `sweep`. |
| Phobos potential and density angular grids | `map --quantity potential`, `density_physical`, or `density_referential`; select `--element N --side inner|outer` for an explicit modeled trace. Use each bundle's frame and rotation metadata. |
| Phobos potential, density, and sensitivity slices | `slice --quantity potential`, `density_physical`, `density_referential`, or `sensitivity_kernel`; choose the bundle for the intended frame and workflow. |
| Phobos sensitivity radial/angular products | `radial` or `map --quantity sensitivity_kernel`; radial harmonics require `--degree L --order M`. |
| Phobos degree power | `power --element N --radial-node J` selects one element-local node; both interface traces are emitted. |
| Ziheng surface radius | `topography RUN/work/Ziheng/PhobosRadius.bundle`, whose quantity is `surface_radius`. |

These commands operate on schema-v1 bundles or the retained CleanBench 7 diagnostic table. They do not recreate retired `.out` files or assert publication-figure provenance.
