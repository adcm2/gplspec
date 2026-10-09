# Field-bundle plotting

`plot_fields.py` is the shared consumer for schema-version-1 field bundles. Run
it from the repository root with `python3 work/plot_fields.py`. Each command
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
python3 work/plot_fields.py radial "$BENCH6_RUN/work/Bench6/lm/harmonic_bundle" \
  --quantity potential --representation harmonic_coefficients --degree 0 --order 0
python3 work/plot_fields.py radial "$BENCH6_RUN/work/Bench6/lm/sensitivity_bundle" \
  --quantity sensitivity_kernel --representation harmonic_coefficients --degree 1 --order 0
python3 work/plot_fields.py map "$PHOBOS_RUN/work/Phobos/Spherical/physical_rotated_grid_bundle" \
  --quantity potential --element 9 --side inner
python3 work/plot_fields.py slice "$PHOBOS_RUN/work/Phobos/Spherical/physical_rotated_slice_bundle" \
  --quantity potential
python3 work/plot_fields.py slice "$PHOBOS_RUN/work/Phobos/Spherical/rotated_density_slice_bundle" \
  --quantity density_physical
python3 work/plot_fields.py power "$PHOBOS_RUN/work/Phobos/Spherical/field_bundle" \
  --element 10 --radial-node 0
python3 work/plot_fields.py topography "$ZIHENG_RUN/work/Ziheng/PhobosRadius.bundle"
```

For the Bench 6 harmonic and sensitivity plots, the producer requires both
unchanged PREM and S40RTS inputs. Its recorded settings are `npoly=8`,
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
recorded workflow settings (`npoly=5`, `lmax=2`, `maxstep=0.01`,
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
python3 work/plot_fields.py radial "$RUN_DIR/work/Bench2/field_bundle" \
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
python3 work/plot_fields.py map \
  "$PHOBOS_RUN/work/Phobos/Spherical/physical_rotated_grid_bundle" \
  --quantity potential --element 9 --side inner \
  --output /tmp/phobos-potential.png
```

Radial plots retain both one-sided endpoint samples and segment lines by
element/material, so they do not draw ramps across a discontinuity. Physical
radius is used when present; harmonic and reference-sphere profiles use
reference radius. Harmonic plots label the selected real or imaginary `(l,m)`
coefficient and reference-radius axis; use `--component real` or
`--component imag` (the default is `real`). Angular/radial profiles label the
sampled ray in radians and physical radius where present. `--difference` adds an absolute-difference panel only after
element/material, side, reference-radius, and physical-coordinate agreement is
checked (coordinate tolerance `rtol=1e-12`, `atol=1e-14`); omit it for
perturbation geometries that intentionally differ. `auto` can compare sampled
angular grids with exact radial profiles on a shared ray, or harmonic
coefficients with harmonic coefficients. Different scientific representations
should be plotted in separate invocations.

Maps require an explicit element and trace side. They show coordinates in the
bundle's recorded frame: a rotated output frame when rotation angles are
present, otherwise the body-fixed frame. Slice interpolation uses only the two
boundary rings belonging to each element, with a center fan for the central
element. No interpolation crosses a material interface. Axes and colorbars
carry SI units from the bundle contract.

## Retired plot scripts and replacement

The old scripts read positional columns from legacy `.out` files. Their plot
roles now use field-bundle quantities and explicit element/ray selectors:

| Retired script | Replacement |
| --- | --- |
| `CleanBench1Plot.py` | `radial RUN/work/Bench1/field_bundle RUN/work/Bench1/integral_bundle RUN/work/Bench1/exact_bundle --quantity potential --representation harmonic_coefficients --degree 0 --order 0`. |
| `CleanBench2Plot.py` | `radial RUN/work/Bench2/field_bundle RUN/work/Bench2/integral_bundle --quantity potential --representation harmonic_coefficients --degree 0 --order 0`. |
| `CleanBench3Plot.py`, `CleanBench4Plot.py`, `CleanBench5Plot.py` | `radial` potential on a shared ray; `auto` supports exact `radial_profile` against a field `angular_grid`. Use `--difference` only for matching geometry. |
| `CleanBench6Plot.py` | `radial RUN/work/Bench6/lm/harmonic_bundle RUN/work/Bench6/lm/integral_bundle --quantity potential --representation harmonic_coefficients --degree L --order M`. |
| `CleanBench6Plot1.py` | The potential comparison uses the `CleanBench6Plot.py` command. For each sensitivity `(L,M)`, run `radial RUN/work/Bench6/lm/sensitivity_bundle --quantity sensitivity_kernel --representation harmonic_coefficients --degree L --order M --component real` and repeat with `--component imag`. |
| `CleanBench6Plot_SK.py` | For each sensitivity `(L,M)`, run `radial RUN/work/Bench6/lm/sensitivity_bundle --quantity sensitivity_kernel --representation harmonic_coefficients --degree L --order M --component real` and repeat with `--component imag`. |
| `CleanBench7Plot.py`, `CleanBench8Plot.py`, `CleanBench9Plot.py`, `CleanBench10Plot.py` | `radial` potential overlays for exact, perturbed, and unperturbed bundles. Omit `--difference` when mappings differ. |
| `CleanBench7Plot1.py` | `radial` potential overlay for zero-radius and finite-radius perturbation runs; compare the shared reference geometry only. |
| `PhobosFullSlice.py`, `PhobosSlicePlot.py`, `PhobosSlicePlot2.py` | `slice --quantity potential`, `density_physical`, `density_referential`, or `sensitivity_kernel`; each bundle retains its distinct quantity/units. |
| `PhobosSliceDensity.py` | `slice --quantity density_physical` or `density_referential`. |
| `PhobosPerturbation.py` | `map` for sampled angular-grid comparisons, or `slice` for equatorial comparison; each perturbation geometry remains a separate field bundle. |
| `PhobosSensitivityKernel.py` | `radial` or `map --quantity sensitivity_kernel`; radial plots select `--degree/--order`. |
| `PhobosSensitivityKernel2.py` | `slice --quantity sensitivity_kernel`. |
| `PhobosSurfacePotential.py`, `PhobosSurfacePotentialRotated.py` | `map --quantity potential --element N --side inner|outer`; run once per requested frame/trace. |
| `PhobosSurfaceTopography.py` | `topography RUN/work/Ziheng/PhobosRadius.bundle`, which carries the explicit `surface_radius` quantity. |
| `PowerSpectrum.py` | `power --element N --radial-node J` selects one element-local radial node from the degree-power representation; all local nodes are emitted, including both interface traces. |

`field_bundle.py` remains the strict reader and `run_field_bundle.py` remains the
producer entry point. Plotting does not regenerate scientific data or alter
saved workflow products.

Legacy plot-script names identify migrated workflow roles only; this mapping
does not assert a specific paper figure or publication provenance.
