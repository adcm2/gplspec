# Field-bundle plotting

`plot_fields.py` is the shared consumer for schema-version-1 field bundles. Run
it from the repository root with `python3 work/plot_fields.py`. Each command
accepts one or more bundle directories and writes a PNG with `--output PATH`.
The reader validates the manifest and streams the CSV; plot payloads retain only
the selected numeric columns.

Examples:

```sh
python3 work/plot_fields.py radial RUN/work/Bench6/lm/harmonic_bundle \
  --quantity potential --representation harmonic_coefficients --degree 0 --order 0
python3 work/plot_fields.py radial RUN/work/Bench6/lm/sensitivity_bundle \
  --quantity sensitivity_kernel --representation harmonic_coefficients --degree 1 --order 0
python3 work/plot_fields.py map RUN/work/Phobos/Spherical/physical_rotated_grid_bundle \
  --quantity potential --element 9 --side inner
python3 work/plot_fields.py slice RUN/work/Phobos/Spherical/physical_rotated_slice_bundle \
  --quantity potential
python3 work/plot_fields.py slice RUN/work/Phobos/Spherical/rotated_density_slice_bundle \
  --quantity density_physical
python3 work/plot_fields.py power RUN/work/Phobos/Spherical/field_bundle \
  --element 10 --side outer
python3 work/plot_fields.py topography RUN/work/Ziheng/PhobosRadius.bundle
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
| `PowerSpectrum.py` | `power --element N --side inner|outer` selects a named endpoint trace from the degree-power representation. |

`field_bundle.py` remains the strict reader and `run_field_bundle.py` remains the
producer entry point. Plotting does not regenerate scientific data or alter
saved workflow products.
