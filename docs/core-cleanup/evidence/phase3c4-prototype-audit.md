# Phase 3C.4 prototype and geometry audit

The source identities below are Git blob IDs from the accepted parent tree
`dcb1a96` (the prior 3C.3 checkpoint). The six files under
`gplspec/src/testremove/` are not CMake targets and are not built by the
supported workflows.

| Source (blob) | Retained scientific content and disposition |
|---|---|
| `Spectral_Element.cpp` (`e03ebdc94d2122e90253a10b11441c3e7475018f`) | Dense scalar radial Poisson weak-form experiment. The local assembly uses a radius-squared derivative contribution, a degree term `l(l+1)`, and exterior diagonal `-(l+1)R`; it also includes PREM mass/potential integrals. The exact normalization/coordinate reconciliation with the active weak form was not established in this unbuilt prototype. The active radial action is covered by the retained operator-action test; independent continuum integrals remain scheduled for 3E. |
| `Spectral_Element_Sparse.cpp` (`67e07c247c3b495e40f54b2acffef712493d3a18`) | Sparse triplet weak-form/solver experiment. Its radial derivative coefficient is `r_phys*r_scale`; its coordinate and normalization reconciliation with the dense experiment and active weak form was not established. The active mapped action is independently covered by the retained operator-action test. |
| `Spectral_Element_Sparse_Test.cpp` (`199df8d4becc6ead06c22805c9e9ee1df9c2a4c1`) | Second sparse assembly/solver experiment with preconditioner and solver alternatives. No independent validated scientific observable is defined. Retained solver and operator tests exercise the supported algorithm contracts. |
| `Spherical_Geometry_Potential.cpp` (`5f7f7bf5863cb41f860633ed80fe6e8d94c0d73b`) | PREM spherical-harmonic solve at `lmax=4`, polynomial order 3, 100 km maximum radial step, and outer radius 1.2. It contrasts the full harmonic solve/source with `PotentialSolver1D` and exports radial/gravity values. Current source quadrature, independent radial density sampling, and operator tests retain the reusable checks; integral comparisons are assigned to 3E. |
| `Aspherical_Geometry_Potential.cpp` (`b6095d385df5687719063c0c28adbc18bb5bbce5`) | PREM plus PERTPREM mapping and S40RTS tomography at `lmax=10`, polynomial order 5, 300 km maximum radial step, and outer radius 1.2. It forms mapped spin-gradient/tensor contributions and compares alternative operator/solver routes. Active mapped-action, tensor-contraction, and solver tests cover the supported scientific contracts; this prototype's own output was not validated as an independent reference. |
| `Test.cpp` (`1d8e3e435d669baf073cfc0c078033bdd8d777ba`) | PREM mapped potential experiment at `lmax=5`, polynomial order 3, 1000 km maximum radial step, and outer radius 1.2, including a manual 3x3 spin/tensor action and solver comparison. The retained direct tensor-contraction and mapped operator-action tests independently check the supported formulas; the prototype output was not adopted as a reference. |

The retired `Earth_General_Models_3D.h` has parent-tree blob
`3cd3489bcce051e5f62e41554380e6d06a2474dc`, retrievable as
`git show dcb1a96:gplspec/src/GeneralModels/Earth_General_Models_3D.h`.
The legacy `General3D` audit used the same `PREM500` input, four-point
Gauss-Lobatto quadrature, maximum element step 0.5 and outer radius 1.2 as the
retained radial model and `Density3D` route. Its node geometry, layer labels,
density and radial elastic samples agree with direct `RadialEarthModel`
evaluation at all element nodes under the predeclared `rtol=1e-12`,
`atol=1e-14*S_q` criterion. Adjacent elements preserve separate material
values at their shared-radius interface; the exterior layer remains zero.
`General3D` stores those radial properties replicated over its angular grid,
which can be reconstructed from `RadialEarthModel`/`spherical_1D`. Its mapping
is zero and its Laplace-tensor storage is zero-initialized; it does not provide
the physical mapped-tensor capability supplied by `Density3D`. No active
consumer or unique scientific capability requires this duplicate owner.

Focused candidate checks passed for the new geometry-ownership test, the
independent source-discretization test, and the operator-action test. The
frozen Stage 05 CSV remains unchanged. Its exact-output script reports small
operation-order differences: the complex `Ax` vector relative difference is
`7.32e-15`; the source `force` vector difference is `1.90e-16`; the boundary
delta difference is `7.13e-17`. The operator action remains inside its
predeclared `1e-11 + 1e-13*||A_ref*x||` bound, and source values pass the direct
quadrature check. These are diagnostics, not changes to the frozen reference
or the scientific tolerances.
