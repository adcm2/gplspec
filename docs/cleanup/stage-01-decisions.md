# Stage-01 deletion and interface decisions

## Candidate classification

| Candidate | Classification | Evidence and decision |
|---|---|---|
| `gplspec/src/testtools.h` | Active implementation | Included by active 1D/3D model and perturbation headers; keep. |
| `gplspec/src/Mapping_Tools.h` | Public compatibility interface and active implementation | Included by `gplspec/Core` and `gplspec/All`; keep its include path and declarations. |
| `Pseudospectral_Matrix_Wrapper.h`, `Pseudospectral_Matrix_Wrapper3D.h` | Active operator implementation | Included by active gravity tools / umbrellas; keep. |
| `SphericalGeometryPreconditioner.h`, `Spherical_Integrator.h` | Active public support and independent numerical reference routines | Exposed by `Core`/`All` and used by the solver path; keep older/reference algorithms. |
| `testremove/Aspherical_Geometry_Potential.cpp`, `Spherical_Geometry_Potential.cpp` | Independent numerical reference programs | Contain separate aspherical and spherical potential calculations; retain as historical numerical references. |
| `testremove/Test.cpp`, `Spectral_Element.cpp`, `Spectral_Element_Sparse.cpp`, `Spectral_Element_Sparse_Test.cpp` | Independent numerical/reference experiments | Contain standalone operator/element formulations and are not in current CMake targets; retain for comparison/history. |
| `testremove/Eigen_Matrix_Wrapper.cpp` | Unused experiment | Standalone Eigen adapter demo; no GPLSpec include, CMake target, or example usage. Remove. |
| `testremove/Learning_GSHTrans.cpp` | Unused experiment | Interactive/random GSHTrans learning driver; no project target or example usage. Remove. |
| `testremove/crtptest.cpp` | Unused experiment | Standalone CRTP prototype with no GPLSpec usage or target. Remove. |
| `testremove/ellipticitytools.h`, `testremove/TestEllipticity` | Unused experiment/stub | Only reference each other; the solver class is empty and neither path is exposed by an umbrella or target. Remove both. |

The `testremove` directory is not an include path in any public umbrella and none of its programs are CMake targets. This repository-local evidence was combined with each file's contents; independent numerical implementations were retained instead of being treated as dead code.

## Other hygiene decisions

- Removed the tracked root CMake cache/build metadata because its cache points to a different checkout; see `implementation_status.md`.
- Keep solver iteration/error reports and constructor mesh/radius information. Remove only explicit debug traces such as numbered `Check` markers, temporary rotation-angle dumps, and row-dimension probes.
- Remove only commented-out umbrella includes that name internal implementation fragments; preserve equations, units, conventions, and explanatory comments.
