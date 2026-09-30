# Deferred numerical observations

No suspected numerical defect was confirmed during stage 00, and no numerical behavior was changed.

The frozen baseline reports zero iterations for the homogeneous sphere, layered sphere, and lateral-density solve at relative tolerance `1e-6`; the smooth mapping solve reports one iteration, the mapping-perturbation base solve one, and the perturbation solve zero. The force and direct operator vectors are finite and nonzero in these fixtures, and all returned values are finite. These iteration counts are preserved observations; assess any concern separately before proposing a numerical correction.

No other numerical issues were identified in the stage-00 characterization.

## Non-blocking build provenance follow-up

The stage-00 manifest records the exact FFTW 3.3.10 and NetCDF 4.9.2 system packages used, but the current CMake `find_package` calls do not enforce those versions. The independent reviewer confirmed this as non-blocking for the fixed baseline; consider explicit version constraints in a future build-provenance task.
