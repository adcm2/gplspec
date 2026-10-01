# Deferred numerical observations

No suspected numerical defect was confirmed during stage 00, and no numerical behavior was changed.

The frozen baseline reports zero iterations for the homogeneous sphere, layered sphere, and lateral-density solve at relative tolerance `1e-6`; the smooth mapping solve reports one iteration, the mapping-perturbation base solve one, and the perturbation solve zero. The force and direct operator vectors are finite and nonzero in these fixtures, and all returned values are finite. These iteration counts are preserved observations; assess any concern separately before proposing a numerical correction.

No other numerical issues were identified in the stage-00 characterization.

## Non-blocking build provenance follow-up

The stage-00 manifest records the exact FFTW 3.3.10 and NetCDF 4.9.2 system packages used, but the current CMake `find_package` calls do not enforce those versions. The independent reviewer confirmed this as non-blocking for the fixed baseline; consider explicit version constraints in a future build-provenance task.

## Stage-05 operator characterization and coverage limits

The original-source stage-05 probes retain small nonzero differences between the complex sesquilinear forms; no symmetry assumption was imposed and no operator correction was made. The 1D values are `<x,Ay> = 302694.90083236271 - 302420.33613110462i` and `<Ax,y> = 302694.86333236267 - 302420.37363110413i`. The 3D file-perturbation values are `<x,Ay> = 1967.8406259289259 - 1923.06562591585i` and `<Ax,y> = 1967.8406259037456 - 1923.0656259410546i`. They are recorded observations from the original source, not newly diagnosed defects.

The stage-05 operator evidence covers finite deterministic fixtures, including a degree-2 coefficient perturbation. It does not establish exhaustive behavior for arbitrary high degrees, all models, or all platforms. The Sol review found this limitation non-blocking for the cleanup scope.
