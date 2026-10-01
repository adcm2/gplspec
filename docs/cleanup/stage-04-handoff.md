# Stage 04 handoff

- Starting accepted base: `2d4f7dfa3d8ba4adb55a4cdc8f1ca0f0bbb8d501`
- Branch: `cleanup/04-models`
- Fixed code candidate: `96d4ee9321b0adf4c30b583aed98924707049ff7`
- Implementer: `gpt-6-luna` medium; independent reviewer: pending `gpt-6-luna` high.

Stage 04 extracted four bounded internal operations from `Density_Model_Constructors.h`: common quadrature/grid/order setup, nested scalar-field initialization, mapped-radius geometry fill, and referential tomography density variation. Each extraction is a separate commit. Public constructor signatures and output interfaces are unchanged. The mapping helper preserves clamp/taper order and per-sample mapping calls. The tomography helper preserves the caller's distinct `OuterRadius()` versus `PlanetRadius()` depth references. Physical density multiplied by the Jacobian remains separate; the aspherical transformed-radius geometry and constructor-specific `RadialMesh` calls remain separate.

Validation command: configure `/tmp/gplspec-stage02-guard-build` with `-DMY_PROJECT_BUILD_EXAMPLES=ON -DGPLSPEC_BUILD_BASELINE_HARNESS=ON`, then `cmake --build /tmp/gplspec-stage02-guard-build -j2`. This built four public-header smoke targets, `stage02_header_link`, all four focused stage04 executables, the stage00 harness, and all examples. All six CTests passed, and `stage02_header_link` ran. `tests/run_stage00.sh /tmp/gplspec-stage02-guard-build/bin/stage00_reference` passed three comparisons of 20,700 records, zero numerical changes, and byte-identical representative output. `git diff --check` passed.

Frozen hashes are unchanged: CSV `ef20bd7e3016661c60903c290c74c599151fdfc0417443bcf81c4215401f105e`, representative output `a902254d3d416c55f103058c045a30036e7b5edd8b56166e086d43de4de13c3f`, and durable original-source executable `456d69bb388149ced74dc5078fe825e3d79193f64ad254ef13898ff2baf03308`. No expected data was regenerated.

No numerical differences or deferred stage04 findings were observed. Independent Luna review is pending. Do not begin stage 05 until this candidate passes review and is integrated into `cleanup/base`.
