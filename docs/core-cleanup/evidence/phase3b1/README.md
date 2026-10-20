# Phase 3B.1 input and model evidence

This checkpoint separates PREM and Phobos file parsing from the current construction code. The baseline is `15582edab8112dbb63e0101be5b26f5135298644`; the candidate was built from that baseline plus the uncommitted 3B.1 changes. No model input, frozen reference, or saved workflow output was changed.

## Input and constructed-value comparison

`run-input-comparisons.sh` builds temporary probes against both the detached baseline and candidate Release dependency trees. It dumps all normalized PREM500 and `prem.200` layer radii and material interpolants, then compares every row byte-for-byte. Its Phobos probe independently extracts all six TABLEA1 fields from the legacy stream path or the typed reader and records the historical `A-iB` conversion, nonzero-order half factor, and degree-20 truncation. The probe source and the baseline raw result CSVs are retained here; candidate outputs had the same SHA-256 in all three comparisons.

| Fixture | Rows / checks | Baseline and candidate raw CSV SHA-256 |
|---|---|---|
| PREM500 | All normalized layer radii, interfaces, density, VPV/VSV/VPH/VSH, eta and solid flags | `170ec35ed193c5b92da6e9886f9de5473b752ccfc9f3fa165cd68b22dd7fc527` |
| prem.200 | Same full model-value dump for the clean_bench_6 input | `5c815d823a186f826c64aa294a615772f8707e5468c0a9fa61b6c4458c7090b1` |
| TABLEA1.DAT | All 1,081 six-field rows; inferred degree 45; derived coefficients through degree 20 | `b31b3a6f49a252a52cdb320240137dd2a881e71c022d86f421f52a07b20856b7` |

Each old/new pair passed `cmp`. The parser fixture test additionally checks independent PREM headers and raw values, duplicated-radius material interfaces, normalized radius/density, malformed-row rejection, TABLEA1 degree and selected degree-0/1 values. Its CTest is `stage3b_model_input`.

## Build and tests

Candidate build: `/tmp/gplspec-phase3-candidate-release`, Release, `BUILD_TESTING=ON`, `GPLSPEC_BUILD_BASELINE_HARNESS=ON`, `MY_PROJECT_BUILD_EXAMPLES=ON`; dependency revisions are the existing Phase 3 baseline pins. `cmake --build /tmp/gplspec-phase3-candidate-release -j2` completed all examples, tutorials, Stage 02 header/two-translation-unit checks, and stage00/stage05 reference targets. Full CTest passed 16/16 in `ctest.log`. The archived transcript has only trailing whitespace removed (5 lines); raw `/tmp` transcript SHA-256 is `6d6e5a1c224e3709f8bc70467716fea491666d44b9af666917113afaca7cdd0c`, normalized archived log SHA-256 is `9a03fbf07b2483c823e2617dd878fb740fb0c88300e2e2a9eaa3df8b4266fb9e`. The rebuilt Stage 00 runner passed three 20,700-record comparisons with zero changed components and byte-identical `MatrixSolution.out`; transcript `stage00.log` has SHA-256 `b9bdb29ec03d626c7ba6f27b6efd4bf1fd216ec4094ec1f5b7a44009fd1096c1`.

## Fixed workflow comparisons

The coordinator ran the rebuilt programs against the captured Phase 2 workflow outputs at one thread. `workflow-manifest.json` retains candidate output identities and the exact byte-mismatch list. Numeric comparison JSONs and unchanged-baseline repeats are included. The original captured-workflow evidence and retrieval recipe remain in [`../phase3a/fixed-workflows-baseline-15582ed/README.md`](../phase3a/fixed-workflows-baseline-15582ed/README.md).

- `clean_bench_2` exited 0 and both files were byte-identical.
- `clean_bench_6` exited 0. Four outputs differed bytewise and failed the frozen componentwise comparison at some near-zero components; the largest observed candidate difference was `2.384185791015625e-7`. An unchanged-baseline repeat also had strict componentwise failures, with matrix maximum `2.086162567138672e-7` and lm matrix maximum `4.172325134277344e-7`. This demonstrates baseline run variability but does not establish that every candidate delta is explained by it. Sensitivity output was byte-identical.
- `phobos_gravity` and `phobos_heterogeneous` exited 0, but their potential outputs failed the existing strict componentwise comparison. Candidate Cartesian maximum absolute difference was `1.0913936421275139e-11`; the untouched Phobos repeat's Cartesian maximum was `1.4210854715202004e-13`, while its physical output repeat also reached `1.0913936421275139e-11`. Do not claim these candidate differences are all accounted for by baseline variability. Sensitivity differences were around `1e-30` with no comparison failures.
- `OutputZiheng` exited 0. All three numeric outputs were within the existing `1e-15 + 1e-13 relative` component tolerance; `PREM200.NO.INTEGRAL.out` was byte-identical. Its two other outputs differed bytewise.

These comparisons preserve the existing tolerances and report their failures. They are not a paper-accuracy pass; the separate Stage 3A solver residual gate remains unresolved and solver code was not changed. The independent review found one follow-up for 3B.2: `OutputZiheng` parses TABLEA1 once in the `Density3D` constructor and again for output coefficients; route both uses through one owned value in that parser migration, keeping the output transform/conversion math fixed.
