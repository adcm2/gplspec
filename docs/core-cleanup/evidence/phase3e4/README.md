# Phase 3E.4 sensitivity evidence

Candidate changes consolidate the single and vector sensitivity selector assembly in `Gravity_Tools::detail::SensitivitySelectorRhs`. The single selector still uses ConjugateGradient; the vector selector still uses BiCGSTAB. Both retain the preconditioned guess, caller tolerance, historical decimal pi in the boundary RHS, and nested output conversion. Vector duplicates retain last-value-wins. Length, harmonic range, finite weights, and positive finite tolerance checks now reject invalid inputs before solve/indexing.

The focused test is `stage3e_sensitivity`. Reproduce it with the existing candidate build:

```sh
cmake -S . -B /tmp/gplspec-phase3-resume/candidate-build
cmake --build /tmp/gplspec-phase3-resume/candidate-build --target stage3e_sensitivity -j1
ctest --test-dir /tmp/gplspec-phase3-resume/candidate-build -R '^stage3e_sensitivity$' --output-on-failure
cmake --build /tmp/gplspec-phase3-resume/candidate-build --target clean_bench_6 phobos_gravity phobos_heterogeneous -j1
```

The independent fixture uses a homogeneous body of normalized radius `0.63`, outer radius `1.2`, `p=5`, `lmax=2`, maximum step `0.1`, Earth length/time/mass norms, and a degree-two axisymmetric fractional density change. The constructor's `1+0.005*GetValueAt` behavior is explicitly inverted in the test tomography so the full models have `delta rho = epsilon*rho*(r/R)^2*Y20`. It constructs both positive and negative full density models at `epsilon=1e-4` and `5e-5` and solves them at test-only tolerance `1e-13`.

The physical derivative prediction at the outer node is

`-4*pi*G_phys*rho_phys*R_phys^5/(35*Rb_phys^3)`.

An independent radial quadrature contracts the returned sensitivity kernel with the known harmonic density coefficient using element-local positive weights `q.W(j)*(h_e/2)*r_j^2`. The kernel already carries the historical `4*pi_legacy*G` boundary RHS normalization; the source normalization cancels in the adjoint contraction. Since `pi_legacy=3.1415926535`, its expected kernel contraction is the full-pi analytic value times `pi_legacy/pi_full`. The 2.86e-11 relative normalization difference is reported separately and is not hidden by the `1e-4` finite-difference criterion. Kernel values are sensitivity coefficients, not potential values; only the contracted derivative is reported in `m^2/s^2` after `PotentialNorm` conversion.

Release CTest passed 1/1. The analytic physical derivative and independent kernel contraction are both `-2.58017e6 m^2/s^2`. The full-model central differences are `-2.58017e6 m^2/s^2` at both steps, with relative errors `3.2071e-12` and `1.28542e-12`, below the predeclared `1e-4` bound. Zero-amplitude source equality and plus/minus symmetry pass the near-zero floor `1e-8*|Phi_base,20| + 1e-12*G_phys*M/Rb`. Selector tests verify the `l*l+l-m` slot and `(-1)^m` phase for `m=+/-1`, empty selectors, length/range/finite/tolerance rejection, and duplicate last-value-wins. A weighted duplicate vector response matches twice the single-selector kernel within the predeclared `1e-8` bound, exercising the preserved BiCGSTAB route.

The supported caller targets `clean_bench_6`, `phobos_gravity`, and `phobos_heterogeneous` compiled and linked. They were not run in this checkpoint; full supported-workflow execution and final comparisons remain in 3F. `experimental/` and frozen references were unchanged. No paper-figure claim follows from this fixture.
