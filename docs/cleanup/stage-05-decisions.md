# Stage 05 operator decisions

Stage 05 began on `cleanup/05-operators` from accepted integrated stage-04 base `e352efcf09bcf0d0fbbbded62f1b218736e2fe79`. The authorized boundary is behavior-preserving extraction of equivalent operator kernels. Preserve scaled gradients, radial weights, canonical signs, derivative slots, scatter-add order, exterior boundary terms, public layouts, and both wrapper interfaces.

| Requested substep | Initial inventory | Boundary / validation |
|---|---|---|
| Scalar-gradient coefficient construction | The 1D and 3D wrappers have distinct radial scaling and coefficient multipliers. | Do not unify whole gradient construction. A smaller coefficient operation is eligible only if exact old-expression comparisons prove equivalence. |
| Canonical tensor-vector contraction | Both wrappers visibly use the same nine tensor slots and `-,+,-` signs to form three components. | Candidate for a small internal kernel with wrapper-specific access adapters. Test deterministic complex vectors/tensors against copied original statements and exercise both wrappers. |
| `MappingPerturbation` vector-gradient construction | The map-based and file-coefficient constructors appear to assemble the same 3x3 df tensor from nine spatial derivative vectors. | Audit all source expressions, loop bounds, and component ordering before deciding; compare each entry against copied original assignments. |
| Perturbed Laplace tensor | Both nontrivial `MappingPerturbation` constructors visibly apply the same metric, inverse-F, Laplace tensor, scalar trace, and two transpose-ordered subtractions. | Preserve matrix operation and subtraction order exactly; compare complex matrices against the copied original sequence. |
| Radial weak assembly | Wrapper radial multipliers and scale factors differ; both also have caller-specific boundaries and scatter behavior. | Do not merge assembly absent a smaller exact common operation. Direct original-source operator comparisons must include both wrapper paths, source vectors, perturbation operators, and boundary characterization. |

## Coverage gap and evidence plan

The frozen stage-00 reference harness covers five fixtures and records full model/perturbation fields, a deterministic complex 3D operator product, source, solution, and representative output; it constructs only `MatrixReplacement3D`. Its committed CSV/output are immutable. Added `stage05_operator_reference` for the 1D path, compiled identically in a detached original-source worktree at `4ef3a66c62d64408d99989dd51c3ccbdc46d0b0c` and in the candidate. `tests/reference/stage05-1d-original.csv` is captured only from that original-source executable, SHA-256 `a9e89b9ac1250f3f01b695ec5be0fcb79990a85da42e6eb80ae001f5f88a43d6`; `stage05_1d_original_reference` checks the candidate against it. It records all 387 complex operator entries, the full force vector, explicit origin and 13 material-interface action values, boundary-on/off action delta, and two complex sesquilinear forms. The original 1D sample is not exactly Hermitian; the measured values are retained as evidence, with no repair. Boundary delta is nonzero on exactly the nine exterior harmonic entries.

Direct 3D original-source and candidate harness CSVs are byte-identical at SHA-256 `ef20bd7e3016661c60903c290c74c599151fdfc0417443bcf81c4215401f105e`, covering operator action, source, perturbation tensors, and solutions; representative output files compare byte-for-byte. The frozen stage00 regression remains required and unchanged. Original build provenance and current substep inventory remain documented in `baseline-manifest.md` and this file.

The full perturbation vector-gradient and perturbed-Laplace source audit is ongoing. The unused `Mapping_Tools::dxitodf` path is not treated as equivalent to the active `MappingPerturbation` constructor paths.

## Completed extraction 1 — canonical tensor-vector contraction

The 1D and 3D wrappers now call the small internal contraction helper. The copied-original-expression test is bit-exact for 1,000 deterministic dense complex cases. Candidate tests `stage05_tensor_contraction` and `stage05_1d_original_reference`, the stage02 patch guard, and the immutable three-comparison stage00 regression pass. Direct candidate-versus-original-source (4ef3a66) stage00 CSVs are byte-identical and their representative output files match; the 1D output fixture is frozen only from the original executable. This extraction did not alter scalar-gradient multipliers, radial quadrature weights, weak assembly, scatter-add, or boundary formulas.
