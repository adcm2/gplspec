# Stage 04 model-construction decisions

| Category | Current assessment | Semantic boundary / validation
|---|---|---|
| Mesh and quadrature setup | In inventory. Constructors commonly create GLL quadrature and a transform grid, but mesh calls differ in physical/referential radii and constructor inputs. | No shared setup extraction until operation order and `RadialMesh` semantics are confirmed identical.
| Storage initialization | First bounded extraction candidate: repeated nested scalar-field allocation for mapping and Jacobian arrays, with constructor-specific zero/one initial values. | Preserve shape `(elements, quadrature nodes, spatial samples)` and each original initial value. Validate exact nested arrays and frozen model intermediates.
| Density sampling | In inventory. Physical constructors multiply their physical density by the Jacobian; model/tomography constructors sample radial density and may apply tomography without the same factor. | Keep physical versus referential semantics explicit; no shared sampler across these paths without exact equivalence.
| Mapping geometry | In inventory. Radial-mapping constructors clamp radius and apply an exterior taper; the aspherical file constructor builds geometry from transformed radius coefficients. | Distinct geometry construction remains separate unless a smaller operation is proven identical. Compare mapping, Jacobian, inverse-F and Laplace-tensor arrays.

Each extraction is a separate logical commit. Original-source stage-00 references remain fixed; expected outputs are never regenerated from the candidate.
