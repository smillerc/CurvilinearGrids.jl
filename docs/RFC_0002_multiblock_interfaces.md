# RFC 0002: Multi-Block Interface System for Unified Grids

- Status: Draft
- Authors: CurvilinearGrids maintainers
- Depends on: `docs/RFC_0001_grid_unification.md`

## 1. Summary

This RFC adds a multi-block mesh module that connects multiple unified grids
(`MappedGrid` and `DiscreteGrid`) through explicit interfaces.

Core constraints:
- Interfaces are strictly 1-to-1.
- Interface faces must be geometrically abutting.
- Only mesh-mesh interfaces are tracked; non-interface faces are implicit outer boundaries.

Coordinate systems and basis traits may differ across an interface, so transfer
operators must support coordinate/basis transformation.

## 2. Motivation

Large geometries are commonly decomposed into structured blocks. The current
single-grid API cannot represent block topology or interface transfer rules.

The goal is to provide:
- explicit block connectivity,
- robust geometric validation,
- reusable interface transfer maps,
- a foundation for conservative multi-block operators.

## 3. Goals and Non-Goals

Goals:
- Represent multi-block topology over unified grids.
- Enforce strict 1-to-1 abutting interfaces.
- Support scalar/vector/tensor transfer across coordinate-system changes.
- Cache interface mappings for repeated exchange.

Non-goals:
- Arbitrary overset/chimera connectivity.
- One-to-many interface couplings.
- AMR/coarsened non-conformal joins in this RFC.

## 4. Scope

In scope:
- `MappedGrid` and `DiscreteGrid` only.
- Interface topology, validation, and transfer.
- Boundary face classification.

Out of scope:
- `OrthogonalGrid` participation in multi-block connectivity (initially).
- Flux correction/reconstruction schemes beyond baseline interpolation transfer.

## 5. Proposed Module Layout

Add:
- `src/grids/multiblock/multiblock.jl`
- `src/grids/multiblock/types.jl`
- `src/grids/multiblock/validation.jl`
- `src/grids/multiblock/transfer.jl`
- `src/grids/multiblock/cache.jl`

Exports (initial):
- `MultiBlockMesh`
- `BlockFace`
- `BlockInterface`
- `validate_multiblock!`
- `build_interface_caches!`
- `exchange_interface!`
- `exchange_all_interfaces!`

## 6. Data Model

### 6.1 Face Identification

```julia
struct BlockFace{N}
  block_id::Int
  axis::Int          # 1..N
  side::Symbol       # :min or :max
end
```

### 6.2 Interface Definition

```julia
struct BlockInterface{N,T}
  left::BlockFace{N}
  right::BlockFace{N}
  permutation::NTuple{N-1,Int}
  flips::NTuple{N-1,Bool}
  tolerance::T
end
```

`permutation` and `flips` describe index-space alignment between interface
planes.

### 6.3 Multi-Block Container

```julia
struct MultiBlockMesh{N,T,B,I,C}
  blocks::B
  interfaces::I
  cache::C
end
```

`blocks` must contain only `MappedGrid{N,T,...}` or `DiscreteGrid{N,T,...}`.

## 7. Topology and Geometry Invariants

The following are required:

1. Interface uniqueness:
- A face may appear in at most one interface endpoint.
- Non-interface faces are treated as outer boundaries and are not tracked.

2. 1-to-1 connectivity:
- A face can belong to at most one `BlockInterface`.
- Each `BlockInterface` binds exactly two faces.

3. Abutting geometry:
- Interface face pairs must coincide geometrically within tolerance.
- Face normals must be opposite at validation samples.

4. Dimensional consistency:
- Both sides of an interface must have matching `N`.
- Interface plane dimensions must agree in node/cell counts after declared
  permutation/flips.

Validation failure must throw with actionable diagnostics listing offending
faces and violation type.

## 8. Transfer Across Coordinate Systems

Interfaces may connect blocks with different:
- coordinate systems (`coordinate_system(typeof(grid))`),
- basis traits (`basis_trait(typeof(grid))`).

Transfer API must dispatch by field rank:
- `transfer_scalar`
- `transfer_vector`
- `transfer_tensor`

Behavior:
- Scalar transfer: interpolation-only.
- Vector/tensor transfer: interpolation + component transformation into
  receiver-side representation.

## 9. Interface Cache

Each interface cache entry stores:
- receiver interface points (or indices),
- mapped donor computational coordinates,
- interpolation stencil/weights,
- transform metadata (if vector/tensor field).

Cache is invalidated when either connected block is updated.

## 10. Public API (Target)

```julia
validate_multiblock!(mb::MultiBlockMesh)
build_interface_caches!(mb::MultiBlockMesh)

exchange_interface!(mb::MultiBlockMesh, iface_id, field; field_kind=:scalar)
exchange_all_interfaces!(mb::MultiBlockMesh, field; field_kind=:scalar)
```

Preferred declaration syntax:

```julia
interfaces = ((mesh_1, :ilo) => (mesh_2, :jhi),)
mb = MultiBlockMesh((mesh_1, mesh_2), interfaces)
```

Computational-coordinate transfer across one interface:

```julia
ξ₂ = computational_coordinate(mb, iface_id, ξ₁; from=:left)
```

Field-kind options initially:
- `:scalar`
- `:vector`
- `:tensor`

## 11. Updated Implementation Sequence

Phase 1: Topology Core (strict by default)
- Add module/types.
- Implement interface-face uniqueness checks.
- Enforce strict 1-to-1 face pairing.
- Treat non-interface faces as implicit outer boundaries.

Phase 2: Abutting Geometry Validation
- Add AABB precheck for candidate rejection.
- Add sampled geometric coincidence checks.
- Add opposite-normal checks.
- Add detailed error reporting.

Phase 3: Scalar Interface Transfer
- Implement donor-to-receiver interpolation for scalar fields.
- Add per-interface cache for mapped donor coordinates/stencils.
- Invalidate/rebuild cache on block updates.

Phase 4: Coordinate/Basis Transformation
- Add vector/tensor transfer dispatch.
- Implement trait-based component transformation between blocks.
- Add tests for mixed coordinate-system interfaces.

Phase 5: Performance and Robustness
- Benchmark cache build + repeated exchange.
- Add tolerance controls and diagnostics.
- Add optional strictness flags only if needed; defaults remain strict.

## 12. Minimal Test Matrix

1. Topology:
- Duplicate face assignment across interfaces fails.

2. Interface rules:
- Non 1-to-1 interface declaration fails.
- Non-abutting interface pair fails.

3. Transfer:
- Constant scalar field preserved across interface.
- Vector field transfer passes coordinate-system transform sanity checks.

4. Cache:
- Reuse cache on repeated exchange.
- Invalidate cache after `update!` on either block.

## 13. Acceptance Criteria

This RFC is accepted when:
- Multi-block module compiles and is exported.
- Strict 1-to-1 and abutting validation is enforced.
- Scalar transfer works for `MappedGrid`/`DiscreteGrid` interfaces.
- Coordinate/basis transfer hooks exist for vector/tensor fields.
- Test matrix in Section 12 passes.
