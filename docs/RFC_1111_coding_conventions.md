# RFC 1111: Coding Conventions

Status: Draft

Part-Of: RFC 0000 series

Last-Updated: 2026-02-28

## 1. Scope

This RFC defines coding conventions for CurvilinearGrids source code and documentation.

## 2. Normative Language

The key words MUST, SHOULD, and MAY indicate requirement strength.

## 3. Documentation Conventions

### 3.1 Function Docstrings

Function docstrings MUST follow this template (omit empty sections):

```julia
"""
    function_signature(arg1, arg2; kwarg=value)

Brief description of the function's purpose.

More detailed explanation if necessary, covering behavior, edge cases, etc.

# Arguments
  - `arg1`: Description of the first argument.
  - `arg2`: Description of the second argument.

# Keywords
  - `kwarg`: Description of the keyword argument and its default value.

# Returns
Description of what the function returns.
"""
```

### 3.2 Struct Docstrings

Struct docstrings MUST follow this template:

```julia
"""

Brief description of the struct's purpose.

More detailed explanation if necessary, covering behavior, edge cases, etc.

# Fields
  - `field1`: Description of the first data member.
  - `field2`: Description of the second data member.
"""
```

### 3.3 Module Docstrings

Module docstrings SHOULD provide a concise high-level description.

## 4. Formatting

- Source formatting MUST use JuliaFormatter with repository settings from `.JuliaFormatter.toml`.

## 5. Module Imports

1. Module imports MUST be explicit and symbol-scoped.
2. Code MUST prefer `using PackageName: symbol1, symbol2` over broad `using PackageName`.
3. `import PackageName: symbol` MUST be used when extending methods from another module.
4. Wildcard-style imports and implicit namespace pollution SHOULD be avoided in all production source files.
5. Example (required style):

```julia
using CurvilinearGrids: centroids, cellvolume
import CurvilinearGrids: face_metric_coefficient
```

## 6. Performance and Type Stability

1. Keep concrete field and container types.
2. Avoid `Vector{Any}`, `Dict{Symbol,Any}`, and untyped struct fields in hot paths.
3. Prefer parametric structs over abstractly typed fields.
4. Preserve type stability across branches.
5. Avoid type-changing local variables in performance-critical code.
6. Use function barriers to separate setup/parsing from tight kernels.
7. Preallocate reusable buffers in stencil/block loops.
8. Use views or explicit loops to avoid unnecessary temporaries.
9. Use `@inbounds` only when safety is guaranteed by invariants/tests.
10. Keep global state `const` and avoid mutable globals in kernels.
11. Validate critical kernels with `@code_warntype`, `@allocated`, and benchmarks.

## 7. HPC-Oriented Coding Style

### 7.1 Core Principle

HPC kernels MUST prioritize correctness and throughput while remaining legible to domain developers.
Code SHOULD be simple first, but not at the cost of avoidable allocations, type instability, poor memory locality, or branch-heavy inner loops.

### 7.2 Mathematical Notation in Code

1. Loop variables SHOULD follow mathematical conventions when they map directly to discretization indices, e.g. `i, j, k`, `d`.
2. Intermediate quantities in stencil/operator kernels SHOULD use mathematically meaningful names, e.g. `Δx`, `∂u∂x`, `J`, `gdd`, `divergence_x`.
3. Formula implementations SHOULD preserve visible correspondence to governing equations (e.g. flux form vs divergence form) so reviewers can verify physics and discretization quickly.
4. When Unicode symbols are already used in a file, mathematically standard symbols MAY be used (`Δ`, `θ`, `ϕ`) if they improve clarity and do not reduce maintainability.

### 7.3 Loop and Kernel Structure

1. Hot loops MUST keep indexing and arithmetic explicit; avoid hidden allocations in comprehensions/broadcast inside kernel bodies.
2. Inner loops SHOULD be contiguous-stride where possible; loop ordering SHOULD favor cache locality and SIMD vectorization.
3. Kernel code SHOULD separate:
   - geometry/metric fetch,
   - local stencil arithmetic,
   - writeback.
4. Repeated terms SHOULD be hoisted into locals when this reduces recomputation without obscuring the formula.
5. Parallel decomposition (thread/task tiling) SHOULD be coarse enough to amortize scheduling overhead and keep writes disjoint.

### 7.4 Naming and Readability

1. Struct fields MUST use descriptive names that describe physical or algorithmic meaning (`centroid_coordinates`, `face_areas`, `cell_volumes`) instead of opaque abbreviations.
2. In functions with long type-qualified names, local aliases MAY be introduced for readability if they are unambiguous and local in scope.
3. Public API names SHOULD be descriptive and stable; helper names in hot paths MAY be shorter if they remain interpretable in context.
4. Avoid shim indirections that hide direct API usage when the underlying package already provides the required operation.

### 7.5 Simplicity vs Performance

1. Prefer the simplest implementation that meets performance requirements.
2. If a less-obvious optimization is required, include a short comment documenting:
   - what is optimized,
   - why it matters (allocation, cache, SIMD, synchronization),
   - invariants relied upon (`@inbounds`, halo width, shape assumptions).
3. Do not accept readability-only refactors that regress measured performance in hot kernels.
4. Do not accept micro-optimizations that materially reduce clarity unless benchmarks show meaningful benefit.

## 8. Compliance Guidance

- New/updated code SHOULD include docstrings that satisfy Section 3.
- Performance-sensitive changes SHOULD include explicit type-stability/allocation validation evidence in review notes.
