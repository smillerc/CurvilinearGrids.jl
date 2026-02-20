# RFC 1111: Coding Conventions

Status: Draft

Part-Of: RFC 0000 series

Last-Updated: 2026-02-19

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

## 5. Performance and Type Stability

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

## 6. Compliance Guidance

- New/updated code SHOULD include docstrings that satisfy Section 3.
- Performance-sensitive changes SHOULD include explicit type-stability/allocation validation evidence in review notes.

