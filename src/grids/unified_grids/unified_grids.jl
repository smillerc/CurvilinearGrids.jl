# Unified grid shared definitions and includes.

"""
Unified grid model (Phase 1, additive):
- MappedGrid
- DiscreteGrid
- OrthogonalGrid
"""

"""
Base abstract type for all unified-grid implementations.
"""
abstract type AbstractUnifiedGrid end

"""
Base abstract type for unified grids that own mapping/metric caches.
"""
abstract type AbstractMappedOrDiscreteGrid <: AbstractUnifiedGrid end

#
# Independent metric caches
#

"""
    UnifiedMetricCache{D}

Cache wrapper for one metric domain (cell or face).

# Fields
  - `data`: Backing metric storage container.
  - `valid`: Cache validity flag.
  - `mode`: Cache mode (`:eager`, `:lazy`, `:off`).
"""
mutable struct UnifiedMetricCache{D}
  data::D
  valid::Bool
  mode::Symbol
end

"""
    UnifiedMetricCaches{C,F}

Container for independently managed cell and face metric caches.

# Fields
  - `cell`: Cell-center metric cache.
  - `face`: Face metric cache.
"""
mutable struct UnifiedMetricCaches{C,F}
  cell::UnifiedMetricCache{C}
  face::UnifiedMetricCache{F}
end

function _check_cache_mode(cache_mode::Symbol)
  if cache_mode ∉ (:eager, :lazy, :off)
    throw(
      ArgumentError(
        "Invalid cache mode `$cache_mode`. Expected one of `:eager`, `:lazy`, `:off`."
      ),
    )
  end
  return cache_mode
end

function _new_metric_caches(cache_mode::Symbol, cell_data, face_data)
  mode = _check_cache_mode(cache_mode)
  UnifiedMetricCaches(
    UnifiedMetricCache(cell_data, false, mode), UnifiedMetricCache(face_data, false, mode)
  )
end

@inline function _has_metric_storage(grid::AbstractMappedOrDiscreteGrid)
  return grid.metric_caches !== nothing
end

@inline function _require_metric_storage(
  grid::AbstractMappedOrDiscreteGrid, caller::AbstractString
)
  if !_has_metric_storage(grid)
    throw(
      ArgumentError(
        "`$caller` is unavailable because this grid was created without metric storage (`compute_metrics=false, cache_mode=:off`).",
      ),
    )
  end
  return nothing
end

@inline function _has_metric_functions(grid::AbstractMappedOrDiscreteGrid)
  return grid.metric_functions_cache !== nothing
end

@inline function _require_metric_functions(
  grid::AbstractMappedOrDiscreteGrid, caller::AbstractString
)
  if !_has_metric_functions(grid)
    throw(
      ArgumentError(
        "`$caller` is unavailable because this grid has no metric function cache."
      ),
    )
  end
  return nothing
end

#


# Constructor validation helpers.

@inline function _check_unified_basis_trait(basis::BasisTrait)
  if basis isa CartesianBasis || basis isa SphericalBasis
    return basis
  end
  throw(
    ArgumentError(
      "Unsupported basis trait $(typeof(basis)) for unified grids. Use `CartesianBasis()` or `SphericalBasis()`.",
    ),
  )
end

include("metrics.jl")
include("mapped_grid.jl")
include("discrete_grid.jl")
include("generation.jl")
include("api.jl")
