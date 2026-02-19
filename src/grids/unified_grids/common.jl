"""
Unified grid model (Phase 1, additive):
- MappedGrid
- DiscreteGrid
- OrthogonalGrid
"""

abstract type AbstractUnifiedGrid end
abstract type AbstractMappedOrDiscreteGrid <: AbstractUnifiedGrid end

#
# Traits
#

abstract type CoordinateSystemTrait end
struct CartesianCS <: CoordinateSystemTrait end
struct CylindricalCS <: CoordinateSystemTrait end
struct SphericalCS <: CoordinateSystemTrait end
struct AxisymmetricCS{Axis} <: CoordinateSystemTrait end
struct CurvilinearCS <: CoordinateSystemTrait end

abstract type BasisTrait end
struct CartesianBasis <: BasisTrait end
struct ContravariantBasis <: BasisTrait end
struct CovariantBasis <: BasisTrait end
struct SphericalBasis <: BasisTrait end

#
# Independent metric caches
#

mutable struct UnifiedMetricCache
  data::Any
  valid::Bool
  mode::Symbol
end

mutable struct UnifiedMetricCaches
  cell::UnifiedMetricCache
  face::UnifiedMetricCache
end

function _check_cache_mode(cache_mode::Symbol)
  if cache_mode ∉ (:eager, :lazy, :off)
    throw(
      ArgumentError(
        "Invalid cache mode `$cache_mode`. Expected one of `:eager`, `:lazy`, `:off`.",
      ),
    )
  end
  return cache_mode
end

function _new_metric_caches(cache_mode::Symbol)
  mode = _check_cache_mode(cache_mode)
  UnifiedMetricCaches(
    UnifiedMetricCache(nothing, false, mode), UnifiedMetricCache(nothing, false, mode)
  )
end

#
# Legacy trait inference
#

_coordinate_system_from_legacy(::CartesianOrthogonalGrid1D) = CartesianCS()
_coordinate_system_from_legacy(::CylindricalOrthogonalGrid1D) = CylindricalCS()
_coordinate_system_from_legacy(::SphericalOrthogonalGrid1D) = SphericalCS()
_coordinate_system_from_legacy(::SphericalGrid3D) = SphericalCS()
_coordinate_system_from_legacy(::SphericalGrid1D) = SphericalCS()
_coordinate_system_from_legacy(::CylindricalGrid1D) = CylindricalCS()
_coordinate_system_from_legacy(::SphericalBasisCurvilinearGrid3D) = SphericalCS()

function _coordinate_system_from_legacy(mesh::AxisymmetricGrid2D)
  if mesh.rotational_axis === :x
    return AxisymmetricCS{:x}()
  else
    return AxisymmetricCS{:y}()
  end
end

_coordinate_system_from_legacy(::AxisymmetricOrthogonalGrid2D) = AxisymmetricCS{:y}()
_coordinate_system_from_legacy(::AbstractCurvilinearGrid) = CurvilinearCS()

_basis_trait_from_legacy(::SphericalBasisCurvilinearGrid3D) = SphericalBasis()
_basis_trait_from_legacy(::AbstractCurvilinearGrid) = ContravariantBasis()
