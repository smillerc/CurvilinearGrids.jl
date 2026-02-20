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
struct SphericalBasis <: BasisTrait end

#
# Independent metric caches
#

mutable struct UnifiedMetricCache{D}
  data::D
  valid::Bool
  mode::Symbol
end

mutable struct UnifiedMetricCaches{C,F}
  cell::UnifiedMetricCache{C}
  face::UnifiedMetricCache{F}
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

function _new_metric_caches(cache_mode::Symbol, cell_data, face_data)
  mode = _check_cache_mode(cache_mode)
  UnifiedMetricCaches(
    UnifiedMetricCache(cell_data, false, mode), UnifiedMetricCache(face_data, false, mode)
  )
end

#
# Unified-grid geometry/metric helpers (AoS)
#

function _allocate_unified_cell_metric_storage(::Val{N}, backend, ::Type{T}, iterators) where {N,T}
  celldims = size(iterators.cell.full)
  metric_type = _metric_eltype(Val(N), T)
  metric_array() = KernelAbstractions.zeros(backend, metric_type, celldims...)
  return (; forward=metric_array(), inverse=metric_array())
end

function _allocate_unified_face_metric_storage(::Val{N}, backend, ::Type{T}, iterators) where {N,T}
  celldims = size(iterators.cell.full)
  metric_type = _metric_eltype(Val(N), T)
  metric_array() = KernelAbstractions.zeros(backend, metric_type, celldims...)
  return ntuple(
    _ -> (; forward=metric_array(), inverse=metric_array(), conserved=metric_array()), N
  )
end

function _allocate_unified_metric_storage(::Val{N}, backend, ::Type{T}, iterators) where {N,T}
  cell = _allocate_unified_cell_metric_storage(Val(N), backend, T, iterators)
  face = _allocate_unified_face_metric_storage(Val(N), backend, T, iterators)

  return cell, face
end

function _ensure_metric_storage!(
  grid::AbstractMappedOrDiscreteGrid, ::Val{N}, ::Type{T}
) where {N,T}
  return grid.metric_caches.cell.data, grid.metric_caches.face.data
end

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

function _allocate_unified_coordinates(::Val{1}, iterators, backend, ::Type{T}) where {T}
  ni_nodes, = size(iterators.node.full)
  ni_cells, = size(iterators.cell.full)

  node_coordinates = (KernelAbstractions.zeros(backend, T, (ni_nodes,)),)
  centroid_coordinates = (KernelAbstractions.zeros(backend, T, (ni_cells,)),)
  return node_coordinates, centroid_coordinates
end

function _allocate_unified_coordinates(::Val{2}, iterators, backend, ::Type{T}) where {T}
  ni_nodes, nj_nodes = size(iterators.node.full)
  ni_cells, nj_cells = size(iterators.cell.full)

  node_coordinates = (
    KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes)),
    KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes)),
  )
  centroid_coordinates = (
    KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells)),
    KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells)),
  )
  return node_coordinates, centroid_coordinates
end

function _allocate_unified_coordinates(::Val{3}, iterators, backend, ::Type{T}) where {T}
  ni_nodes, nj_nodes, nk_nodes = size(iterators.node.full)
  ni_cells, nj_cells, nk_cells = size(iterators.cell.full)

  node_coordinates = (
    KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes, nk_nodes)),
    KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes, nk_nodes)),
    KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes, nk_nodes)),
  )
  centroid_coordinates = (
    KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells, nk_cells)),
    KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells, nk_cells)),
    KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells, nk_cells)),
  )
  return node_coordinates, centroid_coordinates
end

@inline _metric_cache_for_mapping(::Val{1}, mapping_functions, diff_backend) = MetricCache(
  mapping_functions.x, diff_backend
)
@inline _metric_cache_for_mapping(::Val{2}, mapping_functions, diff_backend) = MetricCache(
  mapping_functions.x, mapping_functions.y, diff_backend
)
@inline _metric_cache_for_mapping(::Val{3}, mapping_functions, diff_backend) = MetricCache(
  mapping_functions.x, mapping_functions.y, mapping_functions.z, diff_backend
)

function _build_unified_components(
  ::Val{N},
  mapping_functions,
  celldims::NTuple{N,Int},
  discretization_scheme::Symbol,
  backend,
  diff_backend,
  ::Type{T};
  global_cell_indices=nothing,
) where {N,T}
  GradientDiscretizationScheme, order, _, nhalo, scheme_name = get_gradient_discretization_scheme(
    discretization_scheme
  )

  iterators = get_iterators(celldims, nhalo, global_cell_indices)
  node_coordinates, centroid_coordinates = _allocate_unified_coordinates(
    Val(N), iterators, backend, T
  )
  metric_functions_cache = _metric_cache_for_mapping(Val(N), mapping_functions, diff_backend)
  cell_metric_storage, face_metric_storage = _allocate_unified_metric_storage(
    Val(N), backend, T, iterators
  )
  discr_scheme = GradientDiscretizationScheme(order; use_cache=false)

  return (;
    node_coordinates,
    centroid_coordinates,
    metric_functions_cache,
    cell_metric_storage,
    face_metric_storage,
    nhalo,
    discretization_scheme=discr_scheme,
    discretization_scheme_name=scheme_name,
    iterators,
  )
end

function _compute_unified_node_coordinates!(
  node_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{1},
)
  x = mapping_functions.x
  @threads for I in iterators.node.full
    Iglobal = iterators.global_domain.node.full[I]
    ξ = Iglobal.I .- nhalo
    node_coordinates[1][I] = x(t, ξ..., params)
  end
  return nothing
end

function _compute_unified_node_coordinates!(
  node_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{2},
)
  x = mapping_functions.x
  y = mapping_functions.y
  @threads for I in iterators.node.full
    Iglobal = iterators.global_domain.node.full[I]
    ξη = Iglobal.I .- nhalo
    node_coordinates[1][I] = x(t, ξη..., params)
    node_coordinates[2][I] = y(t, ξη..., params)
  end
  return nothing
end

function _compute_unified_node_coordinates!(
  node_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{3},
)
  x = mapping_functions.x
  y = mapping_functions.y
  z = mapping_functions.z
  @threads for I in iterators.node.full
    Iglobal = iterators.global_domain.node.full[I]
    ξηζ = Iglobal.I .- nhalo
    node_coordinates[1][I] = x(t, ξηζ..., params)
    node_coordinates[2][I] = y(t, ξηζ..., params)
    node_coordinates[3][I] = z(t, ξηζ..., params)
  end
  return nothing
end

function _compute_unified_centroid_coordinates!(
  centroid_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{1},
)
  x = mapping_functions.x
  @threads for I in iterators.cell.full
    Iglobal = iterators.global_domain.cell.full[I]
    ξ = Iglobal.I .- nhalo .+ 0.5
    centroid_coordinates[1][I] = x(t, ξ..., params)
  end
  return nothing
end

function _compute_unified_centroid_coordinates!(
  centroid_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{2},
)
  x = mapping_functions.x
  y = mapping_functions.y
  @threads for I in iterators.cell.full
    Iglobal = iterators.global_domain.cell.full[I]
    ξη = Iglobal.I .- nhalo .+ 0.5
    centroid_coordinates[1][I] = x(t, ξη..., params)
    centroid_coordinates[2][I] = y(t, ξη..., params)
  end
  return nothing
end

function _compute_unified_centroid_coordinates!(
  centroid_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{3},
)
  x = mapping_functions.x
  y = mapping_functions.y
  z = mapping_functions.z
  @threads for I in iterators.cell.full
    Iglobal = iterators.global_domain.cell.full[I]
    ξηζ = Iglobal.I .- nhalo .+ 0.5
    centroid_coordinates[1][I] = x(t, ξηζ..., params)
    centroid_coordinates[2][I] = y(t, ξηζ..., params)
    centroid_coordinates[3][I] = z(t, ξηζ..., params)
  end
  return nothing
end

@inline function _inverse_and_normalized_edge_metrics(
  ::Val{1}, edge, axis::Int, t, ξηζ, params
)
  if axis != 1
    throw(ArgumentError("Invalid 1D face axis: $axis"))
  end

  jinv_edge = edge.Jinv_ᵢ₊½
  norm_edge = edge.norm_Jinv_ᵢ₊½

  jinv_fun = jinv_edge isa NamedTuple ? jinv_edge.ϕᵢ₊½ : jinv_edge
  norm_fun = norm_edge isa NamedTuple ? norm_edge.ϕᵢ₊½ : norm_edge

  G = _as_smatrix(Val(1), jinv_fun(t, ξηζ..., params))
  Ghat = _as_smatrix(Val(1), norm_fun(t, ξηζ..., params))
  return G, Ghat
end

@inline function _inverse_and_normalized_edge_metrics(
  ::Val{2}, edge, axis::Int, t, ξηζ, params
)
  if axis == 1
    G = _as_smatrix(Val(2), edge.Jinv_ᵢ₊½(t, ξηζ..., params))
    Ghat = _as_smatrix(Val(2), edge.norm_Jinv_ᵢ₊½(t, ξηζ..., params))
  elseif axis == 2
    G = _as_smatrix(Val(2), edge.Jinv_ⱼ₊½(t, ξηζ..., params))
    Ghat = _as_smatrix(Val(2), edge.norm_Jinv_ⱼ₊½(t, ξηζ..., params))
  else
    throw(ArgumentError("Invalid 2D face axis: $axis"))
  end
  return G, Ghat
end

@inline function _inverse_and_normalized_edge_metrics(
  ::Val{3}, edge, axis::Int, t, ξηζ, params
)
  if axis == 1
    G = _as_smatrix(Val(3), edge.Jinvᵢ₊½(t, ξηζ..., params))
    Ghat = @SMatrix [
      edge.ξ̂xᵢ₊½(t, ξηζ..., params) edge.ξ̂yᵢ₊½(t, ξηζ..., params) edge.ξ̂zᵢ₊½(t, ξηζ..., params)
      edge.η̂xᵢ₊½(t, ξηζ..., params) edge.η̂yᵢ₊½(t, ξηζ..., params) edge.η̂zᵢ₊½(t, ξηζ..., params)
      edge.ζ̂xᵢ₊½(t, ξηζ..., params) edge.ζ̂yᵢ₊½(t, ξηζ..., params) edge.ζ̂zᵢ₊½(t, ξηζ..., params)
    ]
  elseif axis == 2
    G = _as_smatrix(Val(3), edge.Jinvⱼ₊½(t, ξηζ..., params))
    Ghat = @SMatrix [
      edge.ξ̂xⱼ₊½(t, ξηζ..., params) edge.ξ̂yⱼ₊½(t, ξηζ..., params) edge.ξ̂zⱼ₊½(t, ξηζ..., params)
      edge.η̂xⱼ₊½(t, ξηζ..., params) edge.η̂yⱼ₊½(t, ξηζ..., params) edge.η̂zⱼ₊½(t, ξηζ..., params)
      edge.ζ̂xⱼ₊½(t, ξηζ..., params) edge.ζ̂yⱼ₊½(t, ξηζ..., params) edge.ζ̂zⱼ₊½(t, ξηζ..., params)
    ]
  elseif axis == 3
    G = _as_smatrix(Val(3), edge.Jinvₖ₊½(t, ξηζ..., params))
    Ghat = @SMatrix [
      edge.ξ̂xₖ₊½(t, ξηζ..., params) edge.ξ̂yₖ₊½(t, ξηζ..., params) edge.ξ̂zₖ₊½(t, ξηζ..., params)
      edge.η̂xₖ₊½(t, ξηζ..., params) edge.η̂yₖ₊½(t, ξηζ..., params) edge.η̂zₖ₊½(t, ξηζ..., params)
      edge.ζ̂xₖ₊½(t, ξηζ..., params) edge.ζ̂yₖ₊½(t, ξηζ..., params) edge.ζ̂zₖ₊½(t, ξηζ..., params)
    ]
  else
    throw(ArgumentError("Invalid 3D face axis: $axis"))
  end

  return G, Ghat
end

function _fill_cell_metric_storage!(
  cell_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
) where {N,T}
  @threads for I in iterators.cell.full
    Iglobal = iterators.global_domain.cell.full[I]
    ξηζ = Iglobal.I .- nhalo .+ 0.5

    F = _as_smatrix(Val(N), metric_functions_cache.forward.jacobian(t, ξηζ..., params))
    G = _as_smatrix(Val(N), metric_functions_cache.inverse.Jinv(t, ξηζ..., params))
    J = det(F)

    cell_metric_storage.forward[I] = _metric_from_jacobian(
      SMatrix{N,N,T,N * N}(Tuple(F)), T(J)
    )
    cell_metric_storage.inverse[I] = _metric_from_jacobian(
      SMatrix{N,N,T,N * N}(Tuple(G)), T(J)
    )
  end
  return nothing
end

function _fill_face_metric_storage!(
  face_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
) where {N,T}
  edge = metric_functions_cache.edge

  @threads for I in iterators.cell.full
    Iglobal = iterators.global_domain.cell.full[I]
    ξηζ = Iglobal.I .- nhalo .+ 0.5

    for axis in 1:N
      G, Ghat = _inverse_and_normalized_edge_metrics(Val(N), edge, axis, t, ξηζ, params)
      F = inv(G)
      J = det(F)

      face_metric_storage[axis].forward[I] = _metric_from_jacobian(
        SMatrix{N,N,T,N * N}(Tuple(F)), T(J)
      )
      face_metric_storage[axis].inverse[I] = _metric_from_jacobian(
        SMatrix{N,N,T,N * N}(Tuple(G)), T(J)
      )
      face_metric_storage[axis].conserved[I] = _metric_from_jacobian(
        SMatrix{N,N,T,N * N}(Tuple(Ghat)), T(J)
      )
    end
  end
  return nothing
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
_basis_trait_from_legacy(::AbstractCurvilinearGrid) = CartesianBasis()
