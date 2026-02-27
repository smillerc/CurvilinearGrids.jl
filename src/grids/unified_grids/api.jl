#
# Cache API
#

#
# Trait helpers
#

"""
    coordinate_system(grid::AbstractUnifiedGrid)

Return the coordinate-system trait associated with a unified grid.

# Arguments
  - `grid`: Unified grid instance.

# Returns
Coordinate-system trait instance (for example `CurvilinearCS()` or `SphericalCS()`).
"""
coordinate_system(grid::AbstractUnifiedGrid) = coordinate_system(typeof(grid))
coordinate_system(grid::AbstractOrthogonalGrid) = coordinate_system(typeof(grid))

"""
    basis_trait(grid::Union{MappedGrid,DiscreteGrid})

Return the basis trait associated with a mapped or discrete unified grid.

# Arguments
  - `grid`: Mapped or discrete unified grid instance.

# Returns
Basis trait instance (`CartesianBasis()` or `SphericalBasis()`).
"""
basis_trait(grid::Union{MappedGrid,DiscreteGrid}) = basis_trait(typeof(grid))
function basis_trait(::AbstractOrthogonalGrid)
  throw(ArgumentError("`basis_trait` is undefined for `AbstractOrthogonalGrid`."))
end

coordinate_system(::Type{<:MappedGrid{N,T,CS}}) where {N,T,CS} = CS()
coordinate_system(::Type{<:DiscreteGrid{N,T,CS}}) where {N,T,CS} = CS()
coordinate_system(::Type{<:OrthogonalGrid{N,T,CS}}) where {N,T,CS} = CS()

basis_trait(::Type{<:MappedGrid{N,T,CS,BT}}) where {N,T,CS,BT} = BT()
basis_trait(::Type{<:DiscreteGrid{N,T,CS,BT}}) where {N,T,CS,BT} = BT()
function basis_trait(::Type{<:AbstractOrthogonalGrid})
  throw(ArgumentError("`basis_trait` is undefined for `AbstractOrthogonalGrid`."))
end

Base.eltype(::MappedGrid{N,T}) where {N,T} = T
Base.eltype(::DiscreteGrid{N,T}) where {N,T} = T

# Backward-compatible identity constructor for orthogonal grids.
OrthogonalGrid(grid::OrthogonalGrid) = grid

"""
    invalidate_cell_metrics!(grid::AbstractMappedOrDiscreteGrid)

Mark cell-metric cache entries as stale.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Returns
`nothing`.
"""
function invalidate_cell_metrics!(grid::AbstractMappedOrDiscreteGrid)
  if !_has_metric_storage(grid)
    return nothing
  end
  grid.metric_caches.cell.valid = false
  return nothing
end

"""
    invalidate_face_metrics!(grid::AbstractMappedOrDiscreteGrid)

Mark face-metric cache entries as stale.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Returns
`nothing`.
"""
function invalidate_face_metrics!(grid::AbstractMappedOrDiscreteGrid)
  if !_has_metric_storage(grid)
    return nothing
  end
  grid.metric_caches.face.valid = false
  return nothing
end

"""
    refresh_cell_metrics!(grid::AbstractMappedOrDiscreteGrid; include_halo_region=false)

Recompute and return cell metrics for a mapped/discrete unified grid.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Keywords
  - `include_halo_region`: Reserved compatibility flag. Default: `false`.

# Returns
Cell metric cache payload.
"""
function refresh_cell_metrics!(
  grid::AbstractMappedOrDiscreteGrid; include_halo_region::Bool=false
)
  _require_metric_storage(grid, "refresh_cell_metrics!")
  _refresh_cell_metrics!(grid; include_halo_region=include_halo_region)
end

"""
    refresh_face_metrics!(grid::AbstractMappedOrDiscreteGrid; include_halo_region=false)

Recompute and return face metrics for a mapped/discrete unified grid.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Keywords
  - `include_halo_region`: Reserved compatibility flag. Default: `false`.

# Returns
Face metric cache payload.
"""
function refresh_face_metrics!(
  grid::AbstractMappedOrDiscreteGrid; include_halo_region::Bool=false
)
  _require_metric_storage(grid, "refresh_face_metrics!")
  _refresh_face_metrics!(grid; include_halo_region=include_halo_region)
end

"""
    cell_metrics(grid::AbstractMappedOrDiscreteGrid; refresh=false)

Access cell metric cache data.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Keywords
  - `refresh`: Recompute before returning data. Default: `false`.

# Returns
Cell metric cache payload.
"""
function cell_metrics(grid::AbstractMappedOrDiscreteGrid; refresh::Bool=false)
  _require_metric_storage(grid, "cell_metrics")
  if refresh || !grid.metric_caches.cell.valid
    return refresh_cell_metrics!(grid)
  end
  return grid.metric_caches.cell.data
end

"""
    face_metrics(grid::AbstractMappedOrDiscreteGrid; refresh=false)

Access face metric cache data.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Keywords
  - `refresh`: Recompute before returning data. Default: `false`.

# Returns
Face metric cache payload.
"""
function face_metrics(grid::AbstractMappedOrDiscreteGrid; refresh::Bool=false)
  _require_metric_storage(grid, "face_metrics")
  if refresh || !grid.metric_caches.face.valid
    return refresh_face_metrics!(grid)
  end
  return grid.metric_caches.face.data
end

function cell_metrics(::AbstractOrthogonalGrid; refresh::Bool=false)
  throw(ArgumentError("`cell_metrics` is undefined for `AbstractOrthogonalGrid`."))
end

function face_metrics(::AbstractOrthogonalGrid; refresh::Bool=false)
  throw(ArgumentError("`face_metrics` is undefined for `AbstractOrthogonalGrid`."))
end

#
# Geometry API
#

function coords(grid::Union{MappedGrid{1},DiscreteGrid{1}})
  @views grid.node_coordinates[1][grid.iterators.node.domain]
end

function coords(grid::Union{MappedGrid{2},DiscreteGrid{2}})
  @views (
    grid.node_coordinates[1][grid.iterators.node.domain],
    grid.node_coordinates[2][grid.iterators.node.domain],
  )
end

function coords(grid::Union{MappedGrid{3},DiscreteGrid{3}})
  @views (
    grid.node_coordinates[1][grid.iterators.node.domain],
    grid.node_coordinates[2][grid.iterators.node.domain],
    grid.node_coordinates[3][grid.iterators.node.domain],
  )
end

"""
    cartesian_coordinates(grid::AbstractOrthogonalGrid)

Return orthogonal-grid node coordinates in Cartesian form.

For coordinate systems already stored in Cartesian-like form, this returns
`coords(grid)` unchanged. For spherical 3D grids, this converts `(r, θ, ϕ)` to
`(x, y, z)`.
"""
@inline function cartesian_coordinates(grid::AbstractOrthogonalGrid)
  _cartesian_coordinates(coordinate_system(grid), coords(grid))
end

@inline _cartesian_coordinates(::CoordinateSystemTrait, q) = q
@inline _cartesian_coordinates(::CartesianCS, q) = q
@inline _cartesian_coordinates(::CurvilinearCS, q) = q
@inline _cartesian_coordinates(::AxisymmetricCS, q) = q
@inline _cartesian_coordinates(::CylindricalCS, q) = q

@inline function _cartesian_coordinates(
  ::SphericalCS, q::NTuple{3, <:AbstractVector{T}}
) where {T}
  r, θ, ϕ = q
  R = reshape(r, :, 1, 1)
  sinθ = reshape(sin.(θ), 1, :, 1)
  cosθ = reshape(cos.(θ), 1, :, 1)
  sinϕ = reshape(sin.(ϕ), 1, 1, :)
  cosϕ = reshape(cos.(ϕ), 1, 1, :)
  oneϕ = reshape(one.(ϕ), 1, 1, :)

  x = R .* sinθ .* cosϕ
  y = R .* sinθ .* sinϕ
  z = R .* cosθ .* oneϕ

  return x, y, z
end

function coord(grid::Union{MappedGrid{1},DiscreteGrid{1}}, (i,)::NTuple{1,Int})
  @SVector [grid.node_coordinates[1][i]]
end
function coord(grid::Union{MappedGrid{2},DiscreteGrid{2}}, (i, j)::NTuple{2,Int})
  @SVector [grid.node_coordinates[1][i, j], grid.node_coordinates[2][i, j]]
end
function coord(grid::Union{MappedGrid{3},DiscreteGrid{3}}, (i, j, k)::NTuple{3,Int})
  @SVector [
    grid.node_coordinates[1][i, j, k],
    grid.node_coordinates[2][i, j, k],
    grid.node_coordinates[3][i, j, k],
  ]
end

function centroids(grid::Union{MappedGrid{1},DiscreteGrid{1}})
  @views grid.centroid_coordinates[1][grid.iterators.cell.domain]
end

function centroids(grid::Union{MappedGrid{2},DiscreteGrid{2}})
  @views (
    grid.centroid_coordinates[1][grid.iterators.cell.domain],
    grid.centroid_coordinates[2][grid.iterators.cell.domain],
  )
end

function centroids(grid::Union{MappedGrid{3},DiscreteGrid{3}})
  @views (
    grid.centroid_coordinates[1][grid.iterators.cell.domain],
    grid.centroid_coordinates[2][grid.iterators.cell.domain],
    grid.centroid_coordinates[3][grid.iterators.cell.domain],
  )
end

function centroid(grid::Union{MappedGrid{1},DiscreteGrid{1}}, (i,)::NTuple{1,Int})
  @SVector [grid.centroid_coordinates[1][i]]
end
function centroid(grid::Union{MappedGrid{2},DiscreteGrid{2}}, (i, j)::NTuple{2,Int})
  @SVector [grid.centroid_coordinates[1][i, j], grid.centroid_coordinates[2][i, j]]
end
function centroid(grid::Union{MappedGrid{3},DiscreteGrid{3}}, (i, j, k)::NTuple{3,Int})
  @SVector [
    grid.centroid_coordinates[1][i, j, k],
    grid.centroid_coordinates[2][i, j, k],
    grid.centroid_coordinates[3][i, j, k],
  ]
end

@inline _promote_real_tuple(idx::Tuple{Vararg{Real,N}}) where {N} = promote(idx...)

@inline function _state_for_eval(grid::AbstractMappedOrDiscreteGrid)
  state = grid.state[]
  if !(state isa NamedTuple && haskey(state, :t) && haskey(state, :params))
    throw(ArgumentError("Grid state is missing `(t, params)` and cannot be evaluated."))
  end
  return state.t, state.params
end

@inline function _continuous_coord(
  grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx::Tuple{Vararg{Real,1}}
)
  ξ = _promote_real_tuple(idx)
  t, params = _state_for_eval(grid)
  return @SVector [grid.mapping_functions.x1(t, ξ..., params)]
end

@inline function _continuous_coord(
  grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx::Tuple{Vararg{Real,2}}
)
  ξη = _promote_real_tuple(idx)
  t, params = _state_for_eval(grid)
  return @SVector [
    grid.mapping_functions.x1(t, ξη..., params), grid.mapping_functions.x2(t, ξη..., params)
  ]
end

@inline function _continuous_coord(
  grid::Union{MappedGrid{3},DiscreteGrid{3}}, idx::Tuple{Vararg{Real,3}}
)
  ξηζ = _promote_real_tuple(idx)
  t, params = _state_for_eval(grid)
  return @SVector [
    grid.mapping_functions.x1(t, ξηζ..., params),
    grid.mapping_functions.x2(t, ξηζ..., params),
    grid.mapping_functions.x3(t, ξηζ..., params),
  ]
end

@inline function _continuous_forward_jacobian(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  _require_metric_functions(grid, "forward_cell_metrics")
  ξηζ = _promote_real_tuple(idx)
  t, params = _state_for_eval(grid)
  F = _as_smatrix(Val(N), grid.metric_functions_cache.forward.jacobian(t, ξηζ..., params))
  return SMatrix{N,N,T,N * N}(Tuple(F))
end

@inline function _continuous_inverse_jacobian(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  Jinv = inv(_continuous_forward_jacobian(grid, idx))
  return SMatrix{N,N,T,N * N}(Tuple(Jinv))
end

function coord(
  grid::Union{MappedGrid{N},DiscreteGrid{N}}, idx::Tuple{Vararg{Real,N}}
) where {N}
  _continuous_coord(grid, idx)
end

"""
    InverseCoordinateResult{N,T}

Result payload returned by inverse coordinate mapping routines.

# Fields
  - `coordinate`: Computational coordinate estimate.
  - `converged`: Whether the nonlinear solve converged.
  - `iterations`: Number of Newton iterations performed.
  - `residual_norm`: Infinity norm of the final physical-space residual.
"""
struct InverseCoordinateResult{N,T}
  coordinate::SVector{N,T}
  converged::Bool
  iterations::Int
  residual_norm::T
end

@inline function _computational_bounds(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}
) where {N,T}
  first_node = first(grid.iterators.global_domain.node.domain).I
  last_node = last(grid.iterators.global_domain.node.domain).I

  ξmin = ntuple(i -> T(first_node[i] - grid.nhalo), N)
  ξmax = ntuple(i -> T(last_node[i] - grid.nhalo), N)
  return SVector{N,T}(ξmin), SVector{N,T}(ξmax)
end

@inline function _physical_bounds(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}
) where {N,T}
  mins = MVector{N,T}(undef)
  maxs = MVector{N,T}(undef)
  @inbounds for n in 1:N
    mins[n] = typemax(T)
    maxs[n] = typemin(T)
  end

  @inbounds for I in grid.iterators.node.domain
    for n in 1:N
      x = grid.node_coordinates[n][I]
      mins[n] = min(mins[n], x)
      maxs[n] = max(maxs[n], x)
    end
  end

  return SVector{N,T}(mins), SVector{N,T}(maxs)
end

@inline function _seed_from_bounding_box(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}},
  x_target::SVector{N,T},
  bounds_tolerance::T,
) where {N,T}
  xmin, xmax = _physical_bounds(grid)
  ξmin, ξmax = _computational_bounds(grid)

  @inbounds for n in 1:N
    lo = xmin[n] - bounds_tolerance
    hi = xmax[n] + bounds_tolerance
    if !(lo <= x_target[n] <= hi)
      throw(
        ArgumentError(
          "Physical point lies outside the grid axis-aligned bounding box: x[$n]=$(x_target[n]) ∉ [$lo, $hi].",
        ),
      )
    end
  end

  ξ0 = ntuple(n -> begin
    span = xmax[n] - xmin[n]
    α = span > eps(T) ? (x_target[n] - xmin[n]) / span : T(0.5)
    ξmin[n] + clamp(α, zero(T), one(T)) * (ξmax[n] - ξmin[n])
  end, N)
  return SVector{N,T}(ξ0)
end

@inline function _normalize_inverse_point(::Type{T}, x::Tuple{Vararg{Real,N}}) where {N,T}
  x_promoted = _promote_real_tuple(x)
  SVector{N,T}(ntuple(i -> T(x_promoted[i]), N))
end

@inline function _normalize_inverse_guess(
  ::Type{T}, guess::Tuple{Vararg{Real,N}}
) where {N,T}
  guess_promoted = _promote_real_tuple(guess)
  SVector{N,T}(ntuple(i -> T(guess_promoted[i]), N))
end

@inline function _inverse_coordinate_result(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}},
  x_phys::Tuple{Vararg{Real,N}};
  guess::Union{Nothing,Tuple{Vararg{Real,N}}}=nothing,
  tol::Real=sqrt(eps(T)),
  maxiters::Int=25,
  damping::Bool=true,
  bounds_tolerance::Real=sqrt(eps(T)),
) where {N,T}
  if maxiters <= 0
    throw(ArgumentError("`maxiters` must be positive."))
  end
  tolT = T(tol)
  if tolT <= zero(T)
    throw(ArgumentError("`tol` must be positive."))
  end
  bounds_tolT = T(bounds_tolerance)
  if bounds_tolT < zero(T)
    throw(ArgumentError("`bounds_tolerance` must be non-negative."))
  end

  x_target = _normalize_inverse_point(T, x_phys)
  ξ = if isnothing(guess)
    _seed_from_bounding_box(grid, x_target, bounds_tolT)
  else
    _normalize_inverse_guess(T, guess)
  end

  residual_norm = typemax(T)
  converged = false
  iters = 0
  for iter in 1:maxiters
    iters = iter
    r = _continuous_coord(grid, Tuple(ξ)) - x_target
    residual_norm = norm(r, Inf)
    if residual_norm <= tolT
      converged = true
      break
    end

    Jxξ = _continuous_forward_jacobian(grid, Tuple(ξ))
    Δξ = Jxξ \ r
    ξtrial = ξ - Δξ

    if damping
      λ = one(T)
      improved = false
      while λ >= T(1) / T(64)
        ξcandidate = ξ - λ * Δξ
        rcandidate = _continuous_coord(grid, Tuple(ξcandidate)) - x_target
        if norm(rcandidate, Inf) < residual_norm
          ξtrial = ξcandidate
          improved = true
          break
        end
        λ *= T(0.5)
      end
      if !improved
        ξtrial = ξ - Δξ
      end
    end

    ξ = ξtrial
  end

  if !converged
    r = _continuous_coord(grid, Tuple(ξ)) - x_target
    residual_norm = norm(r, Inf)
    converged = residual_norm <= tolT
  end

  return InverseCoordinateResult{N,T}(ξ, converged, iters, residual_norm)
end

"""
    computational_coordinate(grid, x_phys; kwargs...)

Solve the inverse map `x(ξ) = x_phys` for mapped or discrete unified grids.

The routine uses Newton iterations on the forward mapping with Jacobian
`∂x/∂ξ`, seeded from an axis-aligned physical-space bounding box.

# Arguments
  - `grid`: Mapped or discrete unified grid.
  - `x_phys`: Physical-space point as an `NTuple{N,Real}`.

# Keywords
  - `guess`: Optional initial guess for `ξ`. Default: `nothing` (bounding-box seed).
  - `tol`: Infinity-norm residual tolerance. Default: `sqrt(eps(T))`.
  - `maxiters`: Maximum Newton iterations. Default: `25`.
  - `damping`: Enable simple backtracking damping. Default: `true`.
  - `bounds_tolerance`: Tolerance for bounding-box containment check. Default: `sqrt(eps(T))`.
  - `throw_on_failure`: Throw if nonlinear solve does not converge. Default: `true`.
  - `return_result`: Return `InverseCoordinateResult` instead of only `ξ`. Default: `false`.

# Returns
Either computational coordinate `SVector{N,T}` or `InverseCoordinateResult{N,T}`
when `return_result=true`.
"""
@inline function computational_coordinate(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}},
  x_phys::Tuple{Vararg{Real,N}};
  guess::Union{Nothing,Tuple{Vararg{Real,N}}}=nothing,
  tol::Real=sqrt(eps(T)),
  maxiters::Int=25,
  damping::Bool=true,
  bounds_tolerance::Real=sqrt(eps(T)),
  throw_on_failure::Bool=true,
  return_result::Bool=false,
) where {N,T}
  result = _inverse_coordinate_result(
    grid,
    x_phys;
    guess=guess,
    tol=tol,
    maxiters=maxiters,
    damping=damping,
    bounds_tolerance=bounds_tolerance,
  )
  if throw_on_failure && !result.converged
    throw(
      ErrorException(
        "Inverse mapping failed to converge in $(result.iterations) iterations (residual=$(result.residual_norm)).",
      ),
    )
  end
  return return_result ? result : result.coordinate
end

function computational_coordinate(grid::AbstractOrthogonalGrid, x_phys; kwargs...)
  throw(ArgumentError("`computational_coordinate` is undefined for `AbstractOrthogonalGrid`."))
end

@inline function _cell_forward_metric_at(grid::AbstractMappedOrDiscreteGrid, idx)
  cm = cell_metrics(grid)
  cm.forward[idx...]
end

@inline function _cell_inverse_metric_at(grid::AbstractMappedOrDiscreteGrid, idx)
  cm = cell_metrics(grid)
  cm.inverse[idx...]
end

@inline function _cell_forward_metric_at(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  F = _continuous_forward_jacobian(grid, idx)
  return Metric(F, det(F))
end

@inline function _cell_inverse_metric_at(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  G = _continuous_inverse_jacobian(grid, idx)
  return Metric(G, inv(det(G)))
end

"""
    forward_cell_metrics(grid, idx)

Return the forward cell metric payload at `idx` as a `Metric`.

# Arguments
  - `grid`: Mapped or discrete unified grid.
  - `idx`: `CartesianIndex`, integer tuple (discrete coordinate), or real tuple
    (continuous coordinate).
"""
@inline function forward_cell_metrics(
  grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex
)
  forward_cell_metrics(grid, idx.I)
end
@inline function forward_cell_metrics(
  grid::Union{MappedGrid,DiscreteGrid}, idx::NTuple{N,Int}
) where {N}
  _cell_forward_metric_at(grid, idx)
end
@inline function forward_cell_metrics(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  _cell_forward_metric_at(grid, idx)
end

"""
    inverse_cell_metrics(grid, idx)

Return the inverse cell metric payload at `idx` as a `Metric`.

# Arguments
  - `grid`: Mapped or discrete unified grid.
  - `idx`: `CartesianIndex`, integer tuple (discrete coordinate), or real tuple
    (continuous coordinate).
"""
@inline function inverse_cell_metrics(
  grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex
)
  inverse_cell_metrics(grid, idx.I)
end
@inline function inverse_cell_metrics(
  grid::Union{MappedGrid,DiscreteGrid}, idx::NTuple{N,Int}
) where {N}
  _cell_inverse_metric_at(grid, idx)
end
@inline function inverse_cell_metrics(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  _cell_inverse_metric_at(grid, idx)
end

function forward_cell_metrics(::AbstractOrthogonalGrid, idx)
  throw(ArgumentError("`forward_cell_metrics` is undefined for `AbstractOrthogonalGrid`."))
end

function inverse_cell_metrics(::AbstractOrthogonalGrid, idx)
  throw(ArgumentError("`inverse_cell_metrics` is undefined for `AbstractOrthogonalGrid`."))
end

"""
    cellvolume(grid, idx)

Compute cell volume at a given index.

# Arguments
  - `grid`: Unified grid instance.
  - `idx`: `CartesianIndex` or tuple index.
"""
function cellvolume(grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex)
  cellvolume(grid, idx.I)
end
function cellvolume(
  grid::Union{MappedGrid{N,T,CS,BT},DiscreteGrid{N,T,CS,BT}}, idx::NTuple{N,Int}
) where {N,T,CS,BT}
  _cellvolume_dispatch(CS(), BT(), grid, idx)
end
function cellvolume(
  grid::Union{MappedGrid{N,T,CS,BT},DiscreteGrid{N,T,CS,BT}}, idx::Tuple{Vararg{Real,N}}
) where {N,T,CS,BT}
  _cellvolume_dispatch(CS(), BT(), grid, _promote_real_tuple(idx))
end

@inline function _face_loc_axis_side(::Val{N}, loc::Symbol) where {N}
  l = Symbol(lowercase(String(loc)))
  axis, side = if l in (:ilo, :imin, :xlo, :xmin)
    (1, :lo)
  elseif l in (:ihi, :imax, :xhi, :xmax)
    (1, :hi)
  elseif l in (:jlo, :jmin, :ylo, :ymin)
    (2, :lo)
  elseif l in (:jhi, :jmax, :yhi, :ymax)
    (2, :hi)
  elseif l in (:klo, :kmin, :zlo, :zmin)
    (3, :lo)
  elseif l in (:khi, :kmax, :zhi, :zmax)
    (3, :hi)
  else
    throw(
      ArgumentError(
        "Unsupported face location `$loc`. Expected one of `:ilo/:ihi`, `:jlo/:jhi`, or `:klo/:khi`.",
      ),
    )
  end
  if axis > N
    throw(
      ArgumentError(
        "Face location `$loc` maps to axis $axis, which is invalid for $N-D grid."
      ),
    )
  end
  return axis, side
end

@inline function _face_location_symbol(axis::Int, side::Symbol, ::Val{N}) where {N}
  if axis == 1
    return side === :lo ? :ilo : :ihi
  elseif axis == 2 && N >= 2
    return side === :lo ? :jlo : :jhi
  elseif axis == 3 && N >= 3
    return side === :lo ? :klo : :khi
  end
  throw(ArgumentError("Invalid axis/side combination `(axis=$axis, side=$side)` for N=$N."))
end

@inline function _face_evaluation_coordinate(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}},
  Icell::CartesianIndex{N},
  axis::Int,
  side::Symbol,
) where {N,T}
  Iglobal = grid.iterators.global_domain.cell.full[Icell]
  half = T(0.5)
  base = ntuple(d -> T(Iglobal.I[d] - grid.nhalo) + half, N)
  offset = side === :hi ? half : -half
  return ntuple(d -> d == axis ? base[d] + offset : base[d], N)
end

@inline function _face_mapped_coordinate(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}},
  Icell::CartesianIndex{N},
  axis::Int,
  side::Symbol,
) where {N,T}
  ξηζ_face = _face_evaluation_coordinate(grid, Icell, axis, side)
  q = _continuous_coord(grid, ξηζ_face)
  return SVector{N,T}(ntuple(d -> T(q[d]), N))
end

@inline function _face_forward_jacobian(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}},
  Icell::CartesianIndex{N},
  axis::Int,
  side::Symbol,
) where {N,T}
  _require_metric_functions(grid, "face_area")
  ξηζ_face = _face_evaluation_coordinate(grid, Icell, axis, side)
  t, params = _state_for_eval(grid)
  Fraw = grid.metric_functions_cache.forward.jacobian(t, ξηζ_face..., params)
  return SMatrix{N,N,T,N * N}(Tuple(Fraw))
end

@inline function _face_coord_to_cartesian(
  ::CoordinateSystemTrait, q::SVector{N,T}
) where {N,T}
  return q
end

@inline function _face_coord_to_cartesian(::CartesianCS, q::SVector{N,T}) where {N,T}
  return q
end

@inline function _face_coord_to_cartesian(::CurvilinearCS, q::SVector{N,T}) where {N,T}
  return q
end

@inline function _face_coord_to_cartesian(::AxisymmetricCS, q::SVector{2,T}) where {T}
  return q
end

@inline function _face_coord_to_cartesian(::CylindricalCS, q::SVector{2,T}) where {T}
  return q
end

@inline function _face_coord_to_cartesian(::SphericalCS, q::SVector{3,T}) where {T}
  r, θ, ϕ = q
  sθ = sin(θ)
  cθ = cos(θ)
  sϕ = sin(ϕ)
  cϕ = cos(ϕ)
  return SVector{3,T}(r * sθ * cϕ, r * sθ * sϕ, r * cθ)
end

@inline function _face_coord_to_cartesian_jacobian(
  ::CoordinateSystemTrait, q::SVector{N,T}
) where {N,T}
  return one(SMatrix{N,N,T})
end

@inline function _face_coord_to_cartesian_jacobian(
  ::CartesianCS, q::SVector{N,T}
) where {N,T}
  return one(SMatrix{N,N,T})
end

@inline function _face_coord_to_cartesian_jacobian(
  ::CurvilinearCS, q::SVector{N,T}
) where {N,T}
  return one(SMatrix{N,N,T})
end

@inline function _face_coord_to_cartesian_jacobian(
  ::AxisymmetricCS, q::SVector{2,T}
) where {T}
  return one(SMatrix{2,2,T})
end

@inline function _face_coord_to_cartesian_jacobian(
  ::CylindricalCS, q::SVector{2,T}
) where {T}
  return one(SMatrix{2,2,T})
end

@inline function _face_coord_to_cartesian_jacobian(::SphericalCS, q::SVector{3,T}) where {T}
  r, θ, ϕ = q
  sθ = sin(θ)
  cθ = cos(θ)
  sϕ = sin(ϕ)
  cϕ = cos(ϕ)
  return SMatrix{3,3,T,9}(
    sθ * cϕ,
    sθ * sϕ,
    cθ,
    r * cθ * cϕ,
    r * cθ * sϕ,
    -r * sθ,
    -r * sθ * sϕ,
    r * sθ * cϕ,
    zero(T),
  )
end

@inline function _face_area_vector_from_jacobian(
  jacobian::SMatrix{2,2,T,4}, axis::Int
) where {T}
  if axis == 1
    return SVector{2,T}(jacobian[2, 2], -jacobian[1, 2])
  elseif axis == 2
    return SVector{2,T}(-jacobian[2, 1], jacobian[1, 1])
  end
  throw(ArgumentError("Invalid 2D face axis: $axis"))
end

@inline function _face_area_vector_from_jacobian(
  jacobian::SMatrix{3,3,T,9}, axis::Int
) where {T}
  aξ = SVector{3,T}(jacobian[1, 1], jacobian[2, 1], jacobian[3, 1])
  aη = SVector{3,T}(jacobian[1, 2], jacobian[2, 2], jacobian[3, 2])
  aζ = SVector{3,T}(jacobian[1, 3], jacobian[2, 3], jacobian[3, 3])
  if axis == 1
    return cross(aη, aζ)
  elseif axis == 2
    return cross(aζ, aξ)
  elseif axis == 3
    return cross(aξ, aη)
  end
  throw(ArgumentError("Invalid 3D face axis: $axis"))
end

@inline _face_area_scale(::CoordinateSystemTrait, q::SVector{N,T}) where {N,T} = one(T)
@inline _face_area_scale(::CurvilinearCS, q::SVector{N,T}) where {N,T} = one(T)
@inline _face_area_scale(::CartesianCS, q::SVector{N,T}) where {N,T} = one(T)
@inline _face_area_scale(::SphericalCS, q::SVector{N,T}) where {N,T} = one(T)

@inline _face_rotational_radius(::CylindricalCS, q::SVector{2,T}) where {T} = q[1]
@inline _face_rotational_radius(::AxisymmetricCS{:x}, q::SVector{2,T}) where {T} = q[2]
@inline _face_rotational_radius(::AxisymmetricCS{:y}, q::SVector{2,T}) where {T} = q[1]
function _face_rotational_radius(::AxisymmetricCS{Axis}, ::SVector{2,T}) where {Axis,T}
  throw(
    ArgumentError(
      "Unsupported axisymmetric axis `:$Axis`. Supported values are `:x` and `:y`."
    ),
  )
end

@inline function _face_area_scale(
  cs::Union{CylindricalCS,AxisymmetricCS}, q::SVector{2,T}
) where {T}
  return T(2π) * abs(_face_rotational_radius(cs, q))
end

@inline _face_outward_sign(side::Symbol, ::Type{T}) where {T} =
  side === :hi ? one(T) : -one(T)

@inline function _face_geometry(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}, loc::Symbol
) where {N,T}
  if N != 2 && N != 3
    throw(
      ArgumentError(
        "`face_area` and `outward_face_normal` are only supported for 2D and 3D unified grids.",
      ),
    )
  end
  axis, side = _face_loc_axis_side(Val(N), loc)
  Icell = CartesianIndex(idx)
  q = _face_mapped_coordinate(grid, Icell, axis, side)

  Fq = _face_forward_jacobian(grid, Icell, axis, side)

  cs = coordinate_system(grid)
  A = _face_coord_to_cartesian_jacobian(cs, q)
  Jx = SMatrix{N,N,T,N * N}(Tuple(A * Fq))
  area_vec_plus = _face_area_vector_from_jacobian(Jx, axis)
  outward_vec = area_vec_plus * _face_outward_sign(side, T)
  base_area = norm(outward_vec)
  normal = base_area > zero(T) ? outward_vec / base_area : zero(SVector{N,T})
  area = base_area * _face_area_scale(cs, q)
  x = _face_coord_to_cartesian(cs, q)

  return (; normal=normal, area=area, cartesian_coordinate=x, mapped_coordinate=q)
end

"""
    face_area(grid::Union{MappedGrid,DiscreteGrid}, idx, loc::Symbol)

Return the physical area (3D) or boundary measure (2D) of the face at cell index
`idx` and face selector `loc`.

# Arguments
  - `grid`: Mapped or discrete unified grid.
  - `idx`: Cell index as `CartesianIndex` or `NTuple{N,Int}`.
  - `loc`: Face selector symbol (for example `:ilo`, `:ihi`, `:jlo`, `:jhi`, `:klo`, `:khi`).

"""
@inline function face_area(
  grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex, loc::Symbol
)
  face_area(grid, idx.I, loc)
end
@inline function face_area(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}, loc::Symbol
) where {N,T}
  geom = _face_geometry(grid, idx, loc)
  return geom.area
end

"""
    outward_face_normal(grid::Union{MappedGrid,DiscreteGrid}, idx, loc::Symbol)

Return the outward Cartesian unit normal vector for the face at cell index `idx`
and face selector `loc`.

# Arguments
  - `grid`: Mapped or discrete unified grid.
  - `idx`: Cell index as `CartesianIndex` or `NTuple{N,Int}`.
  - `loc`: Face selector symbol (for example `:ilo`, `:ihi`, `:jlo`, `:jhi`, `:klo`, `:khi`).

"""
@inline function outward_face_normal(
  grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex, loc::Symbol
)
  outward_face_normal(grid, idx.I, loc)
end
@inline function outward_face_normal(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}, loc::Symbol
) where {N,T}
  geom = _face_geometry(grid, idx, loc)
  return geom.normal
end

function face_area(::AbstractOrthogonalGrid, idx, loc::Symbol)
  throw(ArgumentError("`face_area` is undefined for `AbstractOrthogonalGrid`."))
end

function outward_face_normal(::AbstractOrthogonalGrid, idx, loc::Symbol)
  throw(ArgumentError("`outward_face_normal` is undefined for `AbstractOrthogonalGrid`."))
end

@inline function _jacobian_volume_factor(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}
) where {N,T}
  return _cell_forward_metric_at(grid, idx).J
end
@inline function _jacobian_volume_factor(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  return T(det(_continuous_forward_jacobian(grid, idx)))
end

@inline _radial_centroid_1d(grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx::NTuple{1,Int}) = grid.centroid_coordinates[1][idx...]
@inline _radial_centroid_1d(grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx::Tuple{Vararg{Real,1}}) = _continuous_coord(
  grid, idx
)[1]

@inline function _axisymmetric_radius(
  ::AxisymmetricCS{:x}, grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx::NTuple{2,Int}
)
  return grid.centroid_coordinates[2][idx...]
end
@inline function _axisymmetric_radius(
  ::AxisymmetricCS{:x},
  grid::Union{MappedGrid{2},DiscreteGrid{2}},
  idx::Tuple{Vararg{Real,2}},
)
  return _continuous_coord(grid, idx)[2]
end

@inline function _axisymmetric_radius(
  ::AxisymmetricCS{:y}, grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx::NTuple{2,Int}
)
  return grid.centroid_coordinates[1][idx...]
end
@inline function _axisymmetric_radius(
  ::AxisymmetricCS{:y},
  grid::Union{MappedGrid{2},DiscreteGrid{2}},
  idx::Tuple{Vararg{Real,2}},
)
  return _continuous_coord(grid, idx)[1]
end

@inline function _axisymmetric_radius(
  ::AxisymmetricCS{Axis}, grid::AbstractMappedOrDiscreteGrid, idx::Tuple{Vararg{Real,2}}
) where {Axis}
  throw(
    ArgumentError(
      "Unsupported axisymmetric axis `:$Axis`. Supported values are `:x` and `:y`."
    ),
  )
end

@inline function _cellvolume_dispatch(
  ::CoordinateSystemTrait,
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::NTuple{N,Int},
) where {N}
  return _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CoordinateSystemTrait,
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,N}},
) where {N}
  return _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS, ::CartesianBasis, grid::AbstractMappedOrDiscreteGrid, idx::NTuple{1,Int}
)
  r = _radial_centroid_1d(grid, idx)
  return (2π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS,
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,1}},
)
  r = _radial_centroid_1d(grid, idx)
  return (2π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS,
  ::CartesianBasis,
  grid::Union{MappedGrid{2},DiscreteGrid{2}},
  idx::NTuple{2,Int},
)
  r = centroid(grid, idx)[1]
  return (2π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS,
  ::CartesianBasis,
  grid::Union{MappedGrid{2},DiscreteGrid{2}},
  idx::Tuple{Vararg{Real,2}},
)
  r = _continuous_coord(grid, idx)[1]
  return (2π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  cs::AxisymmetricCS{Axis},
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::NTuple{2,Int},
) where {Axis}
  r = _axisymmetric_radius(cs, grid, idx)
  return (2π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  cs::AxisymmetricCS{Axis},
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,2}},
) where {Axis}
  r = _axisymmetric_radius(cs, grid, idx)
  return (2π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, grid::AbstractMappedOrDiscreteGrid, idx::NTuple{1,Int}
)
  r = _radial_centroid_1d(grid, idx)
  return (4π * r^2) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,1}},
)
  r = _radial_centroid_1d(grid, idx)
  return (4π * r^2) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::Union{MappedGrid{3},DiscreteGrid{3}},
  idx::NTuple{3,Int},
)
  c = centroid(grid, idx)
  r = c[1]
  θ = c[2]
  return (r^2 * sin(θ)) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::Union{MappedGrid{3},DiscreteGrid{3}},
  idx::Tuple{Vararg{Real,3}},
)
  c = _continuous_coord(grid, idx)
  r = c[1]
  θ = c[2]
  return (r^2 * sin(θ)) * _jacobian_volume_factor(grid, idx)
end

function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, ::AbstractMappedOrDiscreteGrid, idx::NTuple{2,Int}
)
  throw(ArgumentError("Spherical cellvolume dispatch is not implemented for 2D grids."))
end

function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  ::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,2}},
)
  throw(ArgumentError("Spherical cellvolume dispatch is not implemented for 2D grids."))
end

function _cellvolume_dispatch(
  cs::CoordinateSystemTrait,
  bt::SphericalBasis,
  ::AbstractMappedOrDiscreteGrid,
  ::NTuple{N,Int},
) where {N}
  throw(
    ArgumentError(
      "Unsupported basis/coordinate-system combination: $(typeof(bt)) with $(typeof(cs))."
    ),
  )
end

function _cellvolume_dispatch(
  cs::CoordinateSystemTrait,
  bt::SphericalBasis,
  ::AbstractMappedOrDiscreteGrid,
  ::Tuple{Vararg{Real,N}},
) where {N}
  throw(
    ArgumentError(
      "Unsupported basis/coordinate-system combination: $(typeof(bt)) with $(typeof(cs))."
    ),
  )
end

function cellvolumes(grid::Union{MappedGrid,DiscreteGrid})
  volumes = zeros(eltype(grid), size(grid.iterators.cell.domain))

  for (idx0, idx1) in zip(CartesianIndices(volumes), grid.iterators.cell.domain)
    volumes[idx0] = cellvolume(grid, idx1)
  end

  return volumes
end

function cellvolumes(grid::OrthogonalGrid)
  @views grid.cell_volumes[grid.iterators.cell.domain]
end

cellsize(grid::Union{MappedGrid,DiscreteGrid}) = size(grid.iterators.cell.domain)
cellsize(grid::OrthogonalGrid) = size(grid.iterators.cell.domain)

cellsize_withhalo(grid::Union{MappedGrid,DiscreteGrid}) = size(grid.iterators.cell.full)
cellsize_withhalo(grid::OrthogonalGrid) = size(grid.iterators.cell.full)

"""
    jacobian_matrix(grid, idx)

Return the forward Jacobian matrix at a given cell index.

# Arguments
  - `grid`: Mapped or discrete unified grid.
  - `idx`: `CartesianIndex`, integer tuple, or real tuple.
"""
function jacobian_matrix(grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex)
  jacobian_matrix(grid, idx.I)
end
function jacobian_matrix(grid::Union{MappedGrid,DiscreteGrid}, idx::NTuple{N,Int}) where {N}
  _cell_forward_metric_at(grid, idx).jacobian_matrix
end
function jacobian_matrix(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  _continuous_forward_jacobian(grid, idx)
end
