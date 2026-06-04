#
# Inverse coordinates
#

#
# Inverse coordinate result
#

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
#
# Inverse coordinate helpers
#

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

#
# Public inverse coordinate API
#

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
  throw(
    ArgumentError("`computational_coordinate` is undefined for `AbstractOrthogonalGrid`.")
  )
end
