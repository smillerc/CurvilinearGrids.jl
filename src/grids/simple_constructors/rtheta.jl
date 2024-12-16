"""
    rtheta_grid((r0, θ0), (r1, θ1), (ni_cells, nj_cells), discretization_scheme::Symbol, backend=CPU(), T=Float64) -> CurvilinearGrid2D

Create an equally spaced polar grid based on `r` and `θ`
"""
function rtheta_grid(
  (r0, θ0),
  (r1, θ1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  discretization_scheme::Symbol,
  backend=CPU(),
  T=Float64;
  is_static=true,
)
  ni = ni_cells + 1
  nj = nj_cells + 1

  x = zeros(T, ni, nj)
  y = zeros(T, ni, nj)

  # node positions (non-halo)
  r1d = range(r0, r1; length=ni)
  θ1d = range(θ0, θ1; length=nj)

  @inbounds for j in 1:nj
    for i in 1:ni
      x[i, j] = r1d[i] * cos(θ1d[j])
      y[i, j] = r1d[i] * sin(θ1d[j])
    end
  end

  return CurvilinearGrid2D(
    x, y, discretization_scheme; backend=backend, is_static=is_static
  )
end

"""
    rtheta_grid(r, θ, discretization_scheme::Symbol; backend=CPU()) -> CurvilinearGrid2D

Create polar grid based on vectors of `r` and `θ` coordinates
"""
function rtheta_grid(
  r::AbstractVector{T},
  θ::AbstractVector{T},
  discretization_scheme::Symbol;
  backend=CPU(),
  is_static=true,
) where {T}
  @assert all(r .>= 0) "Radius coordinates must be >= 0"

  ni = length(r)
  nj = length(θ)

  if !all(diff(r) .> 0)
    error("Invalid r vector, spacing between vertices must be > 0 everywhere")
  end

  if !all(diff(θ) .> 0)
    error("Invalid θ vector, spacing between vertices must be > 0 everywhere")
  end

  x = zeros(T, ni, nj)
  y = zeros(T, ni, nj)

  for j in 1:nj
    for i in 1:ni
      x[i, j] = r[i] * cos(θ[j])
      y[i, j] = r[i] * sin(θ[j])
    end
  end

  return CurvilinearGrid2D(
    x, y, discretization_scheme; backend=backend, is_static=is_static
  )
end

"""
    rtheta_grid((r0, θ0), (r1, θ1), (ni_cells, nj_cells), discretization_scheme, snap_to_axis, rotational_axis, backend=CPU(), T=Float64) -> AxisymmetricGrid2D

Create an equally spaced axisymmetric polar grid based on `r` and `θ`.
The axis of rotation is set by `rotational_axis` as `:x` or `:y`
"""
function axisymmetric_rtheta_grid(
  (r0, θ0),
  (r1, θ1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  discretization_scheme::Symbol,
  snap_to_axis::Bool,
  rotational_axis::Symbol;
  backend=CPU(),
  T=Float64,
  is_static=true,
)
  ni = ni_cells + 1
  nj = nj_cells + 1

  x = zeros(T, ni, nj)
  y = zeros(T, ni, nj)

  if rotational_axis === :R || rotational_axis === :x
    axis = :x # equator axis
  else # :z, :Z, etc...
    axis = :y # pole axis
  end

  # node positions (non-halo)
  r1d = range(r0, r1; length=ni)
  θ1d = range(θ0, θ1; length=nj)

  @inbounds for j in 1:nj
    for i in 1:ni
      x[i, j] = r1d[i] * cos(θ1d[j])
      y[i, j] = r1d[i] * sin(θ1d[j])
    end
  end

  return AxisymmetricGrid2D(
    x, y, discretization_scheme, snap_to_axis, axis; backend=backend, is_static=is_static
  )
end

"""
    axisymmetric_rtheta_grid(r, θ, discretization_scheme, snap_to_axis, rotational_axis::Symbol, backend=CPU()) -> AxisymmetricGrid2D

Create polar grid based on vectors of `r` and `θ` coordinates.
The axis of rotation is set by `rotational_axis` as `:x` or `:y`
"""
function axisymmetric_rtheta_grid(
  r::AbstractVector{T},
  θ::AbstractVector{T},
  discretization_scheme::Symbol,
  snap_to_axis::Bool,
  rotational_axis::Symbol;
  backend=CPU(),
  is_static=true,
) where {T}
  @assert all(r .>= 0) "Radius coordinates must be >= 0"

  ni = length(r)
  nj = length(θ)

  if !all(diff(r) .> 0)
    error("Invalid r vector, spacing between vertices must be > 0 everywhere")
  end

  if !all(diff(θ) .> 0)
    error("Invalid θ vector, spacing between vertices must be > 0 everywhere")
  end

  x = zeros(T, ni, nj)
  y = zeros(T, ni, nj)

  if rotational_axis === :R || rotational_axis === :x
    axis = :x # equator axis
  else # :z, :Z, etc...
    axis = :y # pole axis
  end

  for j in 1:nj
    for i in 1:ni
      x[i, j] = r[i] * cos(θ[j])
      y[i, j] = r[i] * sin(θ[j])
    end
  end

  return AxisymmetricGrid2D(
    x, y, discretization_scheme, snap_to_axis, axis; backend=backend, is_static=is_static
  )
end
