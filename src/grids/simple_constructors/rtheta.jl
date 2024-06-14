"""
    RThetaGrid((r0, θ0), (r1, θ1), (ni_cells, nj_cells), nhalo::Int, backend=CPU(), T=Float64)

Create an equally spaced polar grid based on `r` and `θ`
"""
function RThetaGrid(
  (r0, θ0),
  (r1, θ1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  nhalo::Int,
  backend=CPU(),
  T=Float64,
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

  return CurvilinearGrid2D(x, y, nhalo; backend=backend)
end

"""
    RThetaGrid(r, θ, nhalo::Int, backend=CPU())

Create polar grid based on vectors of `r` and `θ` coordinates
"""
function RThetaGrid(
  r::AbstractVector{T}, θ::AbstractVector{T}, nhalo::Int, backend=CPU()
) where {T}
  @assert all(r .>= 0) "Radius coordinates must be >= 0"

  ni = length(r)
  nj = length(θ)

  x = zeros(T, ni, nj)
  y = zeros(T, ni, nj)

  for j in 1:nj
    for i in 1:ni
      x[i, j] = r[i] * cos(θ[j])
      y[i, j] = r[i] * sin(θ[j])
    end
  end

  return CurvilinearGrid2D(x, y, nhalo; backend=backend)
end

"""
    RThetaGrid((r0, θ0), (r1, θ1), (ni_cells, nj_cells), nhalo, snap_to_axis, rotational_axis, backend=CPU(), T=Float64)

Create an equally spaced axisymmetric polar grid based on `r` and `θ`.
The axis of rotation is set by `rotational_axis` as `:x` or `:y`
"""
function RThetaGrid(
  (r0, θ0),
  (r1, θ1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  nhalo::Int,
  snap_to_axis::Bool,
  rotational_axis::Symbol,
  backend=CPU(),
  T=Float64,
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

  return AxisymmetricGrid2D(x, y, nhalo, snap_to_axis, axis; backend=backend)
end

"""
    AxisymmetricRThetaGrid(r, θ, nhalo, snap_to_axis, rotational_axis::Symbol, backend=CPU())

Create polar grid based on vectors of `r` and `θ` coordinates.
The axis of rotation is set by `rotational_axis` as `:x` or `:y`
"""
function AxisymmetricRThetaGrid(
  r::AbstractVector{T},
  θ::AbstractVector{T},
  nhalo::Int,
  snap_to_axis::Bool,
  rotational_axis::Symbol,
  backend=CPU(),
) where {T}
  @assert all(r .>= 0) "Radius coordinates must be >= 0"

  ni = length(r)
  nj = length(θ)

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

  return AxisymmetricGrid2D(x, y, nhalo, snap_to_axis, axis; backend=backend)
end
