function RectlinearGrid(x0, x1, ncells, nhalo::Int, backend=CPU(), T=Float64)
  ni = ncells + 1
  x = collect(T, range(x0, x1; length=ni))
  return CurvilinearGrid1D(x, nhalo; backend=backend)
end

function RectlinearGrid(
  (x0, y0),
  (x1, y1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  nhalo::Int,
  backend=CPU(),
  T=Float64,
)
  @assert !(x0 ≈ x1) "The endpoints x0 and x1 are the same"
  @assert !(y0 ≈ y1) "The endpoints y0 and y1 are the same"

  ni = ni_cells + 1
  nj = nj_cells + 1

  x = zeros(T, ni, nj)
  y = zeros(T, ni, nj)

  x1d = range(x0, x1; length=ni)
  y1d = range(y0, y1; length=nj)

  @inbounds for j in 1:nj
    for i in 1:ni
      x[i, j] = x1d[i]
      y[i, j] = y1d[j]
    end
  end

  return CurvilinearGrid2D(x, y, nhalo; backend=backend)
end

function RectlinearGrid(
  x::AbstractVector{T}, y::AbstractVector{T}, nhalo::Int, backend=CPU()
) where {T}
  ni = length(x)
  nj = length(y)

  @assert ni >= 2 "The x vector must have more than 2 points"
  @assert nj >= 2 "The y vector must have more than 2 points"

  x2d = zeros(T, ni, nj)
  y2d = zeros(T, ni, nj)

  @inbounds for j in 1:nj
    for i in 1:ni
      x2d[i, j] = x[i]
      y2d[i, j] = y[j]
    end
  end

  return CurvilinearGrid2D(x2d, y2d, nhalo; backend=backend)
end

function AxisymmetricRectlinearGrid(
  (x0, y0),
  (x1, y1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  nhalo::Int,
  snap_to_axis::Bool,
  rotational_axis::Symbol,
  backend=CPU(),
  T=Float64,
)
  @assert !(x0 ≈ x1) "The endpoints x0 and x1 are the same"
  @assert !(y0 ≈ y1) "The endpoints y0 and y1 are the same"

  ni = ni_cells + 1
  nj = nj_cells + 1

  x = zeros(T, ni, nj)
  y = zeros(T, ni, nj)

  x1d = range(x0, x1; length=ni)
  y1d = range(y0, y1; length=nj)

  @inbounds for j in 1:nj
    for i in 1:ni
      x[i, j] = x1d[i]
      y[i, j] = y1d[j]
    end
  end

  return AxisymmetricGrid2D(x, y, nhalo, snap_to_axis, rotational_axis; backend=backend)
end

function RectlinearGrid(
  (x0, y0, z0),
  (x1, y1, z1),
  (ni_cells, nj_cells, nk_cells)::NTuple{3,Int},
  nhalo::Int,
  backend=CPU(),
  T=Float64,
)
  ni = ni_cells + 1
  nj = nj_cells + 1
  nk = nk_cells + 1

  x = zeros(T, ni, nj, nk)
  y = zeros(T, ni, nj, nk)
  z = zeros(T, ni, nj, nk)

  # node positions (non-halo)
  x1d = range(x0, x1; length=ni)
  y1d = range(y0, y1; length=nj)
  z1d = range(z0, z1; length=nk)

  @inbounds for k in 1:nk
    for j in 1:nj
      for i in 1:ni
        x[i, j, k] = x1d[i]
        y[i, j, k] = y1d[j]
        z[i, j, k] = z1d[k]
      end
    end
  end

  return CurvilinearGrid3D(x, y, z, nhalo; backend=backend)
end

function RectlinearGrid(
  x::AbstractVector{T},
  y::AbstractVector{T},
  z::AbstractVector{T},
  nhalo::Int,
  backend=CPU(),
) where {T}
  ni = length(x)
  nj = length(y)
  nk = length(z)

  @assert ni >= 2 "The x vector must have more than 2 points"
  @assert nj >= 2 "The y vector must have more than 2 points"
  @assert nk >= 2 "The z vector must have more than 2 points"

  x3d = zeros(T, ni, nj, nk)
  y3d = zeros(T, ni, nj, nk)
  z3d = zeros(T, ni, nj, nk)

  @inbounds for k in 1:nk
    for j in 1:nj
      for i in 1:ni
        x3d[i, j, k] = x[i]
        y3d[i, j, k] = y[j]
        z3d[i, j, k] = z[k]
      end
    end
  end

  return CurvilinearGrid3D(x3d, y3d, z3d, nhalo; backend=backend)
end
