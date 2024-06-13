
using KernelAbstractions

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

function RThetaPhiGrid(
  (r0, θ0, ϕ0),
  (r1, θ1, ϕ1),
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
  r1d = range(r0, r1; length=ni)
  θ1d = range(θ0, θ1; length=nj)
  ϕ1d = range(ϕ0, ϕ1; length=nk)

  @inbounds for k in 1:nk
    for j in 1:nj
      for i in 1:ni
        x[i, j, k] = r1d[i] * sin(θ1d[j]) * cos(ϕ1d[k])
        y[i, j, k] = r1d[i] * sin(θ1d[j]) * sin(ϕ1d[k])
        z[i, j, k] = r1d[i] * cos(θ1d[j])
      end
    end
  end

  return CurvilinearGrid3D(x, y, z, nhalo; backend=backend)
end

function RThetaPhiGrid(
  r::AbstractVector{T},
  θ::AbstractVector{T},
  ϕ::AbstractVector{T},
  nhalo::Int,
  backend=CPU(),
) where {T}
  ni = length(r)
  nj = length(θ)
  nk = length(ϕ)

  x = zeros(T, ni, nj, nk)
  y = zeros(T, ni, nj, nk)
  z = zeros(T, ni, nj, nk)

  @inbounds for k in 1:nk
    for j in 1:nj
      for i in 1:ni
        x[i, j, k] = r[i] * sin(θ[j]) * cos(ϕ[k])
        y[i, j, k] = r[i] * sin(θ[j]) * sin(ϕ[k])
        z[i, j, k] = r[i] * cos(θ[j])
      end
    end
  end

  return CurvilinearGrid3D(x, y, z, nhalo; backend=backend)
end
