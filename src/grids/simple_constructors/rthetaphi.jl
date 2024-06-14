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
