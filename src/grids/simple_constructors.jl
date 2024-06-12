
using KernelAbstractions

function RectlinearGrid(x⃗0, x⃗1, (ni_cells,)::NTuple{1,Int}, nhalo::Int, backend=CPU())
  ni = ni_cells + 1

  ni_nodes = ni + 2nhalo

  x = zeros(ni_nodes, nj_nodes)

  # node positions (non-halo)
  x1d = range(x⃗0[1], x⃗1[1]; length=ni)

  for i in 1:ni
    x[i + nhalo, j + nhalo] = x1d[i]
  end

  return CurvilinearGrid1D(x, nhalo; backend=backend)
end

function RectlinearGrid(
  x⃗0, x⃗1, (ni_cells, nj_cells)::NTuple{2,Int}, nhalo::Int, backend=CPU()
)
  ni = ni_cells + 1
  nj = nj_cells + 1

  ni_nodes = ni + 2nhalo
  nj_nodes = nj + 2nhalo

  x = zeros(ni_nodes, nj_nodes)
  y = zeros(ni_nodes, nj_nodes)

  # node positions (non-halo)
  x1d = range(x⃗0[1], x⃗1[1]; length=ni)
  y1d = range(x⃗0[2], x⃗1[2]; length=nj)

  for j in 1:nj
    for i in 1:ni
      x[i + nhalo, j + nhalo] = x1d[i]
      y[i + nhalo, j + nhalo] = y1d[j]
    end
  end

  return CurvilinearGrid2D(x, y, nhalo; backend=backend)
end

function RectlinearGrid(
  x⃗0, x⃗1, (ni_cells, nj_cells, nk_cells)::NTuple{3,Int}, nhalo::Int, backend=CPU()
)
  ni = ni_cells + 1
  nj = nj_cells + 1
  nk = nk_cells + 1

  ni_nodes = ni + 2nhalo
  nj_nodes = nj + 2nhalo
  nk_nodes = nk + 2nhalo

  x = zeros(ni_nodes, nj_nodes, nk_nodes)
  y = zeros(ni_nodes, nj_nodes, nk_nodes)
  z = zeros(ni_nodes, nj_nodes, nk_nodes)

  # node positions (non-halo)
  x1d = range(x⃗0[1], x⃗1[1]; length=ni)
  y1d = range(x⃗0[2], x⃗1[2]; length=nj)
  z1d = range(x⃗0[3], x⃗1[3]; length=nk)

  for k in 1:nk
    for j in 1:nj
      for i in 1:ni
        x[i + nhalo, j + nhalo, k + nhalo] = x1d[i]
        y[i + nhalo, j + nhalo, k + nhalo] = y1d[j]
        z[i + nhalo, j + nhalo, k + nhalo] = z1d[k]
      end
    end
  end

  return CurvilinearGrid3D(x, y, z, nhalo; backend=backend)
end

function CylindricalGrid(
  x⃗0, x⃗1, (ni_cells, nj_cells)::NTuple{2,Int}, nhalo::Int, backend=CPU()
)
  ni = ni_cells + 1
  nj = nj_cells + 1

  ni_nodes = ni + 2nhalo
  nj_nodes = nj + 2nhalo

  x = zeros(ni_nodes, nj_nodes)
  y = zeros(ni_nodes, nj_nodes)

  # node positions (non-halo)
  r1d = range(x⃗0[1], x⃗1[1]; length=ni)
  θ1d = range(x⃗0[2], x⃗1[2]; length=nj)

  for j in 1:nj
    for i in 1:ni
      x[i + nhalo, j + nhalo] = r1d[i] * cos(θ1d[j])
      y[i + nhalo, j + nhalo] = r1d[i] * sin(θ1d[j])
    end
  end

  return CurvilinearGrid2D(x, y, nhalo; backend=backend)
end

function CylindricalGrid(r, θ, nhalo::Int, backend=CPU()) where {T}
  ni = length(r)
  nj = length(θ)

  ni_nodes = ni + 2nhalo
  nj_nodes = nj + 2nhalo

  x = zeros(ni_nodes, nj_nodes)
  y = zeros(ni_nodes, nj_nodes)

  for j in 1:nj
    for i in 1:ni
      x[i + nhalo, j + nhalo] = r[i] * cos(θ[j])
      y[i + nhalo, j + nhalo] = r[i] * sin(θ[j])
    end
  end

  return CurvilinearGrid2D(x, y, nhalo; backend=backend)
end
