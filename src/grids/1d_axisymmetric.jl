
struct CylindricalGrid1D{CO,CE,NV,EM,CM,DL,CI,DS,CF,JF} <: AbstractCurvilinearGrid
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::Int
  domain_limits::DL
  iterators::CI
  dθ::Float64
  discretization_scheme::DS
  _coordinate_funcs::CF
  _jacobian_matrix_func::JF
end

struct SphericalGrid1D{CO,CE,NV,EM,CM,CI,DS,CF,JF} <: AbstractCurvilinearGrid
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::Int
  # domain_limits::DL
  iterators::CI
  dθ::Float64
  dϕ::Float64
  discretization_scheme::DS
  _coordinate_funcs::CF
  _jacobian_matrix_func::JF
  _nhalo::Int
end

AxisymmetricGrid1D = Union{CylindricalGrid1D,SphericalGrid1D}

function SphericalGrid1D(r::Function, ni::Int, nhalo; T=Float64, backend=CPU())

  # TODO: add asserts in all grid constructors
  @assert ni >= 2 "You need at least 2 grid points to define the mesh (ni < 2)"

  dim = 1
  check_nargs(r, dim, :r)
  test_coord_func(r, dim, :r)

  # Make a full 3d grid with only 1 cell in θ.
  # This is cheap and very useful for certain applications

  # leave the π off so we can use the more accurate cospi/sinpi functions
  # we want the cell to stradle the x-axis

  # δ = 1
  # θ0, θ1 = 1 / 4, 3 / 4 # center θ on the x-axis
  # ϕ0, ϕ1 = -1 / 4, 1 / 4 # center ϕ on the x-axis
  θ0, θ1 = deg2rad.((89.5, 90.5)) ./ π # center θ on the x-axis
  ϕ0, ϕ1 = deg2rad.((-0.5, 0.5)) ./ π # center ϕ on the x-axis

  dϕ = abs(ϕ1 - ϕ0) * π
  dθ = abs(θ1 - θ0) * π
  # if the θ is 0:π and ϕ is 0:2π, the grid metrics will be NaNs, 
  # since the Jacobian will have zero terms (and the inv() will produce NaNs)

  nj, nk = (4, 4)
  θ(j) = @. θ0 + (θ1 - θ0) * ((j - 1) / (nj - 1))
  ϕ(k) = @. ϕ0 + (ϕ1 - ϕ0) * ((k - 1) / (nk - 1))

  X3d(i, j, k) = r(i) * sinpi(θ(j)) * cospi(ϕ(k))
  Y3d(i, j, k) = r(i) * sinpi(θ(j)) * sinpi(ϕ(k))
  Z3d(i, j, k) = r(i) * cospi(θ(j))

  XYZ(i, j, k) = @SVector [
    X3d(i, j, k), #
    Y3d(i, j, k), # 
    Z3d(i, j, k), #
  ]
  function jacobian_matrix_func(ξ, η, ζ, t=0)
    return ForwardDiff.jacobian(x -> XYZ(x[1], x[2], x[3]), @SVector [ξ, η, ζ])
  end

  nnodes = (ni, nj, nk)
  ncells = nnodes .- 1
  # ni_cells = ni - 1
  # nj_cells = 1
  # nk_cells = 1
  # lo = nhalo + 1

  # limits = (
  #   node=(ilo=lo, ihi=ni + nhalo, jlo=lo, jhi=nj + nhalo, klo=lo, khi=nk + nhalo),
  #   cell=(
  #     ilo=lo,
  #     ihi=ni_cells + nhalo,
  #     jlo=lo,
  #     jhi=nj_cells + nhalo,
  #     klo=lo,
  #     khi=nk_cells + nhalo,
  #   ),
  # )

  nodeCI = CartesianIndices(nnodes .+ 2nhalo)
  cellCI = CartesianIndices(ncells .+ 2nhalo)

  node_iter, cell_iter = get_node_cell_iterators(nodeCI, cellCI, nhalo)

  celldims = size(cell_iter.full)
  nodedims = size(node_iter.full)

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  coordinate_funcs = (r=r, XYZ=XYZ, x=X3d, y=Y3d, z=Z3d)

  centroids = (
    r=KernelAbstractions.zeros(backend, T, celldims[1]),
    xyz=StructArray(;
      x=KernelAbstractions.zeros(backend, T, celldims),
      y=KernelAbstractions.zeros(backend, T, celldims),
      z=KernelAbstractions.zeros(backend, T, celldims),
    ),
  )

  _axisym_centroid_coordinates!(
    centroids, coordinate_funcs, cell_iter.full.indices[1], nhalo
  )
  _centroid_coordinates!(centroids.xyz, coordinate_funcs, cell_iter.full, nhalo)

  coords = (
    r=KernelAbstractions.zeros(backend, T, nodedims[1]),
    xyz=StructArray(;
      x=KernelAbstractions.zeros(backend, T, nodedims),
      y=KernelAbstractions.zeros(backend, T, nodedims),
      z=KernelAbstractions.zeros(backend, T, nodedims),
    ),
  )
  _axisym_node_coordinates!(coords, coordinate_funcs, node_iter.full.indices[1], nhalo)
  _node_coordinates!(coords.xyz, coordinate_funcs, node_iter.full, nhalo)

  node_velocities = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))

  discr_scheme = MetricDiscretizationSchemes.MonotoneExplicit6thOrderDiscretization(
    cell_iter.full
  )

  inner_cell_dom = expand(cell_iter.domain, -1)
  inner_node_dom = expand(node_iter.domain, -1)

  domain_iterators = (node=node_iter, cell=cell_iter)
  # domain_iterators = (
  #   node=(full=node_iter.full, domain=inner_node_dom),
  #   cell=(full=cell_iter.full, domain=inner_cell_dom),
  # )

  m = SphericalGrid1D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes[1],
    domain_iterators,
    dθ, # θ wedge angle in radians (with π included)
    dϕ, # ϕ wedge angle in radians (with π included)
    discr_scheme,
    coordinate_funcs,
    jacobian_matrix_func,
    1,
  )

  update_metrics!(m)
  check_for_invalid_metrics(m)
  return m
end

function CylindricalGrid1D(r::Function, ni::Int, nhalo; T=Float64, backend=CPU())

  # TODO: add asserts in all grid constructors
  @assert ni >= 2 "You need at least 2 grid points to define the mesh (ni < 2)"

  dim = 1
  check_nargs(r, dim, :r)
  test_coord_func(r, dim, :r)

  # Make a full 3d grid with only 1 cell in θ.
  # This is cheap and very useful for certain applications

  # leave the π off so we can use the more accurate cospi/sinpi functions
  # we want the cell to stradle the x-axis
  θ0, θ1 = 1 / 4, 3 / 4 # center θ on the x-axis
  z0, z1 = 0, 1
  nj, nk = (2, 2)

  θ(j) = @. θ0 + (θ1 - θ0) * ((j - 1) / (nj - 1))
  z(k) = @. z0 + (z1 - z0) * ((k - 1) / (nk - 1))

  dθ = abs(θ1 - θ0) * π
  # if the θ is 0:π and ϕ is 0:2π, the grid metrics will be NaNs, 
  # since the Jacobian will have zero terms (and the inv() will produce NaNs)

  X3d(i, j, k) = r(i) * cospi(θ(j))
  Y3d(i, j, k) = r(i) * sinpi(θ(j))
  Z3d(i, j, k) = z(k)

  XYZ(i, j, k) = @SVector [
    X3d(i, j, k), #
    Y3d(i, j, k), # 
    Z3d(i, j, k), #
  ]
  function jacobian_matrix_func(ξ, η, ζ, t=0)
    return ForwardDiff.jacobian(x -> XYZ(x[1], x[2], x[3]), @SVector [ξ, η, ζ])
  end

  nnodes = (ni, nj, nk)
  ncells = nnodes .- 1
  ni_cells = ni - 1
  nj_cells = 1
  nk_cells = 1
  lo = nhalo + 1

  limits = (
    node=(ilo=lo, ihi=ni + nhalo, jlo=lo, jhi=nj + nhalo, klo=lo, khi=nk + nhalo),
    cell=(
      ilo=lo,
      ihi=ni_cells + nhalo,
      jlo=lo,
      jhi=nj_cells + nhalo,
      klo=lo,
      khi=nk_cells + nhalo,
    ),
  )

  nodeCI = CartesianIndices(nnodes .+ 2nhalo)
  cellCI = CartesianIndices(ncells .+ 2nhalo)

  node_iter, cell_iter = get_node_cell_iterators(nodeCI, cellCI, nhalo)

  celldims = size(cell_iter.full)
  nodedims = size(node_iter.full)

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  coordinate_funcs = (r=r, XYZ=XYZ, x=X3d, y=Y3d, z=Z3d)

  centroids = (
    r=KernelAbstractions.zeros(backend, T, celldims[1]),
    xyz=StructArray(;
      x=KernelAbstractions.zeros(backend, T, celldims),
      y=KernelAbstractions.zeros(backend, T, celldims),
      z=KernelAbstractions.zeros(backend, T, celldims),
    ),
  )

  # centroids = StructArray((r=KernelAbstractions.zeros(backend, T, celldims[1]),))
  _axisym_centroid_coordinates!(
    centroids, coordinate_funcs, cell_iter.full.indices[1], nhalo
  )
  _centroid_coordinates!(centroids.xyz, coordinate_funcs, cell_iter.full, nhalo)

  coords = (
    r=KernelAbstractions.zeros(backend, T, nodedims[1]),
    xyz=StructArray(;
      x=KernelAbstractions.zeros(backend, T, nodedims),
      y=KernelAbstractions.zeros(backend, T, nodedims),
      z=KernelAbstractions.zeros(backend, T, nodedims),
    ),
  )
  _axisym_node_coordinates!(coords, coordinate_funcs, node_iter.full.indices[1], nhalo)
  _node_coordinates!(coords.xyz, coordinate_funcs, node_iter.full, nhalo)

  node_velocities = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))

  discr_scheme = MetricDiscretizationSchemes.MonotoneExplicit6thOrderDiscretization(
    cell_iter.full
  )

  domain_iterators = (node=node_iter, cell=cell_iter)
  m = CylindricalGrid1D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes[1],
    limits,
    domain_iterators,
    dθ, # θ wedge angle in radians (with π included)
    discr_scheme,
    coordinate_funcs,
    jacobian_matrix_func,
  )

  update_metrics!(m)
  # check_for_invalid_metrics(m)
  return m
end

function update_metrics!(m::AxisymmetricGrid1D, t::Real=0)
  # cell metrics

  domain = m.iterators.cell.domain

  @inbounds for idx in m.iterators.cell.full
    cell_idx = idx.I .+ 0.5
    @unpack J, ξ, η, ζ = metrics(m, cell_idx, t)
    m.cell_center_metrics.ξ.x[idx] = ξ.x
    m.cell_center_metrics.ξ.y[idx] = ξ.y
    m.cell_center_metrics.ξ.z[idx] = ξ.z
    m.cell_center_metrics.ξ.t[idx] = ξ.t

    m.cell_center_metrics.η.x[idx] = η.x
    m.cell_center_metrics.η.y[idx] = η.y
    m.cell_center_metrics.η.z[idx] = η.z
    m.cell_center_metrics.η.t[idx] = η.t

    m.cell_center_metrics.ζ.x[idx] = ζ.x
    m.cell_center_metrics.ζ.y[idx] = ζ.y
    m.cell_center_metrics.ζ.z[idx] = ζ.z
    m.cell_center_metrics.ζ.t[idx] = ζ.t

    m.cell_center_metrics.J[idx] = J
  end

  θ_domain = expand(m.iterators.cell.domain, 2, +1)
  θϕ_domain = expand(θ_domain, 3, +1)
  full_m1 = expand(m.iterators.cell.full, -1)
  m.centroid_coordinates.xyz[full_m1]
  @show full_m1
  MetricDiscretizationSchemes.update_metrics!(
    m.discretization_scheme,
    m.centroid_coordinates.xyz,
    m.cell_center_metrics,
    m.edge_metrics,
    full_m1,
  )

  # i₊½ conserved metrics
  @inbounds for idx in m.iterators.cell.full
    i, j, k = idx.I .+ 0.5 # centroid index

    # get the conserved metrics at (i₊½, j, k)
    # @unpack ξ̂, η̂, ζ̂, J = conservative_metrics(m, (i + 1 / 2, j, k), t)

    # m.edge_metrics.i₊½.ξ̂.x[idx] = ξ̂.x
    # m.edge_metrics.i₊½.ξ̂.y[idx] = ξ̂.y
    # m.edge_metrics.i₊½.ξ̂.z[idx] = ξ̂.z
    # m.edge_metrics.i₊½.ξ̂.t[idx] = ξ̂.t

    # m.edge_metrics.i₊½.η̂.x[idx] = η̂.x
    # m.edge_metrics.i₊½.η̂.y[idx] = η̂.y
    # m.edge_metrics.i₊½.η̂.z[idx] = η̂.z
    # m.edge_metrics.i₊½.η̂.t[idx] = η̂.t

    # m.edge_metrics.i₊½.ζ̂.x[idx] = ζ̂.x
    # m.edge_metrics.i₊½.ζ̂.y[idx] = ζ̂.y
    # m.edge_metrics.i₊½.ζ̂.z[idx] = ζ̂.z
    # m.edge_metrics.i₊½.ζ̂.t[idx] = ζ̂.t

    m.edge_metrics.i₊½.J[idx] = jacobian(m, (i + 1 / 2, j, k), t)
  end

  # j₊½ conserved metrics
  @inbounds for idx in m.iterators.cell.full
    i, j, k = idx.I .+ 0.5 # centroid index

    # get the conserved metrics at (i, j₊½, k)
    # @unpack ξ̂, η̂, ζ̂, J = conservative_metrics(m, (i, j + 1 / 2, k), t)

    # m.edge_metrics.j₊½.ξ̂.x[idx] = ξ̂.x
    # m.edge_metrics.j₊½.ξ̂.y[idx] = ξ̂.y
    # m.edge_metrics.j₊½.ξ̂.z[idx] = ξ̂.z
    # m.edge_metrics.j₊½.ξ̂.t[idx] = ξ̂.t

    # m.edge_metrics.j₊½.η̂.x[idx] = η̂.x
    # m.edge_metrics.j₊½.η̂.y[idx] = η̂.y
    # m.edge_metrics.j₊½.η̂.z[idx] = η̂.z
    # m.edge_metrics.j₊½.η̂.t[idx] = η̂.t

    # m.edge_metrics.j₊½.ζ̂.x[idx] = ζ̂.x
    # m.edge_metrics.j₊½.ζ̂.y[idx] = ζ̂.y
    # m.edge_metrics.j₊½.ζ̂.z[idx] = ζ̂.z
    # m.edge_metrics.j₊½.ζ̂.t[idx] = ζ̂.t

    m.edge_metrics.j₊½.J[idx] = jacobian(m, (i, j + 1 / 2, k), t)
  end

  # k₊½ conserved metrics
  @inbounds for idx in m.iterators.cell.full
    i, j, k = idx.I .+ 0.5 # centroid index

    # get the conserved metrics at (i, j, k₊½)
    # @unpack ξ̂, η̂, ζ̂, J = conservative_metrics(m, (i, j, k + 1 / 2), t)

    # m.edge_metrics.k₊½.ξ̂.x[idx] = ξ̂.x
    # m.edge_metrics.k₊½.ξ̂.y[idx] = ξ̂.y
    # m.edge_metrics.k₊½.ξ̂.z[idx] = ξ̂.z
    # m.edge_metrics.k₊½.ξ̂.t[idx] = ξ̂.t

    # m.edge_metrics.k₊½.η̂.x[idx] = η̂.x
    # m.edge_metrics.k₊½.η̂.y[idx] = η̂.y
    # m.edge_metrics.k₊½.η̂.z[idx] = η̂.z
    # m.edge_metrics.k₊½.η̂.t[idx] = η̂.t

    # m.edge_metrics.k₊½.ζ̂.x[idx] = ζ̂.x
    # m.edge_metrics.k₊½.ζ̂.y[idx] = ζ̂.y
    # m.edge_metrics.k₊½.ζ̂.z[idx] = ζ̂.z
    # m.edge_metrics.k₊½.ζ̂.t[idx] = ζ̂.t

    m.edge_metrics.k₊½.J[idx] = jacobian(m, (i, j, k + 1 / 2), t)
  end

  return nothing
end

# ------------------------------------------------------------------
# Grid Metrics
# ------------------------------------------------------------------

@inline function metrics(m::AxisymmetricGrid1D, (i,)::NTuple{1,Real}, t::Real=0)
  idx_1d = (i, 1, 1)
  return metrics(m, idx_1d, t)
end

@inline function metrics(m::AxisymmetricGrid1D, (i, j, k)::NTuple{3,Real}, t::Real=0)
  _jacobian_matrix = checkeps(
    m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo)
  )
  J = det(_jacobian_matrix)

  inv_jacobian_matrix = inv(_jacobian_matrix)

  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ξz = inv_jacobian_matrix[1, 3]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]
  ηz = inv_jacobian_matrix[2, 3]
  ζx = inv_jacobian_matrix[3, 1]
  ζy = inv_jacobian_matrix[3, 2]
  ζz = inv_jacobian_matrix[3, 3]

  vx, vy, vz = grid_velocities(m, (i, j, k), t)
  ξt = -(vx * ξx + vy * ξy + vz * ξz)
  ηt = -(vx * ηx + vy * ηy + vz * ηz)
  ζt = -(vx * ζx + vy * ζy + vz * ζz)

  ξ = Metric3D(ξx, ξy, ξz, ξt)
  η = Metric3D(ηx, ηy, ηz, ηt)
  ζ = Metric3D(ζx, ζy, ζz, ζt)

  return (; ξ, η, ζ, J)
end

# ------------------------------------------------------------------
# Conservative Grid Metrics; e.g. ξ̂x = ξx * J
# ------------------------------------------------------------------

@inline function conservative_metrics(
  m::AxisymmetricGrid1D, (i,)::NTuple{1,Real}, t::Real=0
)
  idx_1d = (i, 1, 1)
  return conservative_metrics(m, idx_1d, t)
end

@inline function conservative_metrics(
  m::AxisymmetricGrid1D, (i, j, k)::NTuple{3,Real}, t::Real=0
)
  _jacobian_matrix = checkeps(
    m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo)
  )
  J = abs(det(_jacobian_matrix))

  inv_jacobian_matrix = inv(_jacobian_matrix)

  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ξz = inv_jacobian_matrix[1, 3]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]
  ηz = inv_jacobian_matrix[2, 3]
  ζx = inv_jacobian_matrix[3, 1]
  ζy = inv_jacobian_matrix[3, 2]
  ζz = inv_jacobian_matrix[3, 3]

  vx, vy, vz = grid_velocities(m, (i, j, k), t)
  ξt = -(vx * ξx + vy * ξy + vz * ξz)
  ηt = -(vx * ηx + vy * ηy + vz * ηz)
  ζt = -(vx * ζx + vy * ζy + vz * ζz)

  # xξ = _jacobian_matrix[1, 1]
  # yξ = _jacobian_matrix[2, 1]
  # zξ = _jacobian_matrix[3, 1]
  # xη = _jacobian_matrix[1, 2]
  # yη = _jacobian_matrix[2, 2]
  # zη = _jacobian_matrix[3, 2]
  # xζ = _jacobian_matrix[1, 3]
  # yζ = _jacobian_matrix[2, 3]
  # zζ = _jacobian_matrix[3, 3]

  ξ̂ = Metric3D(ξx * J, ξy * J, ξz * J, ξt * J)
  η̂ = Metric3D(ηx * J, ηy * J, ηz * J, ηt * J)
  ζ̂ = Metric3D(ζx * J, ζy * J, ζz * J, ζt * J)

  return (; ξ̂, η̂, ζ̂, J)
end

# ------------------------------------------------------------------
# Jacobian related functions
# ------------------------------------------------------------------
@inline function jacobian_matrix(m::AxisymmetricGrid1D, (i,)::NTuple{1,Real}, t::Real=0)
  idx_1d = (i, 1, 1)
  return jacobian_matrix(m, idx_1d, t)
end

@inline function jacobian_matrix(
  m::AxisymmetricGrid1D, (i, j, k)::NTuple{3,Real}, t::Real=0
)
  return m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo, t) |> checkeps
end

@inline function jacobian(m::AxisymmetricGrid1D, (i,)::NTuple{1,Real}, t::Real=0)
  idx_1d = (i, 1, 1)
  return jacobian(m, idx_1d, t)
end

@inline function jacobian(m::AxisymmetricGrid1D, (i, j, k)::NTuple{3,Real}, t::Real=0)
  return jacobian_matrix(m, (i, j, k), t) |> det |> abs
end

# ------------------------------------------------------------------
# Velocity Functions
# ------------------------------------------------------------------

@inline function grid_velocities(::AxisymmetricGrid1D, (i, j, k)::NTuple{3,Real}, t::Real=0)
  return zeros(3)
end
# @inline centroid_velocities(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)
# @inline node_velocities(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)

# ------------------------------------------------------------------
# Coordinate Functions
# ------------------------------------------------------------------

function _axisym_node_coordinates!(coordinates, coordinate_functions, domain, nhalo)

  # Populate the node coordinates
  for idx in domain
    i = idx - nhalo
    coordinates.r[idx] = coordinate_functions.r(i)
  end

  return nothing
end

function _axisym_centroid_coordinates!(centroids, coordinate_functions, domain, nhalo)

  # Populate the centroid coordinates
  for idx in domain
    i = idx - nhalo + 0.5
    centroids.r[idx] = coordinate_functions.r(i)
  end

  return nothing
end