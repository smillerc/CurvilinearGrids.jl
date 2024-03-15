
"""
CurvilinearGrid3D

# Fields
 - `x`: Node function; e.g., x(i,j,k)
 - `y`: Node function; e.g., y(i,j,k)
 - `z`: Node function; e.g., z(i,j,k)
 - `jacobian_matrix_func`: Function to compute the jacobian matrix, e.g., J(i,j,k)
 - `conserv_metric_func`: Function to compute the conservative metrics
 - `nhalo`: Number of halo cells for all dims
 - `nnodes`: Number of nodes/vertices
 - `limits`: Cell loop limits based on halo cells
"""
struct CurvilinearGrid3D{CO,CE,NV,EM,CM,DL,CI,DS,CF,JF} <: AbstractCurvilinearGrid
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::NTuple{3,Int}
  domain_limits::DL
  iterators::CI
  discretization_scheme::DS
  _coordinate_funcs::CF
  _jacobian_matrix_func::JF # jacobian_matrix(ξ,η,ζ)
end

function CurvilinearGrid3D(
  x::Function,
  y::Function,
  z::Function,
  (n_ξ, n_η, n_ζ),
  nhalo,
  discretization_scheme=:MEG6;
  T=Float64,
  backend=CPU(),
)
  dim = 3
  check_nargs(x, dim, :x)
  check_nargs(y, dim, :y)
  check_nargs(z, dim, :z)

  test_coord_func(x, dim, :x)
  test_coord_func(y, dim, :y)
  test_coord_func(z, dim, :z)

  coord(i, j, k) = @SVector [x(i, j, k), y(i, j, k), z(i, j, k)]
  function jacobian_matrix_func(i, j, k, t)
    return ForwardDiff.jacobian(x -> coord(x[1], x[2], x[3]), @SVector [i, j, k])
  end

  nnodes = (n_ξ, n_η, n_ζ)
  ncells = nnodes .- 1
  ni_cells = n_ξ - 1
  nj_cells = n_η - 1
  nk_cells = n_ζ - 1
  lo = nhalo + 1

  limits = (
    node=(ilo=lo, ihi=n_ξ + nhalo, jlo=lo, jhi=n_η + nhalo, klo=lo, khi=n_ζ + nhalo),
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

  domain_iterators = get_node_cell_iterators(nodeCI, cellCI, nhalo)

  if discretization_scheme === :MEG6 || discretization_scheme === :MonotoneExplicit6thOrder
    if nhalo < 4
      # if nhalo != 4
      error("`nhalo` must = 4 when using the MEG6 discretization scheme")
    end

    discr_scheme = MetricDiscretizationSchemes.MonotoneExplicit6thOrderDiscretization(
      domain_iterators.cell.full
    )

  else
    error("Unknown discretization scheme to compute the conserved metrics")
  end

  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  coordinate_funcs = (; x, y, z)
  centroids = StructArray((
    x=KernelAbstractions.zeros(backend, T, celldims),
    y=KernelAbstractions.zeros(backend, T, celldims),
    z=KernelAbstractions.zeros(backend, T, celldims),
  ))
  _centroid_coordinates!(centroids, coordinate_funcs, domain_iterators.cell.full, nhalo)

  coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
    z=KernelAbstractions.zeros(backend, T, nodedims),
  ))
  _node_coordinates!(coords, coordinate_funcs, domain_iterators.node.full, nhalo)

  node_velocities = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
    z=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  m = CurvilinearGrid3D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes,
    limits,
    domain_iterators,
    discr_scheme,
    coordinate_funcs,
    jacobian_matrix_func,
  )

  update_metrics!(m)
  check_for_invalid_metrics(m)
  return m
end

function update_metrics!(m::CurvilinearGrid3D, t::Real=0)
  # Update the metrics within the non-halo region, e.g., the domain
  domain = m.iterators.cell.domain

  # cell metrics
  @inbounds for idx in m.iterators.cell.full
    cell_idx = idx.I .+ 0.5
    # @unpack J, ξ, η, ζ, x, y, z = metrics(m, cell_idx, 0)
    @unpack J, ξ, η, ζ = metrics(m, cell_idx, t)

    m.cell_center_metrics.ξ.x₁[idx] = ξ.x₁
    m.cell_center_metrics.ξ.x₂[idx] = ξ.x₂
    m.cell_center_metrics.ξ.x₃[idx] = ξ.x₃
    m.cell_center_metrics.ξ.t[idx] = ξ.t

    m.cell_center_metrics.η.x₁[idx] = η.x₁
    m.cell_center_metrics.η.x₂[idx] = η.x₂
    m.cell_center_metrics.η.x₃[idx] = η.x₃
    m.cell_center_metrics.η.t[idx] = η.t

    m.cell_center_metrics.ζ.x₁[idx] = ζ.x₁
    m.cell_center_metrics.ζ.x₂[idx] = ζ.x₂
    m.cell_center_metrics.ζ.x₃[idx] = ζ.x₃
    m.cell_center_metrics.ζ.t[idx] = ζ.t

    # m.cell_center_inv_metrics.xξ[idx] = x.ξ
    # m.cell_center_inv_metrics.yξ[idx] = y.ξ
    # m.cell_center_inv_metrics.zξ[idx] = z.ξ
    # m.cell_center_inv_metrics.xη[idx] = x.η
    # m.cell_center_inv_metrics.yη[idx] = y.η
    # m.cell_center_inv_metrics.zη[idx] = z.η
    # m.cell_center_inv_metrics.xζ[idx] = x.ζ
    # m.cell_center_inv_metrics.yζ[idx] = y.ζ
    # m.cell_center_inv_metrics.zζ[idx] = z.ζ

    m.cell_center_metrics.J[idx] = J
  end

  MetricDiscretizationSchemes.update_metrics!(
    m.discretization_scheme,
    m.centroid_coordinates,
    m.cell_center_metrics,
    m.edge_metrics,
    domain,
  )

  return nothing
end

# ------------------------------------------------------------------
# Grid Metrics
# ------------------------------------------------------------------

# Get the grid metrics
@inline function metrics(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real=0)
  _jacobian_matrix = checkeps(
    m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo, t)
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
  m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real=0
)
  _jacobian_matrix = checkeps(
    m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo, t)
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

  return (; ξ̂, η̂, ζ̂)
end

# ------------------------------------------------------------------
# Jacobian related functions
# ------------------------------------------------------------------

function jacobian_matrix(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real=0)
  # return checkeps(m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo))
  return m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo, t)
end

function jacobian(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real=0)
  return det(jacobian_matrix(m, (i, j, k), t))
end

# ------------------------------------------------------------------
# Velocity Functions
# ------------------------------------------------------------------

@inline grid_velocities(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real=0) =
  (0.0, 0.0, 0.0)
# @inline centroid_velocities(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t) = (0.0, 0.0, 0.0)
# @inline node_velocities(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t) = (0.0, 0.0, 0.0)

# ------------------------------------------------------------------
# Coordinate Functions
# ------------------------------------------------------------------

function _node_coordinates!(
  coordinates::StructArray{T,3}, coordinate_functions, domain, nhalo
) where {T}

  # Populate the node coordinates
  @inbounds for idx in domain
    cell_idx = @. idx.I - nhalo
    coordinates.x[idx] = coordinate_functions.x(cell_idx...)
    coordinates.y[idx] = coordinate_functions.y(cell_idx...)
    coordinates.z[idx] = coordinate_functions.z(cell_idx...)
  end

  return nothing
end

function _centroid_coordinates!(
  centroids::StructArray{T,3}, coordinate_functions, domain, nhalo
) where {T}

  # Populate the centroid coordinates
  @inbounds for idx in domain
    cell_idx = @. idx.I - nhalo + 0.5
    centroids.x[idx] = coordinate_functions.x(cell_idx...)
    centroids.y[idx] = coordinate_functions.y(cell_idx...)
    centroids.z[idx] = coordinate_functions.z(cell_idx...)
  end

  return nothing
end
