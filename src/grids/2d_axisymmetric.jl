
struct RZAxisymmetricGrid2D{CO,CE,NV,EM,CM,DL,CI,CF,JF} <: AbstractCurvilinearGrid
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::NTuple{2,Int}
  domain_limits::DL
  iterators::CI
  _coordinate_funcs::CF
  _jacobian_matrix_func::JF
end

function RZAxisymmetricGrid2D(
  r::Function, z::Function, (ni, nj), nhalo; T=Float64, backend=CPU()
)

  # Ensure that the r and z functions are set up properly, i.e.,
  # they are defined as r(i,j) = ... and z(i,j) = ...
  dim = 2
  check_nargs(r, dim, :r)
  check_nargs(z, dim, :z)
  test_coord_func(r, dim, :r)
  test_coord_func(z, dim, :z)

  # Make a full 3d grid with only 1 cell in θ.
  # This is cheap and very useful for certain applications

  θ1 = 2#π # leave the π off so we can use the more accurate cospi/sinpi functions
  θ(j) = θ1 * (j - 1)

  R3d(i, j, k) = r(i, k) * cospi(θ(j))
  Θ3d(i, j, k) = r(i, k) * sinpi(θ(j))
  Z3d(i, j, k) = z(i, k)

  RΘZ(i, j, k) = @SVector [R3d(i, j, k), Θ3d(i, j, k), Z3d(i, j, k)]
  function jacobian_matrix_func(i, j, k)
    return ForwardDiff.jacobian(x -> RΘZ(x[1], x[2], x[3]), @SVector [i, j, k])
  end

  # jacobian_matrix_func = _setup_jacobian_func(x, y)
  nnodes = (ni, nj)
  ncells = nnodes .- 1
  ni_cells, nj_cells = ncells
  lo = nhalo + 1
  limits = (
    node=(ilo=lo, ihi=ni + nhalo, jlo=lo, jhi=nj + nhalo),
    cell=(ilo=lo, ihi=ni_cells + nhalo, jlo=lo, jhi=nj_cells + nhalo),
  )

  nodeCI = CartesianIndices(nnodes .+ 2nhalo)
  cellCI = CartesianIndices(ncells .+ 2nhalo)

  domain_iterators = get_node_cell_iterators(nodeCI, cellCI, nhalo)
  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  cell_center_metrics = (
    J=KernelAbstractions.zeros(backend, T, celldims),
    ξ=StructArray((
      x=KernelAbstractions.zeros(backend, T, celldims),
      y=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    η=StructArray((
      x=KernelAbstractions.zeros(backend, T, celldims),
      y=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
  )

  edge_metrics = (
    i₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ̂=StructArray((
        x=KernelAbstractions.zeros(backend, T, celldims),
        y=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η̂=StructArray((
        x=KernelAbstractions.zeros(backend, T, celldims),
        y=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
    ),
    j₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ̂=StructArray((
        x=KernelAbstractions.zeros(backend, T, celldims),
        y=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η̂=StructArray((
        x=KernelAbstractions.zeros(backend, T, celldims),
        y=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
    ),
  )

  coordinate_funcs = (; x=r, y=z)
  centroids = StructArray((
    x=KernelAbstractions.zeros(backend, T, celldims),
    y=KernelAbstractions.zeros(backend, T, celldims),
  ))
  _centroid_coordinates!(centroids, coordinate_funcs, domain_iterators.cell.full, nhalo)

  coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
  ))
  _node_coordinates!(coords, coordinate_funcs, domain_iterators.node.full, nhalo)

  node_velocities = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  m = RZAxisymmetricGrid2D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes,
    limits,
    domain_iterators,
    coordinate_funcs,
    jacobian_matrix_func,
  )

  update_metrics!(m)
  # check_for_invalid_metrics(m)
  return m
end

function update_metrics!(m::RZAxisymmetricGrid2D, t=0)

  # cell metrics
  @inbounds for idx in m.iterators.cell.full
    cell_idx = idx.I .+ 0.5
    # @unpack J, ξ, η, x, y = metrics(m, cell_idx, t)
    @unpack J, ξ, η = metrics(m, cell_idx, t)

    m.cell_center_metrics.ξ.x[idx] = ξ.x
    m.cell_center_metrics.ξ.y[idx] = ξ.y
    m.cell_center_metrics.ξ.t[idx] = ξ.t
    m.cell_center_metrics.η.x[idx] = η.x
    m.cell_center_metrics.η.y[idx] = η.y
    m.cell_center_metrics.η.t[idx] = η.t

    # m.cell_center_inv_metrics.xξ[idx] = x.ξ
    # m.cell_center_inv_metrics.yξ[idx] = y.ξ
    # m.cell_center_inv_metrics.xη[idx] = x.η
    # m.cell_center_inv_metrics.yη[idx] = y.η

    m.cell_center_metrics.J[idx] = J
  end

  # i₊½ conserved metrics
  @inbounds for idx in m.iterators.cell.full
    i, j = idx.I .+ 0.5 # centroid index

    # get the conserved metrics at (i₊½, j)
    @unpack ξ̂, η̂, J = conservative_metrics(m, (i + 1 / 2, j), t)

    m.edge_metrics.i₊½.ξ̂.x[idx] = ξ̂.x
    m.edge_metrics.i₊½.ξ̂.y[idx] = ξ̂.y
    m.edge_metrics.i₊½.ξ̂.t[idx] = ξ̂.t
    m.edge_metrics.i₊½.η̂.x[idx] = η̂.x
    m.edge_metrics.i₊½.η̂.y[idx] = η̂.y
    m.edge_metrics.i₊½.η̂.t[idx] = η̂.t
    m.edge_metrics.i₊½.J[idx] = J
  end

  # j₊½ conserved metrics
  @inbounds for idx in m.iterators.cell.full
    i, j = idx.I .+ 0.5 # centroid index

    # get the conserved metrics at (i, j₊½)
    @unpack ξ̂, η̂, J = conservative_metrics(m, (i, j + 1 / 2), t)

    m.edge_metrics.j₊½.ξ̂.x[idx] = ξ̂.x
    m.edge_metrics.j₊½.ξ̂.y[idx] = ξ̂.y
    m.edge_metrics.j₊½.ξ̂.t[idx] = ξ̂.t
    m.edge_metrics.j₊½.η̂.x[idx] = η̂.x
    m.edge_metrics.j₊½.η̂.y[idx] = η̂.y
    m.edge_metrics.j₊½.η̂.t[idx] = η̂.t
    m.edge_metrics.j₊½.J[idx] = J
  end

  return nothing
end

# ------------------------------------------------------------------
# Grid Metrics
# ------------------------------------------------------------------
metrics(m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real}) = metrics(m, (i, j), 0.0)

@inline function metrics(m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real}, t::Real)

  # Get the full 3d jacobian matrix. The 2nd coordinate doesn't matter
  # since it's symmetric about θ
  _jacobian_matrix = checkeps(m._jacobian_matrix_func(i - m.nhalo, 1, j - m.nhalo))
  inv_jacobian_matrix = inv(_jacobian_matrix)

  # Only extract the ∂(r,z) terms
  ξr = inv_jacobian_matrix[1, 1]
  ξz = inv_jacobian_matrix[1, 3]
  ηr = inv_jacobian_matrix[3, 1]
  ηz = inv_jacobian_matrix[3, 3]

  # In this scenario, J is the volume of the node/cell at (i,j),
  # and it includes the revolution term. This is important!
  J = det(_jacobian_matrix)

  vr, vz = grid_velocities(m, (i, j), t)
  ξt = -(vr * ξr + vz * ξz)
  ηt = -(vr * ηr + vz * ηz)

  ξ = Metric2D(ξr, ξz, ξt)
  η = Metric2D(ηr, ηz, ηt)

  return (; ξ, η, J)
end

@inline function planar_metrics(m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real}, t::Real)

  # Get the full 3d jacobian matrix. The 2nd coordinate doesn't matter
  # since it's symmetric about θ
  _jacobian_matrix = checkeps(m._jacobian_matrix_func(i - m.nhalo, 1, j - m.nhalo))
  inv_jacobian_matrix = inv(_jacobian_matrix)

  # Only extract the ∂(r,z) terms
  ξr = inv_jacobian_matrix[1, 1]
  ξz = inv_jacobian_matrix[1, 3]
  ηr = inv_jacobian_matrix[3, 1]
  ηz = inv_jacobian_matrix[3, 3]
  rξ = _jacobian_matrix[1, 1]
  zξ = _jacobian_matrix[1, 3]
  rη = _jacobian_matrix[3, 1]
  zη = _jacobian_matrix[3, 3]

  # In this scenario, J is the AREA of the node/cell at (i,j),
  # and it DOES NOT include the revolution term. This is important!
  J = rξ * zη - zξ * rη

  vr, vz = grid_velocities(m, (i, j), t)
  ξt = -(vr * ξr + vz * ξz)
  ηt = -(vr * ηr + vz * ηz)

  ξ = Metric2D(ξr, ξz, ξt)
  η = Metric2D(ηr, ηz, ηt)

  return (; ξ, η, J)
end

@inline function metrics(m::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real}, t::Real)
  _jacobian_matrix = checkeps(
    m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo)
  )
  inv_jacobian_matrix = inv(_jacobian_matrix)

  # In this scenario, J is the true VOLUME of the node/cell at (i,j,k).
  # This is an important distinction!
  # One may ask how a node has a volume, but then you may be accused of
  # being a pain in the ass... 😆 Just think of it as the the volume
  # occupied by the polyhedral element defined by the centroids
  # of the surrounding cells.
  J = det(_jacobian_matrix)

  ξx₁ = inv_jacobian_matrix[1, 1]
  ξx₂ = inv_jacobian_matrix[1, 2]
  ξx₃ = inv_jacobian_matrix[1, 3]
  ηx₁ = inv_jacobian_matrix[2, 1]
  ηx₂ = inv_jacobian_matrix[2, 2]
  ηx₃ = inv_jacobian_matrix[2, 3]
  ζx₁ = inv_jacobian_matrix[3, 1]
  ζx₂ = inv_jacobian_matrix[3, 2]
  ζx₃ = inv_jacobian_matrix[3, 3]

  vx₁, vx₂, vx₃ = grid_velocities(m, (i, j, k), t)
  ξt = -(vx₁ * ξx₁ + vx₂ * ξx₂ + vx₃ * ξx₃)
  ηt = -(vx₁ * ηx₁ + vx₂ * ηx₂ + vx₃ * ηx₃)
  ζt = -(vx₁ * ζx₁ + vx₂ * ζx₂ + vx₃ * ζx₃)

  ξ = Metric3D(ξx₁, ξx₂, ξx₃, ξt)
  η = Metric3D(ηx₁, ηx₂, ηx₃, ηt)
  ζ = Metric3D(ζx₁, ζx₂, ζx₃, ζt)

  return (; ξ, η, ζ, J)
end

# ------------------------------------------------------------------
# Conservative Grid Metrics; e.g. ξ̂x = ξx * J
# ------------------------------------------------------------------

@inline function conservative_metrics(
  m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real}, t::Real
)

  # Get the full 3d jacobian matrix. The 2nd coordinate doesn't matter
  # since it's symmetric about θ
  _jacobian_matrix = checkeps(m._jacobian_matrix_func(i - m.nhalo, 1, j - m.nhalo))
  inv_jacobian_matrix = inv(_jacobian_matrix)

  # Only extract the ∂(r,z) terms
  ξr = inv_jacobian_matrix[1, 1]
  ξz = inv_jacobian_matrix[1, 3]
  ηr = inv_jacobian_matrix[3, 1]
  ηz = inv_jacobian_matrix[3, 3]

  # In this scenario, J is the volume of the node/cell at (i,j),
  # and it includes the revolution term. This is important!
  J = det(_jacobian_matrix)

  vr, vz = grid_velocities(m, (i, j), t)
  ξt = -(vr * ξr + vz * ξz)
  ηt = -(vr * ηr + vz * ηz)

  ξ̂ = Metric2D(ξr * J, ξz * J, ξt * J)
  η̂ = Metric2D(ηr * J, ηz * J, ηt * J)

  return (; ξ̂, η̂, J)
end

# ------------------------------------------------------------------
# Jacobian related functions
# ------------------------------------------------------------------
function jacobian_matrix(m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real})
  _jacobian_matrix = jacobian_matrix(m, (i - m.nhalo, 1, j - m.nhalo))
  T = eltype(_jacobian_matrix)

  # extract only the 2d portion
  _jacobian_matrix_2d = SMatrix{2,2}(
    _jacobian_matrix[1, 1],
    _jacobian_matrix[3, 1],
    _jacobian_matrix[1, 3],
    _jacobian_matrix[3, 3],
  )

  return checkeps(_jacobian_matrix_2d)
end

function jacobian_matrix(m::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real})
  return checkeps(m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo))
end

function jacobian(m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real})
  _jacobian_matrix = checkeps(jacobian_matrix(m, (i - m.nhalo, j - m.nhalo)))
  return det(_jacobian_matrix)
end

function jacobian(m::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real})
  _jacobian_matrix = checkeps(jacobian_matrix(m, (i - m.nhalo, j - m.nhalo, k - m.nhalo)))
  return det(_jacobian_matrix)
end

# ------------------------------------------------------------------
# Velocity Functions
# ------------------------------------------------------------------

@inline grid_velocities(m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)
@inline grid_velocities(m::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real}, t) =
  (0.0, 0.0, 0.0)
# @inline centroid_velocities(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)
# @inline node_velocities(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)

area(m::RZAxisymmetricGrid2D, (i, k)::NTuple{2,Real}) = jacobian(m, (i, k))
function area(::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real})
  return error(
    """
    You're trying to get the area of a 3d index in a RZAxisymmetricGrid2D,
    which doesn't make physical sense! Use one of the following:
    1. `area(grid, (i,j))`, which will get you the non-rotate area of the node/cell
    2. `volume(grid, (i,j,k))` to get the true "rotate" volume of the node/cell
    """
  )
end

volume(m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real}) = jacobian(m, (i, j))
volume(m::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real}) = jacobian(m, (i, k))
