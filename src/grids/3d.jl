
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
struct CurvilinearGrid3D{CO,CE,NV,EM,CM,DL,CI,TI,DS} <: AbstractCurvilinearGrid
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::NTuple{3,Int}
  domain_limits::DL
  iterators::CI
  tiles::TI
  discretization_scheme::DS
  onbc::@NamedTuple{ilo::Bool, ihi::Bool, jlo::Bool, jhi::Bool, klo::Bool, khi::Bool}
  is_static::Bool
  is_orthogonal::Bool
end

"""
    CurvilinearGrid3D(
  x::AbstractArray{T,3},
  y::AbstractArray{T,3},
  z::AbstractArray{T,3},
  nhalo::Int,
  discretization_scheme=:MEG6;
  backend=CPU(),
) where {T}

"""
function CurvilinearGrid3D(
  x::AbstractArray{T,3},
  y::AbstractArray{T,3},
  z::AbstractArray{T,3},
  nhalo::Int,
  discretization_scheme=:MEG6;
  backend=CPU(),
  on_bc=nothing,
  is_static=false,
  is_orthogonal=false,
  tiles=nothing,
) where {T}

  #
  @assert size(x) == size(y) == size(y)

  nnodes = size(x)
  ni, nj, nk = nnodes

  ncells = nnodes .- 1
  ni_cells = ni - 1
  nj_cells = nj - 1
  nk_cells = nk - 1
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

  domain_iterators = get_node_cell_iterators(nodeCI, cellCI, nhalo)

  if discretization_scheme === :MEG6 || discretization_scheme === :MonotoneExplicit6thOrder
    if nhalo < 4
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

  coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
    z=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  @views begin
    copy!(coords.x[domain_iterators.node.domain], x)
    copy!(coords.y[domain_iterators.node.domain], y)
    copy!(coords.z[domain_iterators.node.domain], z)
  end

  centroids = StructArray((
    x=KernelAbstractions.zeros(backend, T, celldims),
    y=KernelAbstractions.zeros(backend, T, celldims),
    z=KernelAbstractions.zeros(backend, T, celldims),
  ))
  # _centroid_coordinates!(centroids, coords, domain_iterators.cell.domain)

  node_velocities = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
    z=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  if isnothing(on_bc)
    _on_bc = (ilo=true, ihi=true, jlo=true, jhi=true, klo=true, khi=true)
  else
    _on_bc = on_bc
  end

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
    tiles,
    discr_scheme,
    _on_bc,
    is_static,
    is_orthogonal,
  )

  update!(m; force=true)
  return m
end

"""Update metrics after grid coordinates change"""
function update!(mesh::CurvilinearGrid3D; force=false)
  if !mesh.is_static || force
    _centroid_coordinates!(
      mesh.centroid_coordinates, mesh.node_coordinates, mesh.iterators.cell.domain
    )
    update_metrics!(mesh)
    _check_valid_metrics(mesh)
  else
    @warn("Attempting to update grid metrics when grid.is_static = true!")
  end
  return nothing
end

function update_metrics!(m::CurvilinearGrid3D, t::Real=0)
  # Update the metrics within the non-halo region, e.g., the domain
  domain = m.iterators.cell.domain

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
# @inline function metrics(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real=0)
#   _jacobian_matrix = checkeps(
#     m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo, t)
#   )
#   J = det(_jacobian_matrix)

#   inv_jacobian_matrix = inv(_jacobian_matrix)

#   ξx = inv_jacobian_matrix[1, 1]
#   ξy = inv_jacobian_matrix[1, 2]
#   ξz = inv_jacobian_matrix[1, 3]
#   ηx = inv_jacobian_matrix[2, 1]
#   ηy = inv_jacobian_matrix[2, 2]
#   ηz = inv_jacobian_matrix[2, 3]
#   ζx = inv_jacobian_matrix[3, 1]
#   ζy = inv_jacobian_matrix[3, 2]
#   ζz = inv_jacobian_matrix[3, 3]

#   vx, vy, vz = grid_velocities(m, (i, j, k), t)
#   ξt = -(vx * ξx + vy * ξy + vz * ξz)
#   ηt = -(vx * ηx + vy * ηy + vz * ηz)
#   ζt = -(vx * ζx + vy * ζy + vz * ζz)

#   ξ = Metric3D(ξx, ξy, ξz, ξt)
#   η = Metric3D(ηx, ηy, ηz, ηt)
#   ζ = Metric3D(ζx, ζy, ζz, ζt)

#   return (; ξ, η, ζ, J)
# end

# ------------------------------------------------------------------
# Conservative Grid Metrics; e.g. ξ̂x = ξx * J
# ------------------------------------------------------------------

# @inline function conservative_metrics(
#   m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real=0
# )
#   _jacobian_matrix = checkeps(
#     m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo, t)
#   )
#   J = det(_jacobian_matrix)

#   inv_jacobian_matrix = inv(_jacobian_matrix)

#   ξx = inv_jacobian_matrix[1, 1]
#   ξy = inv_jacobian_matrix[1, 2]
#   ξz = inv_jacobian_matrix[1, 3]
#   ηx = inv_jacobian_matrix[2, 1]
#   ηy = inv_jacobian_matrix[2, 2]
#   ηz = inv_jacobian_matrix[2, 3]
#   ζx = inv_jacobian_matrix[3, 1]
#   ζy = inv_jacobian_matrix[3, 2]
#   ζz = inv_jacobian_matrix[3, 3]

#   vx, vy, vz = grid_velocities(m, (i, j, k), t)
#   ξt = -(vx * ξx + vy * ξy + vz * ξz)
#   ηt = -(vx * ηx + vy * ηy + vz * ηz)
#   ζt = -(vx * ζx + vy * ζy + vz * ζz)

#   # xξ = _jacobian_matrix[1, 1]
#   # yξ = _jacobian_matrix[2, 1]
#   # zξ = _jacobian_matrix[3, 1]
#   # xη = _jacobian_matrix[1, 2]
#   # yη = _jacobian_matrix[2, 2]
#   # zη = _jacobian_matrix[3, 2]
#   # xζ = _jacobian_matrix[1, 3]
#   # yζ = _jacobian_matrix[2, 3]
#   # zζ = _jacobian_matrix[3, 3]

#   ξ̂ = Metric3D(ξx * J, ξy * J, ξz * J, ξt * J)
#   η̂ = Metric3D(ηx * J, ηy * J, ηz * J, ηt * J)
#   ζ̂ = Metric3D(ζx * J, ζy * J, ζz * J, ζt * J)

#   return (; ξ̂, η̂, ζ̂, J)
# end

# ------------------------------------------------------------------
# Jacobian related functions
# ------------------------------------------------------------------

# function jacobian_matrix(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real=0)
#   # return checkeps(m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo))
#   return m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo, t)
# end

# function jacobian(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real=0)
#   return det(jacobian_matrix(m, (i, j, k), t))
# end

# ------------------------------------------------------------------
# Velocity Functions
# ------------------------------------------------------------------

@inline grid_velocities(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real=0) =
  (0.0, 0.0, 0.0)
# @inline centroid_velocities(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t) = (0.0, 0.0, 0.0)
# @inline node_velocities(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t) = (0.0, 0.0, 0.0)

# ------------------------------------------------------------------
# Private Functions
# ------------------------------------------------------------------

function _centroid_coordinates!(
  centroids::StructArray{T,3}, coords::StructArray{T,3}, domain
) where {T}

  # Populate the centroid coordinates
  x = coords.x
  y = coords.y
  z = coords.z

  @batch for idx in domain
    i, j, k = idx.I
    #! format: off
    centroids.x[idx] = 0.125(
      x[i, j, k    ] + x[i + 1, j, k    ] + x[i + 1, j + 1, k    ] + x[i, j + 1, k    ] +
      x[i, j, k + 1] + x[i + 1, j, k + 1] + x[i + 1, j + 1, k + 1] + x[i, j + 1, k + 1]
    )
    
    centroids.y[idx] = 0.125(
      y[i, j, k    ] + y[i + 1, j, k    ] + y[i + 1, j + 1, k    ] + y[i, j + 1, k    ] +
      y[i, j, k + 1] + y[i + 1, j, k + 1] + y[i + 1, j + 1, k + 1] + y[i, j + 1, k + 1]
    )
    
    centroids.z[idx] = 0.125(
      z[i, j, k    ] + z[i + 1, j, k    ] + z[i + 1, j + 1, k    ] + z[i, j + 1, k    ] +
      z[i, j, k + 1] + z[i + 1, j, k + 1] + z[i + 1, j + 1, k + 1] + z[i, j + 1, k + 1]
    )
    #! format: on

  end

  return nothing
end

function _check_valid_metrics(mesh::CurvilinearGrid3D)
  domain = mesh.iterators.cell.domain
  # i₊½_domain = expand(domain, 1, -1)
  # j₊½_domain = expand(domain, 2, -1)

  @views begin
    centroid_metrics_valid =
      all(isfinite.(mesh.cell_center_metrics.J[domain])) &&
      all(mesh.cell_center_metrics.J[domain] .> 0)
    #     all(isfinite.(mesh.cell_center_metrics.ξ.x₁[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.ξ.x₂[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.η.x₁[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.η.x₂[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.x₁.ξ[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.x₁.η[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.x₂.ξ[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.x₂.η[domain]))

    # edge_metrics_valid =
    # all(isfinite.(mesh.edge_metrics.i₊½.J[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₁[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₂[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.t[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.η̂.x₁[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.η̂.x₂[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.η̂.t[i₊½_domain])) &&
    # all(isfinite.(mesh.edge_metrics.j₊½.J[j₊½_domain])) # &&
    #     all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.x₁[j₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.x₂[j₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.t[j₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.j₊½.η̂.x₁[j₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.j₊½.η̂.x₂[j₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.j₊½.η̂.t[j₊½_domain]))
  end

  # if !edge_metrics_valid
  #   error("Invalid edge metrics found")
  # end

  if !centroid_metrics_valid
    error("Invalid centroid metrics found")
  end

  return nothing
end
