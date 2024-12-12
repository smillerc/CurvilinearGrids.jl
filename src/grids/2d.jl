
struct CurvilinearGrid2D{CO,CE,NV,EM,CM,DL,CI,TI,DS} <: AbstractCurvilinearGrid2D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::NTuple{2,Int}
  domain_limits::DL
  iterators::CI
  tiles::TI
  discretization_scheme::DS
  onbc::@NamedTuple{ilo::Bool, ihi::Bool, jlo::Bool, jhi::Bool}
  is_static::Bool
  is_orthogonal::Bool
end

"""
AxisymmetricGrid2D
"""
struct AxisymmetricGrid2D{CO,CE,NV,EM,CM,DL,CI,TI,DS} <: AbstractCurvilinearGrid2D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::NTuple{2,Int}
  domain_limits::DL
  iterators::CI
  tiles::TI
  discretization_scheme::DS
  snap_to_axis::Bool
  rotational_axis::Symbol
  onbc::@NamedTuple{ilo::Bool, ihi::Bool, jlo::Bool, jhi::Bool}
  is_static::Bool
  is_orthogonal::Bool
end

"""
    CurvilinearGrid2D(x::AbstractArray{T,2}, y::AbstractArray{T,2}, nhalo::Int, discretization_scheme=:MEG6; backend=CPU()) where {T}

Create a 2d grid with `x` and `y` coordinates. The input coordinates do not include halo / ghost data since the geometry is undefined in these regions.
The `nhalo` argument defines the number of halo cells along each dimension.
"""
function CurvilinearGrid2D(
  x::AbstractArray{T,2},
  y::AbstractArray{T,2},
  nhalo::Int;
  discretization_scheme=:MEG6,
  backend=CPU(),
  on_bc=nothing,
  is_static=false,
  is_orthogonal=false,
  tiles=nothing,
  make_uniform=false,
) where {T}

  #
  @assert size(x) == size(y)

  nnodes = size(x)
  ni, nj = nnodes

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

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  @views begin
    copy!(coords.x[domain_iterators.node.domain], x)
    copy!(coords.y[domain_iterators.node.domain], y)
  end

  centroids = StructArray((
    x=KernelAbstractions.zeros(backend, T, celldims),
    y=KernelAbstractions.zeros(backend, T, celldims),
  ))

  _centroid_coordinates!(centroids, coords, domain_iterators.cell.domain)

  node_velocities = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  # if discretization_scheme === :MEG6 || discretization_scheme === :MonotoneExplicit6thOrder
  #   if nhalo < 4
  #     error("`nhalo` must = 4 when using the MEG6 discretization scheme")
  #   end

  discr_scheme = MetricDiscretizationSchemes.MontoneExplicitGradientScheme6thOrder(;
    use_cache=true, celldims=size(domain_iterators.cell.full), backend=CPU(), T=Float64
  )

  # else
  #   error("Unknown discretization scheme to compute the conserved metrics")
  # end

  if isnothing(on_bc)
    _on_bc = (ilo=true, ihi=true, jlo=true, jhi=true)
  else
    _on_bc = on_bc
  end

  m = CurvilinearGrid2D(
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

  if make_uniform
    first_idx = first(m.iterators.cell.domain)

    fill!(m.cell_center_metrics.J, m.cell_center_metrics.J[first_idx])

    fill!(m.cell_center_metrics.x₁.ξ, m.cell_center_metrics.x₁.ξ[first_idx])
    fill!(m.cell_center_metrics.x₂.ξ, m.cell_center_metrics.x₂.ξ[first_idx])
    fill!(m.cell_center_metrics.x₁.η, m.cell_center_metrics.x₁.η[first_idx])
    fill!(m.cell_center_metrics.x₂.η, m.cell_center_metrics.x₂.η[first_idx])

    fill!(m.cell_center_metrics.ξ.x₁, m.cell_center_metrics.ξ.x₁[first_idx])
    fill!(m.cell_center_metrics.ξ.x₂, m.cell_center_metrics.ξ.x₂[first_idx])
    fill!(m.cell_center_metrics.η.x₁, m.cell_center_metrics.η.x₁[first_idx])
    fill!(m.cell_center_metrics.η.x₂, m.cell_center_metrics.η.x₂[first_idx])

    fill!(m.edge_metrics.i₊½.J, m.edge_metrics.i₊½.J[first_idx])

    fill!(m.edge_metrics.i₊½.ξ̂.x₁, m.edge_metrics.i₊½.ξ̂.x₁[first_idx])
    fill!(m.edge_metrics.i₊½.ξ̂.x₂, m.edge_metrics.i₊½.ξ̂.x₂[first_idx])
    fill!(m.edge_metrics.i₊½.ξ̂.t, m.edge_metrics.i₊½.ξ̂.t[first_idx])
    fill!(m.edge_metrics.i₊½.η̂.x₁, m.edge_metrics.i₊½.η̂.x₁[first_idx])
    fill!(m.edge_metrics.i₊½.η̂.x₂, m.edge_metrics.i₊½.η̂.x₂[first_idx])
    fill!(m.edge_metrics.i₊½.η̂.t, m.edge_metrics.i₊½.η̂.t[first_idx])

    fill!(m.edge_metrics.j₊½.J, m.edge_metrics.j₊½.J[first_idx])

    fill!(m.edge_metrics.j₊½.ξ̂.x₁, m.edge_metrics.j₊½.ξ̂.x₁[first_idx])
    fill!(m.edge_metrics.j₊½.ξ̂.x₂, m.edge_metrics.j₊½.ξ̂.x₂[first_idx])
    fill!(m.edge_metrics.j₊½.ξ̂.t, m.edge_metrics.j₊½.ξ̂.t[first_idx])
    fill!(m.edge_metrics.j₊½.η̂.x₁, m.edge_metrics.j₊½.η̂.x₁[first_idx])
    fill!(m.edge_metrics.j₊½.η̂.x₂, m.edge_metrics.j₊½.η̂.x₂[first_idx])
    fill!(m.edge_metrics.j₊½.η̂.t, m.edge_metrics.j₊½.η̂.t[first_idx])
  end

  return m
end

"""
    AxisymmetricGrid2D(x::AbstractMatrix{T}, y::AbstractMatrix{T},  nhalo::Int,  snap_to_axis::Bool,  rotational_axis::Symbol;  T=Float64, backend=CPU())

Create an axisymmetric 2d grid with `x` and `y` coordinates, with the symmetry axis `rotational_axis = :x` or `rotational_axis = :y`. Enforce coordinates to stay on the axis via `snap_to_axis=True`.
The input coordinates do not include halo / ghost data since the geometry is undefined in these regions. The `nhalo` argument defines the number of halo cells along each dimension.
"""
function AxisymmetricGrid2D(
  x::AbstractMatrix{T},
  y::AbstractMatrix{T},
  nhalo::Int,
  snap_to_axis::Bool,
  rotational_axis::Symbol;
  discretization_scheme=:MEG6,
  backend=CPU(),
  is_static=false,
  on_bc=nothing,
  is_orthogonal=false,
  tiles=nothing,
  make_uniform=false,
) where {T}

  #
  @assert size(x) == size(y)
  @assert rotational_axis === :x || rotational_axis === :y

  nnodes = size(x)
  ni, nj = nnodes

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

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  @views begin
    copy!(coords.x[domain_iterators.node.domain], x)
    copy!(coords.y[domain_iterators.node.domain], y)
  end

  centroids = StructArray((
    x=KernelAbstractions.zeros(backend, T, celldims),
    y=KernelAbstractions.zeros(backend, T, celldims),
  ))

  # _centroid_coordinates!(centroids, coords, domain_iterators.cell.domain)

  node_velocities = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  # if discretization_scheme === :MEG6 || discretization_scheme === :MonotoneExplicit6thOrder
  #   if nhalo < 4
  #     error("`nhalo` must = 4 when using the MEG6 discretization scheme")
  #   end

  discr_scheme = MetricDiscretizationSchemes.MonotoneExplicit6thOrderDiscretization(
    domain_iterators.cell.full
  )

  # else
  #   error("Unknown discretization scheme to compute the conserved metrics")
  # end

  if isnothing(on_bc)
    _on_bc = (ilo=true, ihi=true, jlo=true, jhi=true)
  else
    _on_bc = on_bc
  end

  m = AxisymmetricGrid2D(
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
    snap_to_axis,
    rotational_axis,
    _on_bc,
    is_static,
    is_orthogonal,
  )

  update!(m)

  if make_uniform
    first_idx = first(m.iterators.cell.domain)

    fill!(m.cell_center_metrics.J, m.cell_center_metrics.J[first_idx])

    fill!(m.cell_center_metrics.x₁.ξ, m.cell_center_metrics.x₁.ξ[first_idx])
    fill!(m.cell_center_metrics.x₂.ξ, m.cell_center_metrics.x₂.ξ[first_idx])
    fill!(m.cell_center_metrics.x₁.η, m.cell_center_metrics.x₁.η[first_idx])
    fill!(m.cell_center_metrics.x₂.η, m.cell_center_metrics.x₂.η[first_idx])

    fill!(m.cell_center_metrics.ξ.x₁, m.cell_center_metrics.ξ.x₁[first_idx])
    fill!(m.cell_center_metrics.ξ.x₂, m.cell_center_metrics.ξ.x₂[first_idx])
    fill!(m.cell_center_metrics.η.x₁, m.cell_center_metrics.η.x₁[first_idx])
    fill!(m.cell_center_metrics.η.x₂, m.cell_center_metrics.η.x₂[first_idx])

    fill!(m.edge_metrics.i₊½.J, m.edge_metrics.i₊½.J[first_idx])

    fill!(m.edge_metrics.i₊½.ξ̂.x₁, m.edge_metrics.i₊½.ξ̂.x₁[first_idx])
    fill!(m.edge_metrics.i₊½.ξ̂.x₂, m.edge_metrics.i₊½.ξ̂.x₂[first_idx])
    fill!(m.edge_metrics.i₊½.ξ̂.t, m.edge_metrics.i₊½.ξ̂.t[first_idx])
    fill!(m.edge_metrics.i₊½.η̂.x₁, m.edge_metrics.i₊½.η̂.x₁[first_idx])
    fill!(m.edge_metrics.i₊½.η̂.x₂, m.edge_metrics.i₊½.η̂.x₂[first_idx])
    fill!(m.edge_metrics.i₊½.η̂.t, m.edge_metrics.i₊½.η̂.t[first_idx])

    fill!(m.edge_metrics.j₊½.J, m.edge_metrics.j₊½.J[first_idx])

    fill!(m.edge_metrics.j₊½.ξ̂.x₁, m.edge_metrics.j₊½.ξ̂.x₁[first_idx])
    fill!(m.edge_metrics.j₊½.ξ̂.x₂, m.edge_metrics.j₊½.ξ̂.x₂[first_idx])
    fill!(m.edge_metrics.j₊½.ξ̂.t, m.edge_metrics.j₊½.ξ̂.t[first_idx])
    fill!(m.edge_metrics.j₊½.η̂.x₁, m.edge_metrics.j₊½.η̂.x₁[first_idx])
    fill!(m.edge_metrics.j₊½.η̂.x₂, m.edge_metrics.j₊½.η̂.x₂[first_idx])
    fill!(m.edge_metrics.j₊½.η̂.t, m.edge_metrics.j₊½.η̂.t[first_idx])
  end

  return m
end

"""Update metrics after grid coordinates change"""
function update!(mesh::CurvilinearGrid2D; force=false)
  if !mesh.is_static || force
    _centroid_coordinates!(
      mesh.centroid_coordinates, mesh.node_coordinates, mesh.iterators.cell.domain
    )
    update_metrics!(mesh)
    _check_valid_metrics(mesh)
  else
    error("Attempting to update grid metrics when grid.is_static = true!")
  end
  return nothing
end

"""Update metrics after grid coordinates change"""
function update!(mesh::AxisymmetricGrid2D)
  _centroid_coordinates!(
    mesh.centroid_coordinates, mesh.node_coordinates, mesh.iterators.cell.domain
  )
  update_metrics!(mesh)

  if mesh.snap_to_axis
    _snap_nodes_to_axis(mesh)
  else
    _check_nodes_along_axis(mesh)
  end

  _check_valid_metrics(mesh)
  return nothing
end

function update_metrics!(mesh::AbstractCurvilinearGrid2D, t::Real=0)
  # Update the metrics within the non-halo region, e.g., the domain
  domain = mesh.iterators.cell.domain

  MetricDiscretizationSchemes.update_metrics!(
    mesh.discretization_scheme,
    mesh.centroid_coordinates,
    mesh.cell_center_metrics,
    mesh.edge_metrics,
    domain,
  )

  return nothing
end

"""
    jacobian_matrix(mesh::CurvilinearGrid2D, idx)

The cell-centroid Jacobian matrix (the forward transformation: ∂x/∂ξ, ∂y/∂ξ, ... ). Use `inv(jacobian_matrix(mesh, idx))` to
get the inverse transformation (∂ξ/∂x, ∂ξ/∂y, ...)
"""
function jacobian_matrix(mesh::CurvilinearGrid2D, idx)
  xξ = mesh.cell_center_metrics.x₁.ξ
  yξ = mesh.cell_center_metrics.x₂.ξ
  xη = mesh.cell_center_metrics.x₁.η
  yη = mesh.cell_center_metrics.x₂.η

  return @SMatrix [
    xξ[idx] xη[idx]
    yξ[idx] yη[idx]
  ]
end

"""
    jacobian(mesh::CurvilinearGrid2D, idx)

The cell-centroid Jacobian (determinant of the Jacobian matrix)
"""
jacobian(mesh::CurvilinearGrid2D, idx) = det(jacobian_matrix(mesh, idx))

# ------------------------------------------------------------------
# Velocity Functions
# ------------------------------------------------------------------

@inline grid_velocities(::AbstractCurvilinearGrid2D, (i, j)::NTuple{2,Real}, t::Real=0) =
  (0.0, 0.0)
# @inline centroid_velocities(mesh::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)
# @inline node_velocities(mesh::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)

# ------------------------------------------------------------------
# Private Functions
# ------------------------------------------------------------------

function _centroid_coordinates!(
  centroids::StructArray{T,2}, coords::StructArray{T,2}, domain
) where {T}
  x = coords.x
  y = coords.y
  # Populate the centroid coordinates
  for idx in domain
    i, j = idx.I
    centroids.x[idx] = 0.25(x[i, j] + x[i + 1, j] + x[i + 1, j + 1] + x[i, j + 1])
    centroids.y[idx] = 0.25(y[i, j] + y[i + 1, j] + y[i + 1, j + 1] + y[i, j + 1])
  end

  return nothing
end

function _check_valid_metrics(mesh::AbstractCurvilinearGrid2D)
  domain = mesh.iterators.cell.domain
  i₊½_domain = expand(domain, 1, -1)
  j₊½_domain = expand(domain, 2, -1)

  @views begin
    centroid_metrics_valid =
      all(isfinite.(mesh.cell_center_metrics.forward.J[domain])) &&
      all(mesh.cell_center_metrics.forward.J[domain] .> 0)
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

  # if !centroid_metrics_valid
  #   error("Invalid centroid metrics found")
  # end

  return nothing
end

function _check_nodes_along_axis(mesh::AxisymmetricGrid2D)
  domain = mesh.iterators.node.domain

  @views begin
    if mesh.rotational_axis === :x
      axis_domain = domain[:, 1]

      if any(!iszero(mesh.node_coordinates.y[axis_domain]))
        error(
          "Nodes not aligned to axis of symmetry (`snap_to_axis = true`). Set snap_to_axis = false to disable checks",
        )
      end
    else # mesh.rotational_axis === :y
      axis_domain = domain[1, :]

      if any(!iszero(mesh.node_coordinates.x[axis_domain]))
        error(
          "Nodes not aligned to axis of symmetry (`snap_to_axis = true`). Set snap_to_axis = false to disable checks",
        )
      end
    end
  end
end

function _snap_nodes_to_axis(mesh::AxisymmetricGrid2D)
  domain = mesh.iterators.node.domain

  @views begin
    if mesh.rotational_axis === :x
      axis_domain = domain[:, 1]
      mesh.node_coordinates.y[axis_domain] .= 0

    else # mesh.rotational_axis === :y
      axis_domain = domain[1, :]
      mesh.node_coordinates.x[axis_domain] .= 0
    end
  end
end
