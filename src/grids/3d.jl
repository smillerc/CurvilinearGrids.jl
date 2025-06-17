
struct CurvilinearGrid3D{CO,CE,NV,EM,CM,DL,CI,DS} <: AbstractCurvilinearGrid3D
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
  onbc::@NamedTuple{ilo::Bool, ihi::Bool, jlo::Bool, jhi::Bool, klo::Bool, khi::Bool}
  is_static::Bool
  is_orthogonal::Bool
  discretization_scheme_name::Symbol
end

"""
    CurvilinearGrid3D(x, y, z, nhalo::Int; backend=CPU(), discretization_scheme=:MEG6, on_bc=nothing, is_static=false, is_orthogonal=false, make_uniform=false)

Construct a curvilinear grid in 3D using 3D arrays of x/y/z coordinates.
"""
function CurvilinearGrid3D(
  x::AbstractArray{T,3},
  y::AbstractArray{T,3},
  z::AbstractArray{T,3},
  discretization_scheme::Symbol;
  backend=CPU(),
  on_bc=nothing,
  is_static=false,
  is_orthogonal=false,
  make_uniform=false,
  init_metrics=true,
  kwargs...,
) where {T}

  #
  use_symmetric_conservative_metric_scheme = false

  scheme_name = Symbol(uppercase("$discretization_scheme"))
  if scheme_name === :MEG6 ||
    discretization_scheme == :MontoneExplicitGradientScheme6thOrder
    MetricDiscretizationScheme = MontoneExplicitGradientScheme6thOrder
    nhalo = 5
  elseif scheme_name === :MEG6_SYMMETRIC ||
    discretization_scheme == :MontoneExplicitGradientScheme6thOrder
    MetricDiscretizationScheme = MontoneExplicitGradientScheme6thOrder
    nhalo = 5
    use_symmetric_conservative_metric_scheme = true
  else
    error("Only MontoneExplicitGradientScheme6thOrder or MEG6 is supported for now")
  end

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

  discr_scheme = MetricDiscretizationScheme(;
    use_cache=true,
    celldims=size(domain_iterators.cell.full),
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=use_symmetric_conservative_metric_scheme,
  )

  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  # if init_metrics
  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)
  # else
  #   cell_center_metrics = nothing
  #   edge_metrics = nothing
  # end

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
    discr_scheme,
    _on_bc,
    is_static,
    is_orthogonal,
    scheme_name,
  )

  if init_metrics
    update!(m; force=true)
  end

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

function jacobian_matrix(mesh::CurvilinearGrid3D, (i, j, k))
  xξ = mesh.cell_center_metrics.x₁.ξ
  yξ = mesh.cell_center_metrics.x₂.ξ
  zξ = mesh.cell_center_metrics.x₂.ξ
  xη = mesh.cell_center_metrics.x₁.η
  yη = mesh.cell_center_metrics.x₂.η
  zη = mesh.cell_center_metrics.x₂.η
  xζ = mesh.cell_center_metrics.x₁.η
  yζ = mesh.cell_center_metrics.x₂.η
  zζ = mesh.cell_center_metrics.x₂.η

  return @SMatrix [
    xξ[i, j, k] xη[i, j, k] xζ[i, j, k]
    yξ[i, j, k] yη[i, j, k] yζ[i, j, k]
    zξ[i, j, k] zη[i, j, k] zζ[i, j, k]
  ]
end

"""
    forward_cell_metrics(mesh::AbstractCurvilinearGrid3D, idx)

Get the forward cell metrics for a given cell index `idx`
"""
function forward_cell_metrics(mesh::AbstractCurvilinearGrid3D, idx)
  (;
    x=mesh.cell_center_metrics.x₁[idx],
    y=mesh.cell_center_metrics.x₂[idx],
    z=mesh.cell_center_metrics.x₃[idx],
  )
end

"""
    inverse_cell_metrics(mesh::AbstractCurvilinearGrid3D, idx)

Get the inverse cell metrics for a given cell index `idx`
"""
function inverse_cell_metrics(mesh::AbstractCurvilinearGrid3D, idx)
  (;
    ξ=mesh.cell_center_metrics.ξ[idx],
    η=mesh.cell_center_metrics.η[idx],
    ζ=mesh.cell_center_metrics.ζ[idx],
    ξ̂=mesh.cell_center_metrics.ξ̂[idx],
    η̂=mesh.cell_center_metrics.η̂[idx],
    ζ̂=mesh.cell_center_metrics.ζ̂[idx],
  )
end

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
  i₊½_domain = expand(domain, 1, -1)
  j₊½_domain = expand(domain, 2, -1)
  k₊½_domain = expand(domain, 3, -1)

  @views begin
    centroid_metrics_valid =
      all(isfinite.(mesh.cell_center_metrics.J[domain])) &&
      all(mesh.cell_center_metrics.J[domain] .> 0) &&
      all(isfinite.(mesh.cell_center_metrics.ξ.x₁[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.ξ.x₂[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.ξ.x₃[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.η.x₁[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.η.x₂[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.η.x₃[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.ζ.x₁[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.ζ.x₂[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.ζ.x₃[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₁.ξ[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₁.η[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₁.ζ[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₂.ξ[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₂.η[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₂.ζ[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₃.ξ[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₃.η[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₃.ζ[domain]))

    edge_metrics_valid =
      all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₁[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₂[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₃[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.t[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.η̂.x₁[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.η̂.x₂[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.η̂.x₃[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.η̂.t[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ζ̂.x₁[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ζ̂.x₂[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ζ̂.x₃[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ζ̂.t[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.x₁[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.x₂[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.x₃[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.t[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.η̂.x₁[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.η̂.x₂[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.η̂.x₃[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.η̂.t[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ζ̂.x₁[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ζ̂.x₂[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ζ̂.x₃[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ζ̂.t[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ξ̂.x₁[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ξ̂.x₂[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ξ̂.x₃[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ξ̂.t[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.η̂.x₁[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.η̂.x₂[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.η̂.x₃[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.η̂.t[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ζ̂.x₁[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ζ̂.x₂[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ζ̂.x₃[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ζ̂.t[k₊½_domain]))
  end

  if !edge_metrics_valid
    error("Invalid edge metrics found")
  end

  if !centroid_metrics_valid
    error("Invalid centroid metrics found")
  end

  return nothing
end
