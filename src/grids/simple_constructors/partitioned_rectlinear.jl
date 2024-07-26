
"""
    PartitionedRectlinearGrid(
  x0, x1, ncells, nhalo::Int, partition_fraction::Int, rank::Int; backend=CPU(), T=Float64
)

TBW
"""
function PartitionedRectlinearGrid(
  x0, x1, ncells, nhalo::Int, partition_fraction::Int, rank::Int; backend=CPU(), T=Float64
)
  @assert !(x0 ≈ x1) "The endpoints x0 and x1 are the same"

  ni = ncells + 1

  cell_domain = CartesianIndices((n_cells))

  if any(partition_fraction < 1)
    error("All entries in `partition_fraction` must be >= 1")
  end

  partition_fraction = max(partition_fraction, 1) # ensure no zeros or negative numbers
  on_bc, tiled_node_limits, node_subdomain, cell_subdomain = get_subdomain_limits(
    cell_domain, partition_fraction, rank
  )

  x = zeros(T, size(node_subdomain))

  x1d = range(x0, x1; length=ni)

  @inbounds for (gidx, lidx) in zip(node_subdomain, CartesianIndices(size(node_subdomain)))
    i, = lidx.I
    gi, = gidx.I
    x[i] = x1d[gi]
  end

  return CurvilinearGrid1D(x, nhalo; backend=backend, on_bc=on_bc, tiles=tiled_node_limits)
end

"""
    PartitionedRectlinearGrid(
  (x0, y0),
  (x1, y1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  nhalo::Int,
  partition_fraction::NTuple{2,Int},
  rank::Int;
  backend=CPU(),
  T=Float64,
)

TBW
"""
function PartitionedRectlinearGrid(
  (x0, y0),
  (x1, y1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  nhalo::Int,
  partition_fraction::NTuple{2,Int},
  rank::Int;
  backend=CPU(),
  T=Float64,
)
  @assert !(x0 ≈ x1) "The endpoints x0 and x1 are the same"
  @assert !(y0 ≈ y1) "The endpoints y0 and y1 are the same"

  ni = ni_cells + 1
  nj = nj_cells + 1

  cell_domain = CartesianIndices((ni_cells, nj_cells))

  if any(partition_fraction .< 1)
    error("All entries in `partition_fraction` must be >= 1")
  end

  partition_fraction = max.(partition_fraction, 1) # ensure no zeros or negative numbers
  on_bc, tiled_node_limits, node_subdomain, cell_subdomain = get_subdomain_limits(
    cell_domain, partition_fraction, rank
  )

  x = zeros(T, size(node_subdomain))
  y = zeros(T, size(node_subdomain))

  x1d = range(x0, x1; length=ni)
  y1d = range(y0, y1; length=nj)

  @inbounds for (gidx, lidx) in zip(node_subdomain, CartesianIndices(size(node_subdomain)))
    i, j = lidx.I
    gi, gj = gidx.I
    x[i, j] = x1d[gi]
    y[i, j] = y1d[gj]
  end

  return CurvilinearGrid2D(
    x, y, nhalo; backend=backend, on_bc=on_bc, tiles=tiled_node_limits
  )
end

"""
    PartitionedRectlinearGrid(
  x::AbstractVector{T},
  y::AbstractVector{T},
  nhalo::Int,
  partition_fraction::NTuple{2,Int},
  rank::Int;
  backend=CPU(),
) where {T}

TBW
"""
function PartitionedRectlinearGrid(
  x::AbstractVector{T},
  y::AbstractVector{T},
  nhalo::Int,
  partition_fraction::NTuple{2,Int},
  rank::Int;
  backend=CPU(),
) where {T}
  ni = length(x)
  nj = length(y)

  @assert ni >= 2 "The x vector must have more than 2 points"
  @assert nj >= 2 "The y vector must have more than 2 points"

  cell_domain = CartesianIndices((ni - 1, nj - 1))

  if any(partition_fraction .< 1)
    error("All entries in `partition_fraction` must be >= 1")
  end

  partition_fraction = max.(partition_fraction, 1) # ensure no zeros or negative numbers
  on_bc, tiled_node_limits, node_subdomain, cell_subdomain = get_subdomain_limits(
    cell_domain, partition_fraction, rank
  )

  x2d = zeros(T, size(node_subdomain))
  y2d = zeros(T, size(node_subdomain))

  @inbounds for (gidx, lidx) in zip(node_subdomain, CartesianIndices(size(node_subdomain)))
    i, j = lidx.I
    gi, gj = gidx.I
    x2d[i, j] = x[gi]
    y2d[i, j] = y[gj]
  end

  return CurvilinearGrid2D(
    x2d, y2d, nhalo; backend=backend, on_bc=on_bc, tiles=tiled_node_limits
  )
end

"""
    PartitionedAxisymmetricRectlinearGrid(
  (x0, y0),
  (x1, y1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  nhalo::Int,
  snap_to_axis::Bool,
  rotational_axis::Symbol,
  partition_fraction::NTuple{2,Int},
  rank::Int;
  backend=CPU(),
  T=Float64,
)

TBW
"""
function PartitionedAxisymmetricRectlinearGrid(
  (x0, y0),
  (x1, y1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  nhalo::Int,
  snap_to_axis::Bool,
  rotational_axis::Symbol,
  partition_fraction::NTuple{2,Int},
  rank::Int;
  backend=CPU(),
  T=Float64,
)
  @assert !(x0 ≈ x1) "The endpoints x0 and x1 are the same"
  @assert !(y0 ≈ y1) "The endpoints y0 and y1 are the same"

  ni = ni_cells + 1
  nj = nj_cells + 1

  cell_domain = CartesianIndices((ni_cells, nj_cells))
  if any(partition_fraction .< 1)
    error("All entries in `partition_fraction` must be >= 1")
  end

  partition_fraction = max.(partition_fraction, 1) # ensure no zeros or negative numbers
  on_bc, tiled_node_limits, node_subdomain, cell_subdomain = get_subdomain_limits(
    cell_domain, partition_fraction, rank
  )

  x = zeros(T, size(node_subdomain))
  y = zeros(T, size(node_subdomain))

  x1d = range(x0, x1; length=ni)
  y1d = range(y0, y1; length=nj)

  @inbounds for (gidx, lidx) in zip(node_subdomain, CartesianIndices(size(node_subdomain)))
    i, j = lidx.I
    gi, gj = gidx.I
    x[i, j] = x1d[gi]
    y[i, j] = y1d[gj]
  end

  return AxisymmetricGrid2D(
    x,
    y,
    nhalo,
    snap_to_axis,
    rotational_axis;
    backend=backend,
    on_bc=on_bc,
    tiles=tiled_node_limits,
  )
end

function PartitionedRectlinearGrid(
  (x0, y0, z0),
  (x1, y1, z1),
  (ni_cells, nj_cells, nk_cells)::NTuple{3,Int},
  nhalo::Int,
  partition_fraction::NTuple{3,Int},
  rank::Int;
  backend=CPU(),
  T=Float64,
)
  ni = ni_cells + 1
  nj = nj_cells + 1
  nk = nk_cells + 1

  cell_domain = CartesianIndices((ni_cells, nj_cells))
  if any(partition_fraction .< 1)
    error("All entries in `partition_fraction` must be >= 1")
  end

  partition_fraction = max.(partition_fraction, 1) # ensure no zeros or negative numbers
  on_bc, tiled_node_limits, node_subdomain, cell_subdomain = get_subdomain_limits(
    cell_domain, partition_fraction, rank
  )

  x = zeros(T, size(node_subdomain))
  y = zeros(T, size(node_subdomain))
  z = zeros(T, size(node_subdomain))

  # node positions (non-halo)
  x1d = range(x0, x1; length=ni)
  y1d = range(y0, y1; length=nj)
  z1d = range(z0, z1; length=nk)

  @inbounds for (gidx, lidx) in zip(node_subdomain, CartesianIndices(size(node_subdomain)))
    i, j, k = lidx.I
    gi, gj, gk = gidx.I
    x[i, j, k] = x1d[gi]
    y[i, j, k] = y1d[gj]
    z[i, j, k] = z1d[gk]
  end

  return CurvilinearGrid3D(
    x, y, z, nhalo; backend=backend, on_bc=on_bc, tiles=tiled_node_limits
  )
end

function PartitionedRectlinearGrid(
  x::AbstractVector{T},
  y::AbstractVector{T},
  z::AbstractVector{T},
  nhalo::Int,
  partition_fraction::NTuple{3,Int},
  rank::Int;
  backend=CPU(),
) where {T}
  ni = length(x)
  nj = length(y)
  nk = length(z)

  @assert ni >= 2 "The x vector must have more than 2 points"
  @assert nj >= 2 "The y vector must have more than 2 points"
  @assert nk >= 2 "The z vector must have more than 2 points"

  cell_domain = CartesianIndices((ni - 1, nj - 1, nk - 1))

  if any(partition_fraction .< 1)
    error("All entries in `partition_fraction` must be >= 1")
  end

  partition_fraction = max.(partition_fraction, 1) # ensure no zeros or negative numbers
  on_bc, tiled_node_limits, node_subdomain, cell_subdomain = get_subdomain_limits(
    cell_domain, partition_fraction, rank
  )

  @inbounds for (gidx, lidx) in zip(node_subdomain, CartesianIndices(size(node_subdomain)))
    i, j = lidx.I
    gi, gj = gidx.I
    x3d[i, j, k] = x[gi]
    y3d[i, j, k] = y[gj]
    z3d[i, j, k] = z[gk]
  end

  return CurvilinearGrid3D(
    x3d, y3d, z3d, nhalo; backend=backend, on_bc=on_bc, tiles=tiled_node_limits
  )
end

function subdomain_on_bc(subdomain, domain::CartesianIndices{1})
  return (
    ilo=first(subdomain.indices[1]) == first(domain.indices[1]),
    ihi=last(subdomain.indices[1]) - 1 == last(domain.indices[1]),
  )
end

function subdomain_on_bc(subdomain, domain::CartesianIndices{2})
  return (
    ilo=first(subdomain.indices[1]) == first(domain.indices[1]),
    ihi=last(subdomain.indices[1]) - 1 == last(domain.indices[1]),
    jlo=first(subdomain.indices[2]) == first(domain.indices[2]),
    jhi=last(subdomain.indices[2]) - 1 == last(domain.indices[2]),
  )
end

function subdomain_on_bc(subdomain, domain::CartesianIndices{3})
  return (
    ilo=first(subdomain.indices[1]) == first(domain.indices[1]),
    ihi=last(subdomain.indices[1]) - 1 == last(domain.indices[1]),
    jlo=first(subdomain.indices[2]) == first(domain.indices[2]),
    jhi=last(subdomain.indices[2]) - 1 == last(domain.indices[2]),
    klo=first(subdomain.indices[3]) == first(domain.indices[3]),
    khi=last(subdomain.indices[3]) - 1 == last(domain.indices[3]),
  )
end