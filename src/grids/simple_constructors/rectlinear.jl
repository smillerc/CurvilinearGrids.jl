function RectlinearGrid(
  x0,
  x1,
  ncells,
  nhalo::Int,
  backend=CPU(),
  T=Float64,
  is_static=true,
  make_uniform=false,
  tile_layout=nothing,
  rank::Int=-1,
)
  ni = ncells + 1
  x = collect(T, range(x0, x1; length=ni))
  return CurvilinearGrid1D(x, nhalo; backend=backend)
end

function RectlinearGrid(
  x::AbstractVector{T},
  nhalo::Int,
  backend=CPU(),
  is_static=true,
  make_uniform=false,
  tile_layout=nothing,
  rank::Int=-1,
) where {T}

  #
  ni = length(x)

  if ni < 2
    error("The x vector must have more than 2 points")
  end

  if !isnothing(tile_layout)
    if rank == -1
      error(
        "Tile layout is provided, but rank is invalid; make sure to specify the current MPI rank",
      )
    end

    if rank == 0
      error(
        "Rank is 0, (MPI is zero-based), but needs to be one-based, do rank+1 for this call"
      )
    end

    partition_fraction = max.(tile_layout, 1) # ensure no zeros or negative numbers
    on_bc, tiled_node_limits, node_subdomain, cell_subdomain = get_subdomain_limits(
      cell_domain, partition_fraction, rank
    )

    if length(tiled_node_limits) != prod(partition_fraction)
      error("Unable to partition the mesh to the desired tile layout")
    end

    x_local = zeros(T, size(node_subdomain))

    @inbounds for (gidx, lidx) in
                  zip(node_subdomain, CartesianIndices(size(node_subdomain)))
      i, = lidx.I
      gi, = gidx.I
      x_local[i, j] = x[gi]
    end

    return CurvilinearGrid1D(
      x_local,
      nhalo;
      backend=backend,
      on_bc=on_bc,
      tiles=tiled_node_limits,
      make_uniform=make_uniform,
    )

  else
    return CurvilinearGrid1D(
      x,
      nhalo;
      backend=backend,
      is_orthogonal=true,
      is_static=is_static,
      make_uniform=make_uniform,
    )
  end
end

function RectlinearCylindricalGrid(
  r0,
  r1,
  ncells,
  nhalo::Int;
  snap_to_axis=true,
  backend=CPU(),
  T=Float64,
  is_static=true,
  make_uniform=false,
  tile_layout=nothing,
  rank::Int=-1,
)
  ni = ncells + 1
  r = collect(T, range(r0, r1; length=ni))
  return CylindricalGrid1D(r, nhalo, snap_to_axis; backend=backend)
end

function RectlinearSphericalGrid(
  r0,
  r1,
  ncells,
  nhalo::Int;
  snap_to_axis=true,
  backend=CPU(),
  T=Float64.is_static = true,
  make_uniform=false,
  tile_layout=nothing,
  rank::Int=-1,
)
  ni = ncells + 1
  r = collect(T, range(r0, r1; length=ni))
  return SphericalGrid1D(r, nhalo, snap_to_axis; backend=backend)
end

function RectlinearGrid(
  (x0, y0),
  (x1, y1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  nhalo::Int,
  backend=CPU(),
  T=Float64;
  is_static=true,
  make_uniform=false,
  tile_layout=nothing,
  rank::Int=-1,
)
  if ni_cells < 2 || nj_cells < 2
    error(
      "The number of cells specified must be > 1, given cell dims are $((ni_cells, nj_cells))",
    )
  end

  return RectlinearGrid(
    range(x0, x1; length=ni_cells + 1),
    range(y0, y1; length=nj_cells + 1),
    nhalo,
    backend;
    is_static=is_static,
    make_uniform=make_uniform,
    tile_layout=tile_layout,
    rank=rank,
  )
end

function RectlinearGrid(
  x::AbstractVector{T},
  y::AbstractVector{T},
  nhalo::Int,
  backend=CPU();
  is_static=true,
  make_uniform=false,
  tile_layout=nothing,
  rank::Int=-1,
) where {T}
  ni = length(x)
  nj = length(y)

  if ni < 2
    error("The x vector must have more than 2 points")
  end

  if nj < 2
    error("The y vector must have more than 2 points")
  end

  cell_domain = CartesianIndices((ni - 1, nj - 1))

  if !all(diff(x) .> 0)
    error("Invalid x vector, spacing between vertices must be > 0 everywhere")
  end

  if !all(diff(y) .> 0)
    error("Invalid y vector, spacing between vertices must be > 0 everywhere")
  end

  if !isnothing(tile_layout)
    if rank == -1
      error(
        "Tile layout is provided, but rank is invalid; make sure to specify the current MPI rank",
      )
    end

    if rank == 0
      error(
        "Rank is 0, (MPI is zero-based), but needs to be one-based, do rank+1 for this call"
      )
    end

    partition_fraction = max.(tile_layout, 1) # ensure no zeros or negative numbers
    on_bc, tiled_node_limits, node_subdomain, cell_subdomain = get_subdomain_limits(
      cell_domain, partition_fraction, rank
    )

    if length(tiled_node_limits) != prod(partition_fraction)
      error("Unable to partition the mesh to the desired tile layout")
    end

    x_local = zeros(T, size(node_subdomain))
    y_local = zeros(T, size(node_subdomain))

    @inbounds for (gidx, lidx) in
                  zip(node_subdomain, CartesianIndices(size(node_subdomain)))
      i, j = lidx.I
      gi, gj = gidx.I
      x_local[i, j] = x[gi]
      y_local[i, j] = y[gj]
    end

    return CurvilinearGrid2D(
      x_local,
      y_local,
      nhalo;
      backend=backend,
      on_bc=on_bc,
      tiles=tiled_node_limits,
      make_uniform=make_uniform,
    )

  else
    x2d = zeros(T, ni, nj)
    y2d = zeros(T, ni, nj)

    @inbounds for j in 1:nj
      for i in 1:ni
        x2d[i, j] = x[i]
        y2d[i, j] = y[j]
      end
    end

    return CurvilinearGrid2D(
      x2d,
      y2d,
      nhalo;
      backend=backend,
      is_orthogonal=true,
      is_static=is_static,
      make_uniform=make_uniform,
    )
  end
end

function AxisymmetricRectlinearGrid(
  (x0, y0),
  (x1, y1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  nhalo::Int;
  snap_to_axis::Bool,
  rotational_axis::Symbol,
  backend=CPU(),
  T=Float64,
  is_static=true,
)
  if ni_cells < 2 || nj_cells < 2
    error("The number of cells specified must be > 2")
  end

  if x0 ≈ x1
    error("The endpoints x0 and x1 are the same")
  end
  if y0 ≈ y1
    error("The endpoints y0 and y1 are the same")
  end

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

  return AxisymmetricGrid2D(
    x,
    y,
    nhalo,
    snap_to_axis,
    rotational_axis;
    backend=backend,
    is_orthogonal=true,
    is_static=is_static,
  )
end

function RectlinearGrid(
  (x0, y0, z0),
  (x1, y1, z1),
  (ni_cells, nj_cells, nk_cells)::NTuple{3,Int},
  nhalo::Int;
  backend=CPU(),
  T=Float64,
  is_static=true,
  make_uniform=false,
  tile_layout=nothing,
  rank::Int=-1,
)
  if ni_cells < 2 || nj_cells < 2 || nk_cells < 2
    error(
      "The number of cells specified must be > 1, given cell dims are $((ni_cells, nj_cells, nk_cells))",
    )
  end

  return RectlinearGrid(
    range(x0, x1; length=ni_cells + 1),
    range(y0, y1; length=nj_cells + 1),
    range(z0, z1; length=nk_cells + 1),
    nhalo;
    is_static=is_static,
    backend=backend,
    make_uniform=make_uniform,
    tile_layout=tile_layout,
    rank=rank,
  )
end

function RectlinearGrid(
  x::AbstractVector{T},
  y::AbstractVector{T},
  z::AbstractVector{T},
  nhalo::Int;
  backend=CPU(),
  is_static=true,
  make_uniform=false,
  tile_layout=nothing,
  rank::Int=-1,
) where {T}

  #
  ni = length(x)
  nj = length(y)
  nk = length(z)

  if ni < 2
    error("The x vector must have more than 2 points")
  end

  if nj < 2
    error("The y vector must have more than 2 points")
  end

  if nk < 2
    error("The z vector must have more than 2 points")
  end

  if !all(diff(x) .> 0)
    error("Invalid x vector, spacing between vertices must be > 0 everywhere")
  end

  if !all(diff(y) .> 0)
    error("Invalid y vector, spacing between vertices must be > 0 everywhere")
  end

  if !all(diff(z) .> 0)
    error("Invalid z vector, spacing between vertices must be > 0 everywhere")
  end

  if !isnothing(tile_layout)
    if rank == -1
      error(
        "Tile layout is provided, but rank is invalid; make sure to specify the current MPI rank",
      )
    end

    if rank == 0
      error(
        "Rank is 0, (MPI is zero-based), but needs to be one-based, do rank+1 for this call"
      )
    end

    partition_fraction = max.(tile_layout, 1) # ensure no zeros or negative numbers
    on_bc, tiled_node_limits, node_subdomain, cell_subdomain = get_subdomain_limits(
      cell_domain, partition_fraction, rank
    )

    if length(tiled_node_limits) != prod(partition_fraction)
      error("Unable to partition the mesh to the desired tile layout")
    end

    x_local = zeros(T, size(node_subdomain))
    y_local = zeros(T, size(node_subdomain))
    z_local = zeros(T, size(node_subdomain))

    @inbounds for (gidx, lidx) in
                  zip(node_subdomain, CartesianIndices(size(node_subdomain)))
      i, j, k = lidx.I
      gi, gj, gk = gidx.I
      x_local[i, j, k] = x[gi]
      y_local[i, j, k] = y[gj]
      z_local[i, j, k] = z[gk]
    end

    return CurvilinearGrid2D(
      x_local,
      y_local,
      z_local,
      nhalo;
      backend=backend,
      on_bc=on_bc,
      tiles=tiled_node_limits,
      make_uniform=make_uniform,
    )

  else
    x3d = zeros(T, ni, nj, nk)
    y3d = zeros(T, ni, nj, nk)
    z3d = zeros(T, ni, nj, nk)

    @inbounds for k in 1:nk
      for j in 1:nj
        for i in 1:ni
          x3d[i, j, k] = x[i]
          y3d[i, j, k] = y[j]
          z3d[i, j, k] = z[k]
        end
      end
    end

    return CurvilinearGrid3D(
      x3d,
      y3d,
      z3d,
      nhalo;
      backend=backend,
      is_orthogonal=true,
      is_static=is_static,
      make_uniform=make_uniform,
    )
  end
end
