function PartitionedCurvilinearGrid(
  x::AbstractArray{T,1}, nhalo::Int, partition_fraction::Int, rank::Int; backend=CPU()
) where {T}
  @assert nhalo > 0

  if rank == 0
    error("The rank input for this call must be 1-based, not 0-based!")
  end

  # x, y are coordinate arrays
  cell_domain = CartesianIndices(size(x) .- 1)

  if any(partition_fraction < 1)
    error("All entries in `partition_fraction` must be >= 1")
  end

  partition_fraction = max(partition_fraction, 1) # ensure no zeros or negative numbers
  on_bc, tiled_node_limits, node_subdomain, _ = get_subdomain_limits(
    cell_domain, partition_fraction, rank
  )

  # Create a partitioned CurvilinearGrid2D with the local subdomain coordinates
  x_local = @view x[node_subdomain]

  return CurvilinearGrid1D(
    x_local, nhalo; backend=backend, on_bc=on_bc, tiles=tiled_node_limits
  )
end

function PartitionedCurvilinearGrid(
  x::AbstractArray{T,2},
  y::AbstractArray{T,2},
  nhalo::Int,
  partition_fraction::NTuple{2,Int},
  rank::Int;
  backend=CPU(),
) where {T}
  @assert nhalo > 0

  if rank == 0
    error("The rank input for this call must be 1-based, not 0-based!")
  end

  if size(x) !== size(y)
    error("The x and y arrays must be the same size!")
  end

  # x, y are coordinate arrays
  cell_domain = CartesianIndices(size(x) .- 1)

  if any(partition_fraction .< 1)
    error("All entries in `partition_fraction` must be >= 1")
  end

  partition_fraction = max.(partition_fraction, 1) # ensure no zeros or negative numbers
  on_bc, tiled_node_limits, node_subdomain, _ = get_subdomain_limits(
    cell_domain, partition_fraction, rank
  )

  # Create a partitioned CurvilinearGrid2D with the local subdomain coordinates
  x_local = @view x[node_subdomain]
  y_local = @view y[node_subdomain]

  return CurvilinearGrid2D(
    x_local, y_local, nhalo; backend=backend, on_bc=on_bc, tiles=tiled_node_limits
  )
end

function PartitionedCurvilinearGrid(
  x::AbstractArray{T,3},
  y::AbstractArray{T,3},
  z::AbstractArray{T,3},
  nhalo::Int,
  partition_fraction::NTuple{3,Int},
  rank::Int;
  backend=CPU(),
) where {T}
  @assert nhalo > 0

  if rank == 0
    error("The rank input for this call must be 1-based, not 0-based!")
  end

  if size(x) !== size(y)
    error("The x and y arrays must be the same size!")
  end

  cell_domain = CartesianIndices(size(x) .- 1)

  if any(partition_fraction .< 1)
    error("All entries in `partition_fraction` must be >= 1")
  end

  partition_fraction = max.(partition_fraction, 1) # ensure no zeros or negative numbers
  on_bc, tiled_node_limits, node_subdomain, _ = get_subdomain_limits(
    cell_domain, partition_fraction, rank
  )

  # Create a partitioned CurvilinearGrid2D with the local subdomain coordinates
  x_local = @view x[node_subdomain]
  y_local = @view y[node_subdomain]
  z_local = @view z[node_subdomain]

  return CurvilinearGrid3D(
    x_local, y_local, z_local, nhalo; backend=backend, on_bc=on_bc, tiles=tiled_node_limits
  )
end
