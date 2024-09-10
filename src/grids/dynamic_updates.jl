
"""
    _update_coords(mesh, update_func!, params)

Update the coordinates of the mesh (centroid and node), compute grid velocities,
and update the temporal metric terms

"""
function _update_coords(mesh, update_func!, params)

  #
  Δt = params.Δt

  # call the update function -> updates the node coordinates in-place
  update_func!(mesh.node_coordinates, params)

  x⃗cⁿ⁺¹ = mesh.centroid_coordinates[1] # next position
  x⃗cⁿ = mesh.centroid_coordinates[2] # current position
  x⃗cⁿ⁻¹ = mesh.centroid_coordinates[3] # previous position

  # update and shift all the stages
  for (cⁿ⁺¹, cⁿ, cⁿ⁻¹) in zip(
    StructArrays.components(x⃗cⁿ⁺¹),
    StructArrays.components(x⃗cⁿ),
    StructArrays.components(x⃗cⁿ⁻¹),
  )
    # cⁿ⁺¹ = x⃗cⁿ⁺¹[coord]
    # cⁿ = x⃗cⁿ[coord]
    # cⁿ⁻¹ = x⃗cⁿ⁻¹[coord]

    copy!(cⁿ⁻¹, cⁿ) # update n-1 from n
    copy!(cⁿ, cⁿ⁺¹) # update n from n+1
  end

  # and now compute the n+1 centroid coordinates from the nodes
  _centroid_coordinates!(x⃗cⁿ⁺¹, mesh.node_coordinates, mesh.iterators.cell.domain)

  # estimate the velocities of each centroid (xₜ, yₜ, zₜ)
  _compute_coord_velocities!(
    mesh.coord_velocities, mesh.centroid_coordinates, mesh.iterators.cell.domain, Δt
  )

  # now find ξₜ, ηₜ, ζₜ for all the centroids 
  _update_temporal_metrics!(mesh)

  # the edge interpolation is left to the particular 
  # discretization scheme (outside of this function)

  return nothing
end

#

"""
    compute_coord_velocities(coordinates, domain, Δt)

Compute the velocity of each of the coordinates. This is based on a multi-stage estimate
of the velocity (using the prev, current, and next position)
"""
function _compute_coord_velocities!(coord_velocities, centroid_coordinates, domain, Δt)
  ϕ = 1 / 2

  @assert isfinite(Δt) && Δt > 0
  backend = KernelAbstractions.get_backend(
    first(StructArrays.components(centroid_coordinates[1]))
  )

  kernel = _grid_vel_kernel(backend)

  # We want all threads (particularily on the GPU) to have
  # coalesced reads, so we only send the domain to the GPU
  # and not the halo regions. The offset index insures that
  # the 1st thread index starts at the correct position
  # in the domain
  I0 = first(domain) - oneunit(first(domain))

  x⃗cⁿ⁺¹ = centroid_coordinates[1] # next position
  x⃗cⁿ = centroid_coordinates[2] # current position
  x⃗cⁿ⁻¹ = centroid_coordinates[3] # previous position

  for (dxdt, xⁿ⁺¹, xⁿ, xⁿ⁻¹) in zip(
    StructArrays.components(coord_velocities),
    StructArrays.components(x⃗cⁿ⁺¹),
    StructArrays.components(x⃗cⁿ),
    StructArrays.components(x⃗cⁿ⁻¹),
  )
    kernel(dxdt, xⁿ⁺¹, xⁿ, xⁿ⁻¹, Δt, ϕ, I0; ndrange=size(domain))
  end

  KernelAbstractions.synchronize(backend)
  return nothing
end

@kernel function _grid_vel_kernel(dxdt, xⁿ⁺¹, xⁿ, xⁿ⁻¹, Δt, ϕ, I0)

  # apply offset to account for halo cells
  I = @index(Global, Cartesian)
  I += I0

  dxdt[I] = ((1 + ϕ) * xⁿ⁺¹[I] - (1 + 2ϕ)xⁿ[I] + ϕ * xⁿ⁻¹[I]) / Δt
end

#

function _update_temporal_metrics!(mesh::AbstractCurvilinearGrid1D)
  backend = KernelAbstractions.get_backend(
    first(StructArrays.components(mesh.centroid_coordinates[1]))
  )
  kernel = _temporal_metric_kernel(backend)

  domain = mesh.iterators.cell.domain

  # We want all threads (particularily on the GPU) to have
  # coalesced reads, so we only send the domain to the GPU
  # and not the halo regions. The offset index insures that
  # the 1st thread index starts at the correct position
  # in the domain
  I0 = first(domain) - oneunit(first(domain))

  for dim in (:ξ,)
    ξt = mesh.cell_center_metrics[dim].t
    ξx = mesh.cell_center_metrics[dim].x₁
    dxdt = mesh.coord_velocities.x
    kernel(ξt, ξx, dxdt, I0; ndrange=size(domain))
  end

  KernelAbstractions.synchronize(backend)
  return nothing
end

function _update_temporal_metrics!(mesh::AbstractCurvilinearGrid2D)
  backend = KernelAbstractions.get_backend(
    first(StructArrays.components(mesh.centroid_coordinates[1]))
  )
  kernel = _temporal_metric_kernel(backend)

  domain = mesh.iterators.cell.domain

  # We want all threads (particularily on the GPU) to have
  # coalesced reads, so we only send the domain to the GPU
  # and not the halo regions. The offset index insures that
  # the 1st thread index starts at the correct position
  # in the domain
  I0 = first(domain) - oneunit(first(domain))

  for dim in (:ξ, :η)
    ξt = mesh.cell_center_metrics[dim].t
    ξx = mesh.cell_center_metrics[dim].x₁
    ξy = mesh.cell_center_metrics[dim].x₂
    dxdt = mesh.coord_velocities.x
    dydt = mesh.coord_velocities.y
    kernel(ξt, ξx, ξy, dxdt, dydt, I0; ndrange=size(domain))
  end

  KernelAbstractions.synchronize(backend)
  return nothing
end

function _update_temporal_metrics!(mesh::CurvilinearGrid3D)
  backend = KernelAbstractions.get_backend(
    first(StructArrays.components(mesh.centroid_coordinates[1]))
  )
  kernel = _temporal_metric_kernel(backend)

  domain = mesh.iterators.cell.domain

  # We want all threads (particularily on the GPU) to have
  # coalesced reads, so we only send the domain to the GPU
  # and not the halo regions. The offset index insures that
  # the 1st thread index starts at the correct position
  # in the domain
  I0 = first(domain) - oneunit(first(domain))

  for dim in (:ξ, :η, :ζ)
    ξt = mesh.cell_center_metrics[dim].t
    ξx = mesh.cell_center_metrics[dim].x₁
    ξy = mesh.cell_center_metrics[dim].x₂
    ξz = mesh.cell_center_metrics[dim].x₃
    dxdt = mesh.coord_velocities.x
    dydt = mesh.coord_velocities.y
    dzdt = mesh.coord_velocities.z
    kernel(ξt, ξx, ξy, ξz, dxdt, dydt, dzdt, I0; ndrange=size(domain))
  end

  KernelAbstractions.synchronize(backend)
  return nothing
end

@kernel function _temporal_metric_kernel(ξt, ξx, dxdt, I0)

  # apply offset to account for halo cells
  I = @index(Global, Cartesian)
  I += I0

  ξt[I] = -ξx[I] * dxdt[I]
end

@kernel function _temporal_metric_kernel(ξt, ξx, ξy, dxdt, dydt, I0)

  # apply offset to account for halo cells
  I = @index(Global, Cartesian)
  I += I0

  ξt[I] = -(ξx[I] * dxdt[I] + ξy[I] * dydt[I])
end

@kernel function _temporal_metric_kernel(ξt, ξx, ξy, ξz, dxdt, dydt, dzdt, I0)

  # apply offset to account for halo cells
  I = @index(Global, Cartesian)
  I += I0

  ξt[I] = -(ξx[I] * dxdt[I] + ξy[I] * dydt[I] + ξz[I] * dzdt[I])
end