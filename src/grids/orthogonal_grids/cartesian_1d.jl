struct CartesianOrthogonalGrid1D{NC,CC,CV,I,DL,FA} <: AbstractCurvilinearGrid1D
  node_coordinates::NC
  centroid_coordinates::CC
  cell_volumes::CV
  iterators::I
  domain_limits::DL
  face_areas::FA
  nhalo::Int
end

function CartesianOrthogonalGrid1D(
  _x::AbstractVector{T}, nhalo::Int, backend; halo_coords_included=false
) where {T<:Real}
  x, limits, iters, nodedims, celldims = _prepare_1d_coordinates(
    _x, nhalo, halo_coords_included
  )

  node_coordinates = (; x=KernelAbstractions.zeros(backend, T, nodedims...))
  centroid_coordinates = (; x=KernelAbstractions.zeros(backend, T, celldims...))
  cell_volumes = KernelAbstractions.zeros(backend, T, celldims...)
  face_areas = (; i₊½=KernelAbstractions.zeros(backend, T, celldims...))

  _populate_1d_nodes!(node_coordinates.x, x, iters, halo_coords_included)

  compute_cartesian_centroids!(
    centroid_coordinates, node_coordinates, iters, backend, halo_coords_included
  )
  compute_cartesian_volumes!(
    cell_volumes, node_coordinates, iters, backend, halo_coords_included
  )
  compute_cartesian_face_areas!(
    face_areas, node_coordinates, iters, backend, halo_coords_included, nhalo
  )

  return CartesianOrthogonalGrid1D(
    node_coordinates, centroid_coordinates, cell_volumes, iters, limits, face_areas, nhalo
  )
end

function compute_cartesian_centroids!(
  centroids, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_cartesian_centroids!(backend)(
    centroids.x, node_coordinates.x, domain; ndrange=size(domain)
  )
  return nothing
end

function compute_cartesian_face_areas!(
  face_areas, node_coordinates, iters, backend, halo_coords_included, nhalo
)
  domain = iters.cell.domain
  if halo_coords_included
    offset = nhalo
    i₊½_domain = expand_upper(expand_lower(domain, 1, offset), 1, offset - 1)
  else
    offset = 1
    i₊½_domain = expand_lower(domain, 1, offset)
  end

  _compute_cartesian_face_areas!(backend)(
    face_areas.i₊½, node_coordinates.x, i₊½_domain; ndrange=size(i₊½_domain)
  )
  return nothing
end

function compute_cartesian_volumes!(
  volumes, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_cartesian_volumes!(backend)(
    volumes, node_coordinates.x, domain; ndrange=size(domain)
  )
  return nothing
end

@kernel function _compute_cartesian_centroids!(xc, xnode, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, = I.I
  xc[i] = (xnode[i + 1] + xnode[i]) / 2
end

@kernel function _compute_cartesian_face_areas!(Aᵢ₊½, x, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, = I.I
  Aᵢ₊½[I] = one(x[1])
end

@kernel function _compute_cartesian_volumes!(V, x, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, = I.I
  V[I] = x[i + 1] - x[i]
end
