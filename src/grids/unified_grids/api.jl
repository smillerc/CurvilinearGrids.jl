#
# Cache API
#

#
# Trait helpers
#

coordinate_system(grid::AbstractUnifiedGrid) = grid.coordinate_system_trait
basis_trait(grid::Union{MappedGrid,DiscreteGrid}) = grid.basis_vector_trait
function basis_trait(::OrthogonalGrid)
  throw(ArgumentError("`basis_trait` is undefined for `OrthogonalGrid`."))
end

coordinate_system(::Type{<:MappedGrid{L,CS}}) where {L,CS} = CS()
coordinate_system(::Type{<:DiscreteGrid{L,CS}}) where {L,CS} = CS()
coordinate_system(::Type{<:OrthogonalGrid{L,CS}}) where {L,CS} = CS()

basis_trait(::Type{<:MappedGrid{L,CS,BT}}) where {L,CS,BT} = BT()
basis_trait(::Type{<:DiscreteGrid{L,CS,BT}}) where {L,CS,BT} = BT()
function basis_trait(::Type{<:OrthogonalGrid})
  throw(ArgumentError("`basis_trait` is undefined for `OrthogonalGrid`."))
end

function invalidate_cell_metrics!(grid::AbstractMappedOrDiscreteGrid)
  grid.metric_caches.cell.valid = false
  grid.metric_caches.cell.data = nothing
  return nothing
end

function invalidate_face_metrics!(grid::AbstractMappedOrDiscreteGrid)
  grid.metric_caches.face.valid = false
  grid.metric_caches.face.data = nothing
  return nothing
end

function refresh_cell_metrics!(grid::AbstractMappedOrDiscreteGrid; include_halo_region::Bool=false)
  _refresh_cell_metrics!(grid; include_halo_region=include_halo_region)
end

function refresh_face_metrics!(grid::AbstractMappedOrDiscreteGrid; include_halo_region::Bool=false)
  _refresh_face_metrics!(grid; include_halo_region=include_halo_region)
end

function cell_metrics(grid::AbstractMappedOrDiscreteGrid; refresh::Bool=false)
  if refresh || !grid.metric_caches.cell.valid
    return refresh_cell_metrics!(grid)
  end
  return grid.metric_caches.cell.data
end

function face_metrics(grid::AbstractMappedOrDiscreteGrid; refresh::Bool=false)
  if refresh || !grid.metric_caches.face.valid
    return refresh_face_metrics!(grid)
  end
  return grid.metric_caches.face.data
end

function cell_metrics(::OrthogonalGrid; refresh::Bool=false)
  throw(ArgumentError("`cell_metrics` is undefined for `OrthogonalGrid`."))
end

function face_metrics(::OrthogonalGrid; refresh::Bool=false)
  throw(ArgumentError("`face_metrics` is undefined for `OrthogonalGrid`."))
end

#
# Adapter access and delegated geometry API
#

legacy_grid(grid::OrthogonalGrid) = grid.legacy
legacy_grid(grid::MappedGrid) = grid.core
legacy_grid(grid::DiscreteGrid) = grid.core

coords(grid::MappedGrid) = coords(grid.core)
coords(grid::DiscreteGrid) = coords(grid.core)
coords(grid::OrthogonalGrid) = coords(grid.legacy)

coord(grid::MappedGrid, idx) = coord(grid.core, idx)
coord(grid::DiscreteGrid, idx) = coord(grid.core, idx)
coord(grid::OrthogonalGrid, idx) = coord(grid.legacy, idx)

centroids(grid::MappedGrid) = centroids(grid.core)
centroids(grid::DiscreteGrid) = centroids(grid.core)
centroids(grid::OrthogonalGrid) = centroids(grid.legacy)

centroid(grid::MappedGrid, idx) = centroid(grid.core, idx)
centroid(grid::DiscreteGrid, idx) = centroid(grid.core, idx)
centroid(grid::OrthogonalGrid, idx) = centroid(grid.legacy, idx)

cellvolume(grid::MappedGrid, idx) = cellvolume(grid.core, idx)
cellvolume(grid::DiscreteGrid, idx) = cellvolume(grid.core, idx)
cellvolume(grid::OrthogonalGrid, idx) = cellvolume(grid.legacy, idx)

cellvolumes(grid::MappedGrid) = cellvolumes(grid.core)
cellvolumes(grid::DiscreteGrid) = cellvolumes(grid.core)
cellvolumes(grid::OrthogonalGrid) = cellvolumes(grid.legacy)

cellsize(grid::MappedGrid) = cellsize(grid.core)
cellsize(grid::DiscreteGrid) = cellsize(grid.core)
cellsize(grid::OrthogonalGrid) = cellsize(grid.legacy)

cellsize_withhalo(grid::MappedGrid) = cellsize_withhalo(grid.core)
cellsize_withhalo(grid::DiscreteGrid) = cellsize_withhalo(grid.core)
cellsize_withhalo(grid::OrthogonalGrid) = cellsize_withhalo(grid.legacy)

function jacobian_matrix(grid::AbstractMappedOrDiscreteGrid, idx)
  jacobian_matrix(grid.core, idx)
end

forward_cell_metrics(grid::AbstractMappedOrDiscreteGrid, idx) = forward_cell_metrics(
  grid.core, idx
)
inverse_cell_metrics(grid::AbstractMappedOrDiscreteGrid, idx) = inverse_cell_metrics(
  grid.core, idx
)
