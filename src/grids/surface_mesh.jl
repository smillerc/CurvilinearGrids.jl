using CartesianDomains
using LinearAlgebra
using StaticArrays: SMatrix, SVector
using WriteVTK

@inline function _normalize_boundary_symbol(boundary::Symbol)
  Symbol(lowercase(String(boundary)))
end

@inline function _boundary_axis_side(::Val{N}, boundary::Symbol) where {N}
  b = _normalize_boundary_symbol(boundary)
  axis, side = if b in (:ilo, :imin, :xlo, :xmin)
    (1, :lo)
  elseif b in (:ihi, :imax, :xhi, :xmax)
    (1, :hi)
  elseif b in (:jlo, :jmin, :ylo, :ymin)
    (2, :lo)
  elseif b in (:jhi, :jmax, :yhi, :ymax)
    (2, :hi)
  elseif b in (:klo, :kmin, :zlo, :zmin)
    (3, :lo)
  elseif b in (:khi, :kmax, :zhi, :zmax)
    (3, :hi)
  else
    throw(
      ArgumentError(
        "Unsupported boundary `$boundary`. Expected one of `:ilo/:ihi`, `:jlo/:jhi`, or `:klo/:khi`.",
      ),
    )
  end
  if axis > N
    throw(
      ArgumentError(
        "Boundary `$boundary` maps to axis $axis, which is invalid for $N-D grid."
      ),
    )
  end
  return axis, side
end

@inline function _tangential_ranges(ranges::NTuple{N,UnitRange{Int}}, axis::Int) where {N}
  ntuple(i -> ranges[i < axis ? i : i + 1], N - 1)
end

@inline function _build_surface_index(
  tangential::NTuple{M,Int}, axis::Int, axis_index::Int, ::Val{N}
) where {M,N}
  if M != N - 1
    throw(
      ArgumentError(
        "Invalid tangential tuple length $M for surface index in $N dimensions."
      ),
    )
  end
  full = ntuple(i -> begin
    if i == axis
      return axis_index
    end
    ti = i < axis ? i : i - 1
    return tangential[ti]
  end, N)
  return CartesianIndex(full)
end

@inline function _coord_to_cartesian(::CoordinateSystemTrait, q::SVector{N,T}) where {N,T}
  return q
end

@inline function _coord_to_cartesian(::CartesianCS, q::SVector{N,T}) where {N,T}
  return q
end

@inline function _coord_to_cartesian(::CurvilinearCS, q::SVector{N,T}) where {N,T}
  return q
end

@inline function _coord_to_cartesian(::AxisymmetricCS, q::SVector{2,T}) where {T}
  return q
end

@inline function _coord_to_cartesian(::CylindricalCS, q::SVector{2,T}) where {T}
  return q
end

@inline function _coord_to_cartesian(::SphericalCS, q::SVector{3,T}) where {T}
  r, θ, ϕ = q
  sθ = sin(θ)
  cθ = cos(θ)
  sϕ = sin(ϕ)
  cϕ = cos(ϕ)
  return SVector{3,T}(r * sθ * cϕ, r * sθ * sϕ, r * cθ)
end

@inline function _coord_to_cartesian_jacobian(
  ::CoordinateSystemTrait, q::SVector{N,T}
) where {N,T}
  return one(SMatrix{N,N,T})
end

@inline function _coord_to_cartesian_jacobian(::CartesianCS, q::SVector{N,T}) where {N,T}
  return one(SMatrix{N,N,T})
end

@inline function _coord_to_cartesian_jacobian(::CurvilinearCS, q::SVector{N,T}) where {N,T}
  return one(SMatrix{N,N,T})
end

@inline function _coord_to_cartesian_jacobian(::AxisymmetricCS, q::SVector{2,T}) where {T}
  return one(SMatrix{2,2,T})
end

@inline function _coord_to_cartesian_jacobian(::CylindricalCS, q::SVector{2,T}) where {T}
  return one(SMatrix{2,2,T})
end

@inline function _coord_to_cartesian_jacobian(::SphericalCS, q::SVector{3,T}) where {T}
  r, θ, ϕ = q
  sθ = sin(θ)
  cθ = cos(θ)
  sϕ = sin(ϕ)
  cϕ = cos(ϕ)
  return SMatrix{3,3,T,9}(
    sθ * cϕ,
    sθ * sϕ,
    cθ,
    r * cθ * cϕ,
    r * cθ * sϕ,
    -r * sθ,
    -r * sθ * sϕ,
    r * sθ * cϕ,
    zero(T),
  )
end

@inline function _area_vector_from_covariant(
  jacobian::SMatrix{2,2,T,4}, axis::Int
) where {T}
  if axis == 1
    return SVector{2,T}(jacobian[2, 2], -jacobian[1, 2])
  elseif axis == 2
    return SVector{2,T}(-jacobian[2, 1], jacobian[1, 1])
  end
  throw(ArgumentError("Invalid 2D face axis: $axis"))
end

@inline function _area_vector_from_covariant(
  jacobian::SMatrix{3,3,T,9}, axis::Int
) where {T}
  aξ = SVector{3,T}(jacobian[1, 1], jacobian[2, 1], jacobian[3, 1])
  aη = SVector{3,T}(jacobian[1, 2], jacobian[2, 2], jacobian[3, 2])
  aζ = SVector{3,T}(jacobian[1, 3], jacobian[2, 3], jacobian[3, 3])
  if axis == 1
    return cross(aη, aζ)
  elseif axis == 2
    return cross(aζ, aξ)
  elseif axis == 3
    return cross(aξ, aη)
  end
  throw(ArgumentError("Invalid 3D face axis: $axis"))
end

@inline function _cartesian_forward_jacobian(
  cs::CoordinateSystemTrait, q::SVector{N,T}, jacobian_qξ::SMatrix{N,N,T}
) where {N,T}
  A = _coord_to_cartesian_jacobian(cs, q)
  J = A * jacobian_qξ
  return SMatrix{N,N,T,N * N}(Tuple(J))
end

@inline _surface_area_scale(::CoordinateSystemTrait, q::SVector{N,T}) where {N,T} = one(T)
@inline _surface_area_scale(::CurvilinearCS, q::SVector{N,T}) where {N,T} = one(T)
@inline _surface_area_scale(::CartesianCS, q::SVector{N,T}) where {N,T} = one(T)
@inline _surface_area_scale(::SphericalCS, q::SVector{N,T}) where {N,T} = one(T)

@inline _rotational_radius(::CylindricalCS, q::SVector{2,T}) where {T} = q[1]
@inline _rotational_radius(::AxisymmetricCS{:x}, q::SVector{2,T}) where {T} = q[2]
@inline _rotational_radius(::AxisymmetricCS{:y}, q::SVector{2,T}) where {T} = q[1]
function _rotational_radius(::AxisymmetricCS{Axis}, ::SVector{2,T}) where {Axis,T}
  throw(
    ArgumentError(
      "Unsupported axisymmetric axis `:$Axis` for `SurfaceGrid`. Supported values are `:x` and `:y`.",
    ),
  )
end

@inline function _surface_area_scale(
  cs::Union{CylindricalCS,AxisymmetricCS}, q::SVector{2,T}
) where {T}
  return T(2π) * abs(_rotational_radius(cs, q))
end

@inline function _surface_orientation_sign(
  side::Symbol, outward_normal::Bool, ::Type{T}
) where {T}
  s = side === :hi ? one(T) : -one(T)
  return outward_normal ? s : -s
end

@inline function _normal_and_area(area_vec::SVector{N,T}) where {N,T}
  area = norm(area_vec)
  normal = area > zero(T) ? area_vec / area : zero(SVector{N,T})
  return normal, area
end

"""
    SurfaceGrid{N,T,...}

Boundary-extracted surface representation of a 2D/3D unified mapped/discrete
grid.

# Fields
  - `boundary`: Boundary selector symbol.
  - `outward_normal`: Whether normals are oriented outward (`true`) or inward (`false`).
  - `axis`: Boundary-normal computational axis index.
  - `side`: Boundary side symbol (`:lo` or `:hi`).
  - `node_coordinates`: Boundary node coordinates in Cartesian space.
  - `face_centers`: Boundary face-center coordinates in Cartesian space.
  - `face_normals`: Boundary unit normal vectors in Cartesian space.
  - `face_areas`: Physical boundary measure. For 3D this is true surface area. For
    2D this is edge length, except axisymmetric/cylindrical coordinate systems where
    it is the rotated physical area (`2πr` times edge length).
"""
struct SurfaceGrid{N,T,CS<:CoordinateSystemTrait,BT<:BasisTrait,NC,FC,FN,FA}
  boundary::Symbol
  outward_normal::Bool
  axis::Int
  side::Symbol
  coordinate_system::CS
  basis::BT
  node_coordinates::NC
  face_centers::FC
  face_normals::FN
  face_areas::FA
end

function SurfaceGrid(
  grid::CurvilinearGrids.OrthogonalGrid{N,T,CurvilinearGrids.CartesianCS},
  boundary::Symbol,
  outward_normal::Bool=true,
) where {N,T}
  if N != 2 && N != 3
    throw(ArgumentError("`SurfaceGrid` is only supported for 2D and 3D unified grids."))
  end

  axis, side = _boundary_axis_side(Val(N), boundary)
  node_ranges = Tuple(grid.iterators.node.domain.indices)
  cell_ranges = Tuple(grid.iterators.cell.domain.indices)
  node_axis_idx = side === :lo ? first(node_ranges[axis]) : last(node_ranges[axis])
  cell_axis_idx = side === :lo ? first(cell_ranges[axis]) : last(cell_ranges[axis])
  tangential_node_ranges = _tangential_ranges(node_ranges, axis)
  tangential_cell_ranges = _tangential_ranges(cell_ranges, axis)

  node_dims = ntuple(i -> length(tangential_node_ranges[i]), N - 1)
  cell_dims = ntuple(i -> length(tangential_cell_ranges[i]), N - 1)

  cs = coordinate_system(grid)
  bt = CartesianBasis()
  loc = GridTypes._face_location_symbol(axis, side, Val(N))
  normal_sign = _surface_orientation_sign(side, outward_normal, T)

  node_coordinates = ntuple(_ -> Array{T}(undef, node_dims...), N)
  face_centers = ntuple(_ -> Array{T}(undef, cell_dims...), N)
  face_normals = ntuple(_ -> Array{T}(undef, cell_dims...), N)
  face_areas = Array{T}(undef, cell_dims...)

  for (Ilocal, Itang) in
      zip(CartesianIndices(node_dims), CartesianIndices(tangential_node_ranges))
    Iparent = _build_surface_index(Itang.I, axis, node_axis_idx, Val(N))
    qnode = SVector{N,T}(ntuple(d -> T(grid.node_coordinates[d][Iparent[d]]), N))
    for d in 1:N
      node_coordinates[d][Ilocal] = qnode[d]
    end
  end

  for (Ilocal, Itang) in
      zip(CartesianIndices(cell_dims), CartesianIndices(tangential_cell_ranges))
    Icell = _build_surface_index(Itang.I, axis, cell_axis_idx, Val(N))
    center = CurvilinearGrids.face_coordinate(grid, Icell.I, loc)
    area = CurvilinearGrids.face_area(grid, Icell.I, loc)
    for d in 1:N
      face_centers[d][Ilocal] = center[d]
      face_normals[d][Ilocal] = d == axis ? normal_sign : zero(T)
    end
    face_areas[Ilocal] = area
  end

  return SurfaceGrid{
    N,
    T,
    typeof(cs),
    typeof(bt),
    typeof(node_coordinates),
    typeof(face_centers),
    typeof(face_normals),
    typeof(face_areas),
  }(
    boundary,
    outward_normal,
    axis,
    side,
    cs,
    bt,
    node_coordinates,
    face_centers,
    face_normals,
    face_areas,
  )
end

"""
    SurfaceGrid(grid::Union{MappedGrid,DiscreteGrid}, boundary::Symbol, outward_normal=true)

Extract a boundary surface from a 2D/3D mapped or discrete unified grid.

# Arguments
  - `grid`: Source `MappedGrid` or `DiscreteGrid` (`N=2` or `N=3`).
  - `boundary`: Boundary selector `:{i,j,k}{lo,hi}` (for example `:ilo`, `:jhi`, `:klo`).
  - `outward_normal`: If `true`, normals point outward; if `false`, inward.

# Returns
`SurfaceGrid` containing boundary geometry, unit normals, and face areas.
"""
function SurfaceGrid(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}},
  boundary::Symbol,
  outward_normal::Bool=true,
) where {N,T}
  if N != 2 && N != 3
    throw(ArgumentError("`SurfaceGrid` is only supported for 2D and 3D unified grids."))
  end

  axis, side = _boundary_axis_side(Val(N), boundary)
  node_ranges = Tuple(grid.iterators.node.domain.indices)
  cell_ranges = Tuple(grid.iterators.cell.domain.indices)
  node_axis_idx = side === :lo ? first(node_ranges[axis]) : last(node_ranges[axis])
  cell_axis_idx = side === :lo ? first(cell_ranges[axis]) : last(cell_ranges[axis])
  tangential_node_ranges = _tangential_ranges(node_ranges, axis)
  tangential_cell_ranges = _tangential_ranges(cell_ranges, axis)

  node_dims = ntuple(i -> length(tangential_node_ranges[i]), N - 1)
  cell_dims = ntuple(i -> length(tangential_cell_ranges[i]), N - 1)

  node_arrays = ntuple(d -> Array(grid.node_coordinates[d]), N)
  cs = coordinate_system(grid)
  bt = basis_trait(grid)
  loc = GridTypes._face_location_symbol(axis, side, Val(N))

  node_coordinates = ntuple(_ -> Array{T}(undef, node_dims...), N)
  face_centers = ntuple(_ -> Array{T}(undef, cell_dims...), N)
  face_normals = ntuple(_ -> Array{T}(undef, cell_dims...), N)
  face_areas = Array{T}(undef, cell_dims...)

  for (Ilocal, Itang) in
      zip(CartesianIndices(node_dims), CartesianIndices(tangential_node_ranges))
    Iparent = _build_surface_index(Itang.I, axis, node_axis_idx, Val(N))
    qnode = SVector{N,T}(ntuple(d -> T(node_arrays[d][Iparent]), N))
    xnode = _coord_to_cartesian(cs, qnode)
    for d in 1:N
      node_coordinates[d][Ilocal] = xnode[d]
    end
  end

  for (Ilocal, Itang) in
      zip(CartesianIndices(cell_dims), CartesianIndices(tangential_cell_ranges))
    Icell = _build_surface_index(Itang.I, axis, cell_axis_idx, Val(N))
    geom = GridTypes._face_geometry(grid, Icell.I, loc)
    center = geom.cartesian_coordinate
    normal_out = geom.normal
    normal = outward_normal ? normal_out : -normal_out
    area = geom.area

    for d in 1:N
      face_centers[d][Ilocal] = center[d]
      face_normals[d][Ilocal] = normal[d]
    end
    face_areas[Ilocal] = area
  end

  return SurfaceGrid{
    N,
    T,
    typeof(cs),
    typeof(bt),
    typeof(node_coordinates),
    typeof(face_centers),
    typeof(face_normals),
    typeof(face_areas),
  }(
    boundary,
    outward_normal,
    axis,
    side,
    cs,
    bt,
    node_coordinates,
    face_centers,
    face_normals,
    face_areas,
  )
end

"""
    extract_surface_mesh(
      grid::Union{MappedGrid,DiscreteGrid},
      boundary::Symbol;
      outward_normal=true,
    )

Construct a `SurfaceGrid` from a unified mapped/discrete grid boundary.
"""
function extract_surface_mesh(
  grid::Union{MappedGrid,DiscreteGrid}, boundary::Symbol; outward_normal::Bool=true
)
  return SurfaceGrid(grid, boundary, outward_normal)
end

function extract_surface_mesh(
  grid::CurvilinearGrids.OrthogonalGrid{N,T,CurvilinearGrids.CartesianCS},
  boundary::Symbol;
  outward_normal::Bool=true,
) where {N,T}
  return SurfaceGrid(grid, boundary, outward_normal)
end

"""
    extract_surface_mesh(mesh::AbstractCurvilinearGrid2D, loc::Symbol)

Extract one boundary edge from a legacy 2D mesh in physical coordinates.
"""
function extract_surface_mesh(mesh::AbstractCurvilinearGrid2D, loc::Symbol)
  i, j = (1, 2)
  full = mesh.iterators.node.domain

  if loc === :ilo
    domain = extract_from_lower(full, i, 1)
  elseif loc === :jlo
    domain = extract_from_lower(full, j, 1)
  elseif loc === :ihi
    domain = extract_from_upper(full, i, 1)
  elseif loc === :jhi
    domain = extract_from_upper(full, j, 1)
  else
    throw(ArgumentError("Unsupported 2D boundary `$loc`."))
  end

  x = Array(mesh.node_coordinates.x[domain])
  y = Array(mesh.node_coordinates.y[domain])

  return (x, y)
end

"""
    extract_surface_mesh(mesh::AbstractCurvilinearGrid3D, loc::Symbol)

Extract one boundary surface from a legacy 3D mesh in physical coordinates.
"""
function extract_surface_mesh(mesh::AbstractCurvilinearGrid3D, loc::Symbol)
  i, j, k = (1, 2, 3)
  full = mesh.iterators.node.domain

  if loc === :ilo
    domain = extract_from_lower(full, i, 1)
  elseif loc === :jlo
    domain = extract_from_lower(full, j, 1)
  elseif loc === :klo
    domain = extract_from_lower(full, k, 1)
  elseif loc === :ihi
    domain = extract_from_upper(full, i, 1)
  elseif loc === :jhi
    domain = extract_from_upper(full, j, 1)
  elseif loc === :khi
    domain = extract_from_upper(full, k, 1)
  else
    throw(ArgumentError("Unsupported 3D boundary `$loc`."))
  end

  x = Array(mesh.node_coordinates.x[domain])
  y = Array(mesh.node_coordinates.y[domain])
  z = Array(mesh.node_coordinates.z[domain])

  return (x, y, z)
end

"""
    save_vtk(surface::SurfaceGrid{2}, fn="surface")

Write a 2D boundary `SurfaceGrid` (1D line) to VTK with only node coordinates,
face normals, and face areas.
"""
@inline _surface_component_names(::Val{1}) = ["x1"]
@inline _surface_component_names(::Val{2}) = ["x1", "x2"]
@inline _surface_component_names(::Val{3}) = ["x1", "x2", "x3"]

@inline _surface_scalar_cell_field(data::AbstractVector, ::SurfaceGrid{2}) = reshape(data, :, 1)
@inline _surface_scalar_cell_field(data::AbstractArray, ::SurfaceGrid{3}) = data

@inline _surface_vector_cell_field(data::NTuple{N,<:AbstractVector}, ::SurfaceGrid{2}) where {N} =
  ntuple(i -> reshape(data[i], :, 1), N)
@inline _surface_vector_cell_field(data::NTuple{N,<:AbstractArray}, ::SurfaceGrid{3}) where {N} =
  data

function _write_surface_extra_cell_data!(vtk, surface::SurfaceGrid, extra_cell_data)
  isnothing(extra_cell_data) && return nothing
  for (key, value) in pairs(extra_cell_data)
    if value isa Tuple
      ncomponents = length(value)
      vtk[String(key), VTKCellData(), component_names = _surface_component_names(Val(ncomponents))] =
        _surface_vector_cell_field(value, surface)
    else
      vtk[String(key), VTKCellData()] = _surface_scalar_cell_field(value, surface)
    end
  end
  return nothing
end

function VTKOutput.save_vtk(surface::SurfaceGrid{2}, fn="surface"; extra_cell_data=nothing)
  x_line, y_line = surface.node_coordinates
  npoints = length(x_line)
  x = Array{eltype(x_line)}(undef, npoints, 2)
  y = Array{eltype(y_line)}(undef, npoints, 2)
  @views begin
    x[:, 1] .= x_line
    x[:, 2] .= x_line
    y[:, 1] .= y_line
    y[:, 2] .= y_line
  end
  area = reshape(surface.face_areas, :, 1)
  nx = reshape(surface.face_normals[1], :, 1)
  ny = reshape(surface.face_normals[2], :, 1)

  @info "Writing to $fn.vti"
  vtk_grid(fn, (x, y)) do vtk
    vtk["face_area", VTKCellData()] = area
    vtk["face_normal", VTKCellData(), component_names = ["x1", "x2"]] = (nx, ny)
    _write_surface_extra_cell_data!(vtk, surface, extra_cell_data)
  end
  return nothing
end

function _surface_vtk_points(surface::SurfaceGrid{3})
  x, y, z = surface.node_coordinates
  n1, n2 = size(x)
  points = Array{SVector{3,eltype(x)}}(undef, n1, n2)
  @inbounds for j in 1:n2
    for i in 1:n1
      points[i, j] = SVector(x[i, j], y[i, j], z[i, j])
    end
  end
  return points
end

"""
    save_vtk(surface::SurfaceGrid{3}, fn="surface")

Write a 3D boundary `SurfaceGrid` (2D surface) to VTK with only node
coordinates, face normals, and face areas.
"""
function VTKOutput.save_vtk(surface::SurfaceGrid{3}, fn="surface"; extra_cell_data=nothing)
  points = _surface_vtk_points(surface)
  @info "Writing to $fn.vti"
  vtk_grid(fn, points) do vtk
    vtk["face_area", VTKCellData()] = surface.face_areas
    vtk["face_normal", VTKCellData(), component_names = ["x1", "x2", "x3"]] = (
      surface.face_normals[1], surface.face_normals[2], surface.face_normals[3]
    )
    _write_surface_extra_cell_data!(vtk, surface, extra_cell_data)
  end
  return nothing
end
