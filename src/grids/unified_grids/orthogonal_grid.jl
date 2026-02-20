#
# OrthogonalGrid
#

mutable struct OrthogonalGrid{N,T,L,CS<:CoordinateSystemTrait,GC} <: AbstractUnifiedGrid
  legacy::L
  geometry_cache::GC
end

function OrthogonalGrid(
  legacy::Union{
    CartesianOrthogonalGrid1D,
    CylindricalOrthogonalGrid1D,
    SphericalOrthogonalGrid1D,
    AxisymmetricOrthogonalGrid2D,
    SphericalGrid3D,
  };
  coordinate_system::CoordinateSystemTrait=_coordinate_system_from_legacy(legacy),
  geometry_cache=nothing,
)
  N = ndims(legacy.iterators.cell.full)
  T = eltype(legacy)
  OrthogonalGrid{N,T,typeof(legacy),typeof(coordinate_system),typeof(geometry_cache)}(
    legacy, geometry_cache
  )
end

function OrthogonalGrid(
  x::AbstractVector{T},
  nhalo::Int;
  backend=CPU(),
  coordinate_system::CoordinateSystemTrait=CartesianCS(),
  halo_coords_included=false,
) where {T<:Real}
  legacy = if coordinate_system isa CartesianCS
    CartesianOrthogonalGrid1D(x, nhalo, backend; halo_coords_included=halo_coords_included)
  elseif coordinate_system isa CylindricalCS
    CylindricalOrthogonalGrid1D(x, nhalo, backend; halo_coords_included=halo_coords_included)
  elseif coordinate_system isa SphericalCS
    SphericalOrthogonalGrid1D(x, nhalo, backend; halo_coords_included=halo_coords_included)
  else
    throw(
      ArgumentError(
        "Unsupported coordinate system for 1D orthogonal grid: $(typeof(coordinate_system))",
      ),
    )
  end

  OrthogonalGrid(legacy; coordinate_system=coordinate_system)
end

function OrthogonalGrid(
  r::AbstractVector{T},
  z::AbstractVector{T},
  nhalo::Int;
  backend=CPU(),
  coordinate_system::CoordinateSystemTrait=AxisymmetricCS{:y}(),
  halo_coords_included=false,
) where {T<:Real}
  if !(coordinate_system isa AxisymmetricCS)
    throw(ArgumentError("2D orthogonal constructor currently supports only axisymmetric CS."))
  end
  legacy = AxisymmetricOrthogonalGrid2D(
    r, z, nhalo, backend; halo_coords_included=halo_coords_included
  )
  OrthogonalGrid(legacy; coordinate_system=coordinate_system)
end

function OrthogonalGrid(
  r::AbstractVector{T},
  theta::AbstractVector{T},
  phi::AbstractVector{T},
  nhalo::Int;
  backend=CPU(),
  coordinate_system::CoordinateSystemTrait=SphericalCS(),
  halo_coords_included=false,
) where {T<:Real}
  if !(coordinate_system isa SphericalCS)
    throw(ArgumentError("3D orthogonal constructor currently supports only spherical CS."))
  end
  legacy = SphericalGrid3D(
    r, theta, phi, nhalo, backend; halo_coords_included=halo_coords_included
  )
  OrthogonalGrid(legacy; coordinate_system=coordinate_system)
end
