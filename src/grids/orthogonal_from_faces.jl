function _orthogonal_grid_from_faces(
  ::CartesianCS, faces::NTuple{1,Vector{T}}, halo::Int, backend
) where {T}
  return CartesianOrthogonalGrid1D(faces[1], halo, backend; halo_coords_included=true)
end

function _orthogonal_grid_from_faces(
  ::CylindricalCS, faces::NTuple{1,Vector{T}}, halo::Int, backend
) where {T}
  return CylindricalOrthogonalGrid1D(faces[1], halo, backend; halo_coords_included=true)
end

function _orthogonal_grid_from_faces(
  ::SphericalCS, faces::NTuple{1,Vector{T}}, halo::Int, backend
) where {T}
  return SphericalOrthogonalGrid1D(faces[1], halo, backend; halo_coords_included=true)
end

function _orthogonal_grid_from_faces(
  ::CartesianCS, faces::NTuple{2,Vector{T}}, halo::Int, backend
) where {T}
  return CartesianOrthogonalGrid2D(faces[1], faces[2], halo, backend; halo_coords_included=true)
end

function _orthogonal_grid_from_faces(
  ::AxisymmetricCS{:x}, faces::NTuple{2,Vector{T}}, halo::Int, backend
) where {T}
  return AxisymmetricOrthogonalGrid2D(
    faces[1], faces[2], halo, backend; halo_coords_included=true, rotational_axis=:x
  )
end

function _orthogonal_grid_from_faces(
  ::AxisymmetricCS{:y}, faces::NTuple{2,Vector{T}}, halo::Int, backend
) where {T}
  return AxisymmetricOrthogonalGrid2D(
    faces[1], faces[2], halo, backend; halo_coords_included=true, rotational_axis=:y
  )
end

function _orthogonal_grid_from_faces(
  ::SphericalCS, faces::NTuple{2,Vector{T}}, halo::Int, backend
) where {T}
  return SphericalOrthogonalGrid2D(faces[1], faces[2], halo, backend; halo_coords_included=true)
end

function _orthogonal_grid_from_faces(
  ::CartesianCS, faces::NTuple{3,Vector{T}}, halo::Int, backend
) where {T}
  return CartesianOrthogonalGrid3D(
    faces[1], faces[2], faces[3], halo, backend; halo_coords_included=true
  )
end

function _orthogonal_grid_from_faces(
  ::SphericalCS, faces::NTuple{3,Vector{T}}, halo::Int, backend
) where {T}
  return SphericalGrid3D(faces[1], faces[2], faces[3], halo, backend; halo_coords_included=true)
end

function _orthogonal_grid_from_faces(coordinate_system, faces, halo, backend)
  throw(
    ArgumentError(
      "Unsupported orthogonal grid coordinate system $(typeof(coordinate_system)) for shared local block construction",
    ),
  )
end
