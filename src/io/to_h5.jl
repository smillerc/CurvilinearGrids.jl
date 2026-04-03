module H5Output

using HDF5, Unitful
using Serialization
using ..GridTypes

export write_coordinates, read_coordinates

@inline function _serialize_payload(x)
  io = IOBuffer()
  Serialization.serialize(io, x)
  return take!(io)
end

@inline function _deserialize_payload(data::AbstractVector{UInt8})
  io = IOBuffer(data)
  return Serialization.deserialize(io)
end

@inline function _coordinate_system_tag(cs::CoordinateSystemTrait)
  if cs isa CartesianCS
    return "CartesianCS"
  elseif cs isa CylindricalCS
    return "CylindricalCS"
  elseif cs isa SphericalCS
    return "SphericalCS"
  elseif cs isa CurvilinearCS
    return "CurvilinearCS"
  elseif cs isa AxisymmetricCS{:x}
    return "AxisymmetricCS{:x}"
  elseif cs isa AxisymmetricCS{:y}
    return "AxisymmetricCS{:y}"
  end
  throw(
    ArgumentError("Unsupported coordinate system trait $(typeof(cs)) for serialization.")
  )
end

@inline function _basis_trait_tag(bt::BasisTrait)
  if bt isa CartesianBasis
    return "CartesianBasis"
  elseif bt isa SphericalBasis
    return "SphericalBasis"
  end
  throw(ArgumentError("Unsupported basis trait $(typeof(bt)) for serialization."))
end

@inline function _coordinate_system_from_tag(tag::AbstractString)
  if tag == "CartesianCS"
    return CartesianCS()
  elseif tag == "CylindricalCS"
    return CylindricalCS()
  elseif tag == "SphericalCS"
    return SphericalCS()
  elseif tag == "CurvilinearCS"
    return CurvilinearCS()
  elseif tag == "AxisymmetricCS{:x}"
    return AxisymmetricCS{:x}()
  elseif tag == "AxisymmetricCS{:y}"
    return AxisymmetricCS{:y}()
  end
  throw(ArgumentError("Unsupported serialized coordinate system tag `$tag`."))
end

@inline function _basis_trait_from_tag(tag::AbstractString)
  if tag == "CartesianBasis"
    return CartesianBasis()
  elseif tag == "SphericalBasis"
    return SphericalBasis()
  end
  throw(ArgumentError("Unsupported serialized basis tag `$tag`."))
end

@inline function _write_discretization_scheme(filename::String, mesh)
  if hasproperty(mesh, :discretization_scheme_name)
    h5write(
      filename,
      "discretization_scheme",
      String(getproperty(mesh, :discretization_scheme_name)),
    )
  end
end

function _read_discretization_scheme(filename::String; default::Symbol=:meg6)
  grid_file = h5open(filename, "r")
  scheme = if haskey(grid_file, "discretization_scheme")
    Symbol(read(grid_file, "discretization_scheme"))
  else
    default
  end
  close(grid_file)
  return scheme
end

### Write h5 functions

"""
    write_coordinates(mesh::CurvilinearGrid1D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A})

Write 1D curvilinear grid information to HDF5 with coordinate units recorded in `units`.
"""
function write_coordinates(
  mesh::CurvilinearGrid1D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)
    _write_discretization_scheme(filename, mesh)
  end
end

"""
    write_coordinates(mesh::CurvilinearGrid2D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A})

Write 2D curvilinear grid information to HDF5 with coordinate units recorded in `units`.
"""
function write_coordinates(
  mesh::CurvilinearGrid2D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x, y = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x'))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y'))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)
    _write_discretization_scheme(filename, mesh)
  end
end

"""
    write_coordinates(mesh::CurvilinearGrid3D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A})

Write 3D curvilinear grid information to HDF5 with coordinate units recorded in `units`.
"""
function write_coordinates(
  mesh::CurvilinearGrid3D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x, y, z = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "z", collect(z))
    h5writeattr(filename, "z", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)
    _write_discretization_scheme(filename, mesh)
  end
end

function write_coordinates(
  mesh::DiscreteGrid{1}, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x = coords(mesh)

  h5open(filename, "w") do file
    h5write(filename, "x", collect(x))
    h5writeattr(filename, "x", Dict("Units" => string(units)))
    h5write(filename, "grid_type", "DiscreteGrid1D")
  end
end

function write_coordinates(
  mesh::DiscreteGrid{2}, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x, y = coords(mesh)

  h5open(filename, "w") do file
    h5write(filename, "x", collect(x'))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y'))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "grid_type", "DiscreteGrid2D")
  end
end

function write_coordinates(
  mesh::DiscreteGrid{3}, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x, y, z = coords(mesh)

  h5open(filename, "w") do file
    h5write(filename, "x", collect(x))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "z", collect(z))
    h5writeattr(filename, "z", Dict("Units" => string(units)))

    h5write(filename, "grid_type", "DiscreteGrid3D")
  end
end

function write_coordinates(
  mesh::MappedGrid{1}, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x = coords(mesh)
  state = mesh.state[]
  t = state isa NamedTuple && haskey(state, :t) ? state.t : zero(eltype(mesh))
  params = state isa NamedTuple && haskey(state, :params) ? state.params : (;)
  celldims = collect(size(mesh.iterators.cell.domain))
  mapping_bytes = _serialize_payload(mesh.mapping_functions)
  params_bytes = _serialize_payload(params)

  h5open(filename, "w") do file
    h5write(filename, "x", collect(x))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "grid_type", "MappedGrid1D")
    h5write(filename, "nhalo", mesh.nhalo)
    h5write(filename, "time", t)
    h5write(filename, "celldims", celldims)
    h5write(
      filename, "coordinate_system_trait", _coordinate_system_tag(coordinate_system(mesh))
    )
    h5write(filename, "basis_trait", _basis_trait_tag(basis_trait(mesh)))
    h5write(filename, "mapping_functions_serialized", mapping_bytes)
    h5write(filename, "mapping_params_serialized", params_bytes)
  end
end

function write_coordinates(
  mesh::MappedGrid{2}, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x, y = coords(mesh)
  state = mesh.state[]
  t = state isa NamedTuple && haskey(state, :t) ? state.t : zero(eltype(mesh))
  params = state isa NamedTuple && haskey(state, :params) ? state.params : (;)
  celldims = collect(size(mesh.iterators.cell.domain))
  mapping_bytes = _serialize_payload(mesh.mapping_functions)
  params_bytes = _serialize_payload(params)

  h5open(filename, "w") do file
    h5write(filename, "x", collect(x'))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y'))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "grid_type", "MappedGrid2D")
    h5write(filename, "nhalo", mesh.nhalo)
    h5write(filename, "time", t)
    h5write(filename, "celldims", celldims)
    h5write(
      filename, "coordinate_system_trait", _coordinate_system_tag(coordinate_system(mesh))
    )
    h5write(filename, "basis_trait", _basis_trait_tag(basis_trait(mesh)))
    h5write(filename, "mapping_functions_serialized", mapping_bytes)
    h5write(filename, "mapping_params_serialized", params_bytes)
  end
end

function write_coordinates(
  mesh::MappedGrid{3}, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x, y, z = coords(mesh)
  state = mesh.state[]
  t = state isa NamedTuple && haskey(state, :t) ? state.t : zero(eltype(mesh))
  params = state isa NamedTuple && haskey(state, :params) ? state.params : (;)
  celldims = collect(size(mesh.iterators.cell.domain))
  mapping_bytes = _serialize_payload(mesh.mapping_functions)
  params_bytes = _serialize_payload(params)

  h5open(filename, "w") do file
    h5write(filename, "x", collect(x))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "z", collect(z))
    h5writeattr(filename, "z", Dict("Units" => string(units)))

    h5write(filename, "grid_type", "MappedGrid3D")
    h5write(filename, "nhalo", mesh.nhalo)
    h5write(filename, "time", t)
    h5write(filename, "celldims", celldims)
    h5write(
      filename, "coordinate_system_trait", _coordinate_system_tag(coordinate_system(mesh))
    )
    h5write(filename, "basis_trait", _basis_trait_tag(basis_trait(mesh)))
    h5write(filename, "mapping_functions_serialized", mapping_bytes)
    h5write(filename, "mapping_params_serialized", params_bytes)
  end
end

"""
    write_coordinates(mesh::UniformGrid1D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A})

Write 1D uniform grid information to HDF5 with coordinate units recorded in `units`.
"""
function write_coordinates(
  mesh::UniformGrid1D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)
    _write_discretization_scheme(filename, mesh)
  end
end

"""
    write_coordinates(mesh::UniformGrid2D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A})

Write 2D uniform grid information to HDF5 with coordinate units recorded in `units`.
"""
function write_coordinates(
  mesh::UniformGrid2D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x, y = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x'))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y'))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)
    _write_discretization_scheme(filename, mesh)
  end
end

"""
    write_coordinates(mesh::UniformGrid3D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A})

Write 3D uniform grid information to HDF5 with coordinate units recorded in `units`.
"""
function write_coordinates(
  mesh::UniformGrid3D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x, y, z = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "z", collect(z))
    h5writeattr(filename, "z", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)
    _write_discretization_scheme(filename, mesh)
  end
end

function write_coordinates(
  mesh::OrthogonalGrid{3,T,SphericalCS},
  filename::String,
  units::Unitful.FreeUnits{N,Unitful.𝐋,A},
) where {N,A,T}
  r, θ, ϕ = coords(mesh)
  mesh_type = "SphericalGrid3D"

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "r", collect(r))
    h5writeattr(filename, "r", Dict("Units" => string(units)))

    h5write(filename, "theta", collect(θ))
    h5writeattr(filename, "theta", Dict("Units" => "rad"))

    h5write(filename, "phi", collect(ϕ))
    h5writeattr(filename, "phi", Dict("Units" => "rad"))

    h5write(filename, "nhalo", mesh.nhalo)
    h5writeattr(filename, "nhalo", Dict("Description" => "Number of halo cells"))

    h5write(filename, "grid_type", mesh_type)
  end
end

"""
    write_coordinates(mesh::RectilinearGrid2D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A})

Write 2D rectilinear grid information to HDF5 with coordinate units recorded in `units`.
"""
function write_coordinates(
  mesh::RectilinearGrid2D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x, y = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x'))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y'))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)
    _write_discretization_scheme(filename, mesh)
  end
end

"""
    write_coordinates(mesh::RectilinearGrid3D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A})

Write 3D rectilinear grid information to HDF5 with coordinate units recorded in `units`.
"""
function write_coordinates(
  mesh::RectilinearGrid3D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x, y, z = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "z", collect(z))
    h5writeattr(filename, "z", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)
    _write_discretization_scheme(filename, mesh)
  end
end

"""
    write_coordinates(mesh::CylindricalGrid1D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A})

Write 1D cylindrical grid information to HDF5 with coordinate units recorded in `units`.
"""
function write_coordinates(
  mesh::CylindricalGrid1D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))
  snap_to_axis = mesh.snap_to_axis

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)

    h5write(filename, "snap_to_axis", snap_to_axis)
  end
end

"""
    write_coordinates(mesh::AxisymmetricGrid2D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A})

Write 2D axisymmetric grid information to HDF5 with coordinate units recorded in `units`, preserving the rotational axis information.
"""
function write_coordinates(
  mesh::AxisymmetricGrid2D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x, y = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))
  snap_to_axis = mesh.snap_to_axis
  rotational_axis = String(mesh.rotational_axis)

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x'))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y'))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)

    h5write(filename, "snap_to_axis", snap_to_axis)

    h5write(filename, "rotational_axis", rotational_axis)
  end
end

"""
    write_coordinates(mesh::SphericalGrid1D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A})

Write 1D spherical grid information to HDF5 with coordinate units recorded in `units`, including whether nodes are snapped to the axis.
"""
function write_coordinates(
  mesh::SphericalGrid1D, filename::String, units::Unitful.FreeUnits{N,Unitful.𝐋,A}
) where {N,A}
  x = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))
  snap_to_axis = mesh.snap_to_axis

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)

    h5write(filename, "snap_to_axis", snap_to_axis)
  end
end

### Read h5 functions

"""
    read_coordinates(filename::String, nhalo; compute_metrics::Bool=true, cache_mode::Symbol=:eager)

Load grid information from an HDF5 file written by `write_coordinates` and reconstruct the corresponding mesh using the chosen number
of halo cells `nhalo`.
"""
function read_coordinates(
  filename::String; compute_metrics::Bool=true, cache_mode::Symbol=:eager
)
  read_coordinates(filename, 0; compute_metrics=compute_metrics, cache_mode=cache_mode)
end

function read_coordinates(
  filename::String, nhalo; compute_metrics::Bool=true, cache_mode::Symbol=:eager
)
  grid_file = h5open(filename, "r")
  grid_type = read(grid_file, "grid_type")
  close(grid_file)

  if grid_type == "CurvilinearGrid1D"
    x = read_CurvilinearGrid1D(filename)
    mesh = CurvilinearGrid1D(x, _read_discretization_scheme(filename))
  elseif grid_type == "CurvilinearGrid2D"
    x, y = read_CurvilinearGrid2D(filename)
    mesh = CurvilinearGrid2D(x, y, _read_discretization_scheme(filename))
  elseif grid_type == "CurvilinearGrid3D"
    x, y, z = read_CurvilinearGrid3D(filename)
    mesh = CurvilinearGrid3D(x, y, z, _read_discretization_scheme(filename))
  elseif grid_type == "MappedGrid1D"
    mapping_functions, params, celldims, t, nhalo, cs, bt = _read_MappedGrid_payload(
      filename, Val(1)
    )
    mesh = MappedGrid(
      mapping_functions.x1,
      params,
      celldims,
      nhalo;
      t=t,
      compute_metrics=compute_metrics,
      cache_mode=cache_mode,
      coordinate_system=cs,
      basis=bt,
    )
  elseif grid_type == "MappedGrid2D"
    mapping_functions, params, celldims, t, nhalo, cs, bt = _read_MappedGrid_payload(
      filename, Val(2)
    )
    mesh = MappedGrid(
      mapping_functions.x1,
      mapping_functions.x2,
      params,
      celldims,
      nhalo;
      t=t,
      compute_metrics=compute_metrics,
      cache_mode=cache_mode,
      coordinate_system=cs,
      basis=bt,
    )
  elseif grid_type == "MappedGrid3D"
    mapping_functions, params, celldims, t, nhalo, cs, bt = _read_MappedGrid_payload(
      filename, Val(3)
    )
    mesh = MappedGrid(
      mapping_functions.x1,
      mapping_functions.x2,
      mapping_functions.x3,
      params,
      celldims,
      nhalo;
      t=t,
      compute_metrics=compute_metrics,
      cache_mode=cache_mode,
      coordinate_system=cs,
      basis=bt,
    )
  elseif grid_type == "DiscreteGrid1D"
    x = read_CurvilinearGrid1D(filename)
    mesh = DiscreteGrid(x, nhalo; compute_metrics=compute_metrics, cache_mode=cache_mode)
  elseif grid_type == "DiscreteGrid2D"
    x, y = read_CurvilinearGrid2D(filename)
    mesh = DiscreteGrid(x, y, nhalo; compute_metrics=compute_metrics, cache_mode=cache_mode)
  elseif grid_type == "DiscreteGrid3D"
    x, y, z = read_CurvilinearGrid3D(filename)
    mesh = DiscreteGrid(
      x, y, z, nhalo; compute_metrics=compute_metrics, cache_mode=cache_mode
    )
  elseif grid_type == "UniformGrid1D"
    x = read_UniformGrid1D(filename)
    mesh = UniformGrid1D(x, _read_discretization_scheme(filename))
  elseif grid_type == "UniformGrid2D"
    x0, x1, y0, y1, ∂x = read_UniformGrid2D(filename)
    mesh = UniformGrid2D((x0, y0), (x1, y1), ∂x, _read_discretization_scheme(filename))
  elseif grid_type == "UniformGrid3D"
    x0, x1, y0, y1, z0, z1, ∂x = read_UniformGrid3D(filename)
    mesh = UniformGrid3D(
      (x0, y0, z0), (x1, y1, z1), ∂x, _read_discretization_scheme(filename)
    )
  elseif grid_type == "RectilinearGrid2D"
    x, y = read_RectilinearGrid2D(filename)
    mesh = RectilinearGrid2D(x, y, _read_discretization_scheme(filename))
  elseif grid_type == "RectilinearGrid3D"
    x, y, z = read_RectilinearGrid3D(filename)
    mesh = RectilinearGrid3D(x, y, z, _read_discretization_scheme(filename))
  elseif grid_type == "CylindricalGrid1D"
    x, snap_to_axis = read_CylindricalGrid1D(filename)
    mesh = CylindricalGrid1D(x, :meg6, snap_to_axis)
  elseif grid_type == "AxisymmetricGrid2D"
    x, y, snap_to_axis, rotational_axis = read_AxisymmetricGrid2D(filename)
    mesh = AxisymmetricGrid2D(x, y, :meg6, snap_to_axis, rotational_axis)
  elseif grid_type == "SphericalGrid1D"
    x, snap_to_axis = read_SphericalGrid1D(filename)
    mesh = SphericalGrid1D(x, :meg6, snap_to_axis)
  elseif grid_type == "SphericalGrid3D"
    r, θ, ϕ, nhalo = read_SphericalGrid3D(filename)
    mesh = SphericalGrid3D(r, θ, ϕ, nhalo)
  else
    throw(ArgumentError("Unsupported serialized grid type `$grid_type`."))
  end

  return mesh
end

function _read_MappedGrid_payload(filename::String, ::Val{N}) where {N}
  grid_file = h5open(filename, "r")
  mapping_bytes = read(grid_file, "mapping_functions_serialized")
  params_bytes = read(grid_file, "mapping_params_serialized")
  celldims = Tuple(Int.(read(grid_file, "celldims")))
  t = read(grid_file, "time")
  nhalo = Int(read(grid_file, "nhalo"))
  coordinate_tag = read(grid_file, "coordinate_system_trait")
  basis_tag = read(grid_file, "basis_trait")
  close(grid_file)

  if length(celldims) != N
    throw(
      ArgumentError(
        "Serialized mapped-grid celldims length $(length(celldims)) does not match expected dimension $N.",
      ),
    )
  end

  mapping_functions = _deserialize_payload(Vector{UInt8}(mapping_bytes))
  params = _deserialize_payload(Vector{UInt8}(params_bytes))
  coordinate_system = _coordinate_system_from_tag(String(coordinate_tag))
  basis = _basis_trait_from_tag(String(basis_tag))

  return mapping_functions,
  params, NTuple{N,Int}(celldims), t, nhalo, coordinate_system,
  basis
end

function read_CurvilinearGrid1D(filename::String)
  grid_file = h5open(filename, "r")

  x = read(grid_file, "x")

  close(grid_file)

  return collect(x)
end

function read_CurvilinearGrid2D(filename::String)
  grid_file = h5open(filename, "r")

  x = read(grid_file, "x")
  y = read(grid_file, "y")

  close(grid_file)

  return collect(x'), collect(y')
end

function read_CurvilinearGrid3D(filename::String)
  grid_file = h5open(filename, "r")

  x = read(grid_file, "x")
  y = read(grid_file, "y")
  z = read(grid_file, "z")

  close(grid_file)

  return collect(x), collect(y), collect(z)
end

function read_UniformGrid1D(filename::String)
  grid_file = h5open(filename, "r")

  x = read(grid_file, "x")

  close(grid_file)

  return collect(x)
end

function read_UniformGrid2D(filename::String)
  grid_file = h5open(filename, "r")

  x = read(grid_file, "x")
  y = read(grid_file, "y")
  x0 = minimum(x)
  x1 = maximum(x)
  y0 = minimum(y)
  y1 = maximum(y)
  ∂x = collect(x')[2, 1] - collect(x')[1, 1]

  close(grid_file)

  return x0, x1, y0, y1, ∂x
end

function read_UniformGrid3D(filename::String)
  grid_file = h5open(filename, "r")

  x = read(grid_file, "x")
  y = read(grid_file, "y")
  z = read(grid_file, "z")
  x0 = minimum(x)
  x1 = maximum(x)
  y0 = minimum(y)
  y1 = maximum(y)
  z0 = minimum(z)
  z1 = maximum(z)
  ∂x = x[2, 1, 1] - x[1, 1, 1]

  close(grid_file)

  return x0, x1, y0, y1, z0, z1, ∂x
end

function read_RectilinearGrid2D(filename::String)
  grid_file = h5open(filename, "r")

  x = read(grid_file, "x")
  y = read(grid_file, "y")

  close(grid_file)

  xvec = collect(x')[:, 1]
  yvec = collect(y')[1, :]

  return xvec, yvec
end

function read_RectilinearGrid3D(filename::String)
  grid_file = h5open(filename, "r")

  x = read(grid_file, "x")
  y = read(grid_file, "y")
  z = read(grid_file, "z")

  close(grid_file)

  xvec = x[:, 1, 1]
  yvec = y[1, :, 1]
  zvec = z[1, 1, :]

  return xvec, yvec, zvec
end

function read_CylindricalGrid1D(filename::String)
  grid_file = h5open(filename, "r")

  x = read(grid_file, "x")
  snap_to_axis = read(grid_file, "snap_to_axis")

  close(grid_file)

  return collect(x), snap_to_axis
end

function read_AxisymmetricGrid2D(filename::String)
  grid_file = h5open(filename, "r")

  x = read(grid_file, "x")
  y = read(grid_file, "y")
  snap_to_axis = read(grid_file, "snap_to_axis")
  rotational_axis = Symbol(read(grid_file, "rotational_axis"))

  close(grid_file)

  return collect(x'), collect(y'), snap_to_axis, rotational_axis
end

function read_SphericalGrid1D(filename::String)
  grid_file = h5open(filename, "r")

  x = read(grid_file, "x")
  snap_to_axis = read(grid_file, "snap_to_axis")

  close(grid_file)

  return collect(x), snap_to_axis
end

function read_SphericalGrid3D(filename::String)
  grid_file = h5open(filename, "r")

  r = read(grid_file, "r")
  θ = read(grid_file, "theta")
  ϕ = read(grid_file, "phi")
  nhalo = read(grid_file, "nhalo")

  close(grid_file)

  return collect(r), collect(θ), collect(ϕ), nhalo
end

end
