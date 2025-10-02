module H5Output

using HDF5, Unitful
using ..GridTypes

export write_coordinates, read_coordinates

### Write h5 functions

""" Write 1D Curvilinear grid information to .h5 format """
function write_coordinates(
  mesh::CurvilinearGrid1D, filename::String, units::Unitful.FreeUnits
)
  dim_type = dimension(units)
  if dim_type != Unitful.𝐋
    @error("Passed in units should be a Unitful.𝐋. $units is a Unitful.$(dim_type)")
  end

  x = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)
  end
end

""" Write 2D Curvilinear grid information to .h5 format """
function write_coordinates(
  mesh::CurvilinearGrid2D, filename::String, units::Unitful.FreeUnits
)
  dim_type = dimension(units)
  if dim_type != Unitful.𝐋
    @error("Passed in units should be a Unitful.𝐋. $units is a Unitful.$(dim_type)")
  end

  x, y = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x'))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y'))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)
  end
end

""" Write 3D Curvilinear grid information to .h5 format """
function write_coordinates(
  mesh::CurvilinearGrid3D, filename::String, units::Unitful.FreeUnits
)
  dim_type = dimension(units)
  if dim_type != Unitful.𝐋
    @error("Passed in units should be a Unitful.𝐋. $units is a Unitful.$(dim_type)")
  end

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
  end
end

""" Write 1D Uniform grid information to .h5 format """
function write_coordinates(mesh::UniformGrid1D, filename::String, units::Unitful.FreeUnits)
  dim_type = dimension(units)
  if dim_type != Unitful.𝐋
    @error("Passed in units should be a Unitful.𝐋. $units is a Unitful.$(dim_type)")
  end

  x = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)
  end
end

""" Write 2D Uniform grid information to .h5 format """
function write_coordinates(mesh::UniformGrid2D, filename::String, units::Unitful.FreeUnits)
  dim_type = dimension(units)
  if dim_type != Unitful.𝐋
    @error("Passed in units should be a Unitful.𝐋. $units is a Unitful.$(dim_type)")
  end

  x, y = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x'))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y'))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)
  end
end

""" Write 3D Uniform grid information to .h5 format """
function write_coordinates(mesh::UniformGrid3D, filename::String, units::Unitful.FreeUnits)
  dim_type = dimension(units)
  if dim_type != Unitful.𝐋
    @error("Passed in units should be a Unitful.𝐋. $units is a Unitful.$(dim_type)")
  end

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
  end
end

""" Write 2D Rectilinear grid information to .h5 format """
function write_coordinates(
  mesh::RectilinearGrid2D, filename::String, units::Unitful.FreeUnits
)
  dim_type = dimension(units)
  if dim_type != Unitful.𝐋
    @error("Passed in units should be a Unitful.𝐋. $units is a Unitful.$(dim_type)")
  end

  x, y = coords(mesh)
  mesh_type = String(nameof(typeof(mesh)))

  # Write mesh
  h5open(filename, "w") do file
    h5write(filename, "x", collect(x'))
    h5writeattr(filename, "x", Dict("Units" => string(units)))

    h5write(filename, "y", collect(y'))
    h5writeattr(filename, "y", Dict("Units" => string(units)))

    h5write(filename, "grid_type", mesh_type)
  end
end

""" Write 3D Rectilinear grid information to .h5 format """
function write_coordinates(
  mesh::RectilinearGrid3D, filename::String, units::Unitful.FreeUnits
)
  dim_type = dimension(units)
  if dim_type != Unitful.𝐋
    @error("Passed in units should be a Unitful.𝐋. $units is a Unitful.$(dim_type)")
  end

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
  end
end

""" Write 1D Cylindrical grid information to .h5 format """
function write_coordinates(
  mesh::CylindricalGrid1D, filename::String, units::Unitful.FreeUnits
)
  dim_type = dimension(units)
  if dim_type != Unitful.𝐋
    @error("Passed in units should be a Unitful.𝐋. $units is a Unitful.$(dim_type)")
  end

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

""" Write 2D Axisymmetric grid information to .h5 format """
function write_coordinates(
  mesh::AxisymmetricGrid2D, filename::String, units::Unitful.FreeUnits
)
  dim_type = dimension(units)
  if dim_type != Unitful.𝐋
    @error("Passed in units should be a Unitful.𝐋. $units is a Unitful.$(dim_type)")
  end

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

""" Write 1D Spherical grid information to .h5 format """
function write_coordinates(
  mesh::SphericalGrid1D, filename::String, units::Unitful.FreeUnits
)
  dim_type = dimension(units)
  if dim_type != Unitful.𝐋
    @error("Passed in units should be a Unitful.𝐋. $units is a Unitful.$(dim_type)")
  end

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

""" Read grid info from .h5 format and return mesh """
function read_coordinates(filename::String; discretization_scheme=:meg6)
  grid_file = h5open(filename, "r")
  grid_type = read(grid_file, "grid_type")
  close(grid_file)

  if grid_type == "CurvilinearGrid1D"
    x = read_CurvilinearGrid1D(filename)
    mesh = CurvilinearGrid1D(x, discretization_scheme)
  elseif grid_type == "CurvilinearGrid2D"
    x, y = read_CurvilinearGrid2D(filename)
    mesh = CurvilinearGrid2D(x, y, discretization_scheme)
  elseif grid_type == "CurvilinearGrid3D"
    x, y, z = read_CurvilinearGrid3D(filename)
    mesh = CurvilinearGrid3D(x, y, z, discretization_scheme)
  elseif grid_type == "UniformGrid1D"
    x = read_UniformGrid1D(filename)
    mesh = UniformGrid1D(x, discretization_scheme)
  elseif grid_type == "UniformGrid2D"
    x0, x1, y0, y1, ∂x = read_UniformGrid2D(filename)
    mesh = UniformGrid2D((x0, y0), (x1, y1), ∂x, discretization_scheme)
  elseif grid_type == "UniformGrid3D"
    x0, x1, y0, y1, z0, z1, ∂x = read_UniformGrid3D(filename)
    mesh = UniformGrid3D((x0, y0, z0), (x1, y1, z1), ∂x, discretization_scheme)
  elseif grid_type == "RectilinearGrid2D"
    x, y = read_RectilinearGrid2D(filename)
    mesh = RectilinearGrid2D(x, y, discretization_scheme)
  elseif grid_type == "RectilinearGrid3D"
    x, y, z = read_RectilinearGrid3D(filename)
    mesh = RectilinearGrid3D(x, y, z, discretization_scheme)
  elseif grid_type == "CylindricalGrid1D"
    x, snap_to_axis = read_CylindricalGrid1D(filename)
    mesh = CylindricalGrid1D(x, discretization_scheme, snap_to_axis)
  elseif grid_type == "AxisymmetricGrid2D"
    x, y, snap_to_axis, rotational_axis = read_AxisymmetricGrid2D(filename)
    mesh = AxisymmetricGrid2D(x, y, discretization_scheme, snap_to_axis, rotational_axis)
  elseif grid_type == "SphericalGrid1D"
    x, snap_to_axis = read_SphericalGrid1D(filename)
    mesh = SphericalGrid1D(x, discretization_scheme, snap_to_axis)
  end

  return mesh
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

end