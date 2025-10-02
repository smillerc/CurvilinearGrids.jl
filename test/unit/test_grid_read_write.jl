using Unitful

function initialize_mesh_Curvilinear1D()
  x0 = 0.0u"m"
  x1 = 1.0u"m"

  Nx = 50

  x = LinRange(ustrip.(u"cm", x0), ustrip.(u"cm", x1), Nx)

  return CurvilinearGrid1D(x, :meg6)
end

function initialize_mesh_Curvilinear2D()
  x0 = 0.0u"m"
  x1 = 1.0u"m"
  y0 = 0.0u"m"
  y1 = 0.5u"m"

  Nx = 50
  Ny = 50

  x = LinRange(ustrip.(u"cm", x0), ustrip.(u"cm", x1), Nx)
  y = LinRange(ustrip.(u"cm", y0), ustrip.(u"cm", y1), Ny)

  return CurvilinearGrid2D(x, y, :meg6)
end

function initialize_mesh_Curvilinear3D()
  x0 = 0.0u"m"
  x1 = 1.0u"m"
  y0 = 0.0u"m"
  y1 = 0.5u"m"
  z0 = 0.0u"m"
  z1 = 2.0u"m"

  Nx = 50
  Ny = 50
  Nz = 50
  ∂x = 100 / 51

  x = LinRange(ustrip.(u"cm", x0), ustrip.(u"cm", x1), Nx)
  y = LinRange(ustrip.(u"cm", y0), ustrip.(u"cm", y1), Ny)
  z = LinRange(ustrip.(u"cm", z0), ustrip.(u"cm", z1), Nz)

  snap_to_axis = true
  rotational_axis = :y

  #=
  return UniformGrid3D(
    (ustrip.(u"cm", x0), ustrip.(u"cm", y0), ustrip.(u"cm", z0)),
    (ustrip.(u"cm", x1), ustrip.(u"cm", y1), ustrip.(u"cm", z1)),
    ∂x,
    :meg6,
  )
    =#
  return CurvilinearGrid3D(x, y, z, :meg6)
  #return SphericalGrid1D(x, :meg6, snap_to_axis)
end

function initialize_mesh_Uniform1D()
  x0 = 0.0u"m"
  x1 = 1.0u"m"

  Nx = 50

  return UniformGrid1D((ustrip.(u"cm", x0), ustrip.(u"cm", x1)), Nx, :meg6)
end

function initialize_mesh_Uniform2D()
  x0 = 0.0u"m"
  x1 = 1.0u"m"
  y0 = 0.0u"m"
  y1 = 0.5u"m"

  ∂x = 100 / 51

  return UniformGrid2D(
    (ustrip.(u"cm", x0), ustrip.(u"cm", y0)),
    (ustrip.(u"cm", x1), ustrip.(u"cm", y1)),
    ∂x,
    :meg6,
  )
end

function initialize_mesh_Uniform3D()
  x0 = 0.0u"m"
  x1 = 1.0u"m"
  y0 = 0.0u"m"
  y1 = 0.5u"m"
  z0 = 0.0u"m"
  z1 = 2.0u"m"

  ∂x = 100 / 51

  return UniformGrid3D(
    (ustrip.(u"cm", x0), ustrip.(u"cm", y0), ustrip.(u"cm", z0)),
    (ustrip.(u"cm", x1), ustrip.(u"cm", y1), ustrip.(u"cm", z1)),
    ∂x,
    :meg6,
  )
end

function initialize_mesh_Rectilinear2D()
  x0 = 0.0u"m"
  x1 = 1.0u"m"
  y0 = 0.0u"m"
  y1 = 0.5u"m"

  Nx = 50
  Ny = 50

  x = LinRange(ustrip.(u"cm", x0), ustrip.(u"cm", x1), Nx)
  y = LinRange(ustrip.(u"cm", y0), ustrip.(u"cm", y1), Ny)

  return RectilinearGrid2D(x, y, :meg6)
end

function initialize_mesh_Rectilinear3D()
  x0 = 0.0u"m"
  x1 = 1.0u"m"
  y0 = 0.0u"m"
  y1 = 0.5u"m"
  z0 = 0.0u"m"
  z1 = 2.0u"m"

  Nx = 50
  Ny = 50
  Nz = 50

  x = LinRange(ustrip.(u"cm", x0), ustrip.(u"cm", x1), Nx)
  y = LinRange(ustrip.(u"cm", y0), ustrip.(u"cm", y1), Ny)
  z = LinRange(ustrip.(u"cm", z0), ustrip.(u"cm", z1), Nz)

  return RectilinearGrid3D(x, y, z, :meg6)
end

function initialize_mesh_Cylindrical1D()
  x0 = 0.0u"m"
  x1 = 1.0u"m"

  Nx = 50

  x = LinRange(ustrip.(u"cm", x0), ustrip.(u"cm", x1), Nx)

  snap_to_axis = true

  return CylindricalGrid1D(x, :meg6, snap_to_axis)
end

function initialize_mesh_Axisymmetric2D()
  x0 = 0.0u"m"
  x1 = 1.0u"m"
  y0 = 0.0u"m"
  y1 = 0.5u"m"

  snap_to_axis = true
  rotational_axis = :y

  return AxisymmetricGrid2D(
    (ustrip.(u"cm", x0), ustrip.(u"cm", y0)),
    (ustrip.(u"cm", x1), ustrip.(u"cm", y1)),
    (50, 50),
    :meg6,
    snap_to_axis,
    rotational_axis,
  )
end

function initialize_mesh_Spherical1D()
  x0 = 0.0u"m"
  x1 = 1.0u"m"

  Nx = 50

  x = LinRange(ustrip.(u"cm", x0), ustrip.(u"cm", x1), Nx)

  snap_to_axis = true

  return SphericalGrid1D(x, :meg6, snap_to_axis)
end

function test_write(mesh)
  grid_fn = "$(@__DIR__)/grid.h5"
  units = u"cm"

  write_coordinates(mesh, grid_fn, units)
end

function test_read()
  grid_fn = "$(@__DIR__)/grid.h5"

  mesh = read_coordinates(grid_fn)

  rm(grid_fn)

  return mesh
end

@testset "Read/Write CurvilinearGrid to .h5" begin
  mesh1D_0 = initialize_mesh_Curvilinear1D()
  mesh2D_0 = initialize_mesh_Curvilinear2D()
  mesh3D_0 = initialize_mesh_Curvilinear3D()

  test_write(mesh1D_0)
  mesh1D = test_read()

  test_write(mesh2D_0)
  mesh2D = test_read()

  test_write(mesh3D_0)
  mesh3D = test_read()

  x1_0 = coords(mesh1D_0)
  x1 = coords(mesh1D)

  x2_0, y2_0 = coords(mesh2D_0)
  x2, y2 = coords(mesh2D)

  x3_0, y3_0, z3_0 = coords(mesh3D_0)
  x3, y3, z3 = coords(mesh3D)

  @test x1_0 == x1

  @test x2_0 == x2
  @test y2_0 == y2

  @test x3_0 == x3
  @test y3_0 == y3
  @test z3_0 == z3
end

@testset "Read/Write UniformGrid to .h5" begin
  mesh1D_0 = initialize_mesh_Uniform1D()
  mesh2D_0 = initialize_mesh_Uniform2D()
  mesh3D_0 = initialize_mesh_Uniform3D()

  test_write(mesh1D_0)
  mesh1D = test_read()

  test_write(mesh2D_0)
  mesh2D = test_read()

  test_write(mesh3D_0)
  mesh3D = test_read()

  x1_0 = coords(mesh1D_0)
  x1 = coords(mesh1D)

  x2_0, y2_0 = coords(mesh2D_0)
  x2, y2 = coords(mesh2D)

  x3_0, y3_0, z3_0 = coords(mesh3D_0)
  x3, y3, z3 = coords(mesh3D)

  @test x1_0 == x1

  @test x2_0 == x2
  @test y2_0 == y2

  @test x3_0 == x3
  @test y3_0 == y3
  @test z3_0 == z3
end

@testset "Read/Write RectilinearGrid to .h5" begin
  mesh2D_0 = initialize_mesh_Rectilinear2D()
  mesh3D_0 = initialize_mesh_Rectilinear3D()

  test_write(mesh2D_0)
  mesh2D = test_read()

  test_write(mesh3D_0)
  mesh3D = test_read()

  x2_0, y2_0 = coords(mesh2D_0)
  x2, y2 = coords(mesh2D)

  x3_0, y3_0, z3_0 = coords(mesh3D_0)
  x3, y3, z3 = coords(mesh3D)

  @test x2_0 == x2
  @test y2_0 == y2

  @test x3_0 == x3
  @test y3_0 == y3
  @test z3_0 == z3
end

@testset "Read/Write CylindricalGrid to .h5" begin
  mesh1D_0 = initialize_mesh_Cylindrical1D()

  test_write(mesh1D_0)
  mesh1D = test_read()

  x1_0 = coords(mesh1D_0)
  x1 = coords(mesh1D)

  @test x1_0 == x1
end

@testset "Read/Write AxisymmetricGrid to .h5" begin
  mesh2D_0 = initialize_mesh_Axisymmetric2D()

  test_write(mesh2D_0)
  mesh2D = test_read()

  x2_0, y2_0 = coords(mesh2D_0)
  x2, y2 = coords(mesh2D)

  @test x2_0 == x2
  @test y2_0 == y2
end

@testset "Read/Write SphericalGrid to .h5" begin
  mesh1D_0 = initialize_mesh_Spherical1D()

  test_write(mesh1D_0)
  mesh1D = test_read()

  x1_0 = coords(mesh1D_0)
  x1 = coords(mesh1D)

  @test x1_0 == x1
end