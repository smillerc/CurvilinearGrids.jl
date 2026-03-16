using Unitful
using KernelAbstractions

mapped2d_x(t, ξ, η, p) = ξ + p.αx * sin(η / p.Lη)
mapped2d_y(t, ξ, η, p) = η + p.αy * cos(ξ / p.Lξ)

function initialize_mesh_Discrete1D()
  x0 = 0.0u"m"
  x1 = 1.0u"m"

  Nx = 50

  x = LinRange(ustrip.(u"cm", x0), ustrip.(u"cm", x1), Nx)

  return DiscreteGrid(x, 5; compute_metrics=false, cache_mode=:off)
end

function initialize_mesh_Discrete2D()
  x0 = 0.0u"m"
  x1 = 1.0u"m"
  y0 = 0.0u"m"
  y1 = 0.5u"m"

  Nx = 50
  Ny = 50

  x = LinRange(ustrip.(u"cm", x0), ustrip.(u"cm", x1), Nx)
  y = LinRange(ustrip.(u"cm", y0), ustrip.(u"cm", y1), Ny)

  x2d = [x[i] for i in eachindex(x), j in eachindex(y)]
  y2d = [y[j] for i in eachindex(x), j in eachindex(y)]

  return DiscreteGrid(x2d, y2d, 5; compute_metrics=false, cache_mode=:off)
end

function initialize_mesh_Discrete3D()
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

  x3d = [x[i] for i in eachindex(x), j in eachindex(y), k in eachindex(z)]
  y3d = [y[j] for i in eachindex(x), j in eachindex(y), k in eachindex(z)]
  z3d = [z[k] for i in eachindex(x), j in eachindex(y), k in eachindex(z)]

  return DiscreteGrid(x3d, y3d, z3d, 5; compute_metrics=false, cache_mode=:off)
  #return SphericalGrid1D(x, :meg6, snap_to_axis)
end

function initialize_mesh_Uniform1D()
  return UniformGrid1D((0.0, 100.0), 49, :meg6)
end

function initialize_mesh_Uniform2D()
  return UniformGrid2D((0.0, 0.0), (100.0, 50.0), 2.0, :meg6)
end

function initialize_mesh_Uniform3D()
  return UniformGrid3D((0.0, 0.0, 0.0), (100.0, 50.0, 200.0), 5.0, :meg6)
end

function initialize_mesh_Rectilinear2D()
  x = [0.0, 4.0, 9.5, 17.0, 27.0, 40.0, 56.0, 75.0, 100.0]
  y = [0.0, 3.0, 7.5, 14.0, 23.0, 35.0, 50.0]
  return RectilinearGrid2D(x, y, :meg6)
end

function initialize_mesh_Rectilinear3D()
  x = [0.0, 4.0, 9.5, 17.0, 27.0, 40.0, 56.0, 75.0, 100.0]
  y = [0.0, 3.0, 7.5, 14.0, 23.0, 35.0, 50.0]
  z = [0.0, 6.0, 15.0, 29.0, 50.0, 79.0, 116.0, 160.0, 200.0]
  return RectilinearGrid3D(
    x, y, z, :meg6
  )
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

function initialize_mesh_Spherical3D()
  r = collect(range(100.0, 200.0; length=8))
  θ = collect(range(0.3, 1.1; length=7))
  ϕ = collect(range(-0.2, 0.4; length=6))
  return SphericalGrid3D(r, θ, ϕ, 1, CPU())
end

function initialize_mesh_Mapped2D()
  params = (; αx=0.1, αy=0.15, Lξ=20.0, Lη=16.0)
  return MappedGrid(
    mapped2d_x,
    mapped2d_y,
    params,
    (20, 16),
    5;
    compute_metrics=false,
    cache_mode=:off,
    coordinate_system=CurvilinearCS(),
    basis=CartesianBasis(),
  )
end

function test_write(mesh)
  grid_fn = "$(@__DIR__)/grid.h5"
  units = u"cm"

  write_coordinates(mesh, grid_fn, units)
end

function test_read(; kwargs...)
  grid_fn = "$(@__DIR__)/grid.h5"

  mesh = read_coordinates(grid_fn; kwargs...)

  rm(grid_fn; force=true)

  return mesh
end

@testset "Read/Write DiscreteGrid to .h5" begin
  mesh1D_0 = initialize_mesh_Discrete1D()
  mesh2D_0 = initialize_mesh_Discrete2D()
  mesh3D_0 = initialize_mesh_Discrete3D()

  test_write(mesh1D_0)
  mesh1D = test_read(; compute_metrics=false, cache_mode=:off)

  test_write(mesh2D_0)
  mesh2D = test_read(; compute_metrics=false, cache_mode=:off)

  test_write(mesh3D_0)
  mesh3D = test_read(; compute_metrics=false, cache_mode=:off)

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

@testset "Read/Write SphericalGrid3D to .h5" begin
  mesh3D_0 = initialize_mesh_Spherical3D()

  test_write(mesh3D_0)
  mesh3D = test_read()

  r0, θ0, ϕ0 = coords(mesh3D_0)
  r, θ, ϕ = coords(mesh3D)

  @test r0 == r
  @test θ0 == θ
  @test ϕ0 == ϕ
  @test mesh3D.nhalo == mesh3D_0.nhalo
end

@testset "Read/Write MappedGrid to .h5" begin
  mesh2D_0 = initialize_mesh_Mapped2D()

  test_write(mesh2D_0)
  mesh2D = test_read(; compute_metrics=false, cache_mode=:off)

  x2_0, y2_0 = coords(mesh2D_0)
  x2, y2 = coords(mesh2D)
  @test x2_0 == x2
  @test y2_0 == y2

  ξη = (7.25, 5.5)
  @test coord(mesh2D_0, ξη) ≈ coord(mesh2D, ξη)
  @test typeof(coordinate_system(mesh2D)) == typeof(coordinate_system(mesh2D_0))
  @test typeof(basis_trait(mesh2D)) == typeof(basis_trait(mesh2D_0))
end
