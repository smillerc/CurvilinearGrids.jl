
@testitem "3D Mesh - Rectangular Mesh" begin
  include("common.jl")

  function rect_grid(nx, ny, nz)
    x0, x1 = (0, 2)
    y0, y1 = (1, 3)
    z0, z1 = (-1, 2)

    _x(ξ) = @. x0 + (x1 - x0) * ((ξ - 1) / (nx - 1))
    _y(η) = @. y0 + (y1 - y0) * ((η - 1) / (ny - 1))
    _z(ζ) = @. z0 + (z1 - z0) * ((ζ - 1) / (nz - 1))

    x(ξ, η, ζ) = _x(ξ)
    y(ξ, η, ζ) = _y(η)
    z(ξ, η, ζ) = _z(ζ)

    x((ξ, η, ζ)) = _x(ξ)
    y((ξ, η, ζ)) = _y(η)
    z((ξ, η, ζ)) = _z(ζ)

    return (x, y, z)
  end

  ni, nj, nk = (5, 9, 13)
  nhalo = 0
  x, y, z = rect_grid(ni, nj, nk)
  mesh = CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)

  metrics(mesh, (2, 3, 4)) == (
    ξx=2.0,
    ξy=0.0,
    ξz=0.0,
    ηx=0.0,
    ηy=4.0,
    ηz=0.0,
    ζx=0.0,
    ζy=0.0,
    ζz=4.0,
    ξt=0.0,
    ηt=0.0,
    ζt=0.0,
  )

  conservative_metrics(mesh, (2, 3, 4)) == (
    ξ̂x=0.0625,
    ξ̂y=0.0,
    ξ̂z=0.0,
    η̂x=0.0,
    η̂y=0.125,
    η̂z=0.0,
    ζ̂x=0.0,
    ζ̂y=0.0,
    ζ̂z=0.125,
    ξt=0.0,
    ηt=0.0,
    ζt=0.0,
  )

  bm1 = @benchmark metrics($mesh, (2, 3, 4))
  @test bm1.allocs == 0

  @test jacobian_matrix(mesh, (2, 3, 4)) == @SMatrix [
    0.5 0.0 0.0
    0.0 0.25 0.0
    0.0 0.0 0.25
  ]
  bm2 = @benchmark jacobian_matrix($mesh, (2, 3, 4))
  @test bm2.allocs == 0

  cell_volume = 0.5 * 0.25 * 0.25
  @test jacobian(mesh, (2, 3, 4)) == cell_volume

  bm3 = @benchmark jacobian($mesh, (2, 3, 4))
  @test bm3.allocs == 0

  ilo, ihi, jlo, jhi, klo, khi = mesh.limits
  @test ilo = jlo = klo == 1
  @test ihi == 4
  @test jhi == 8
  @test khi == 12

  # for k in klo:khi
  #   for j in jlo:jhi
  #     for i in ilo:ihi
  #       @test metrics(mesh, (i, j, k)) == m
  #     end
  #   end
  # end

  @test coord(mesh, 1, 1, 1) == [0, 1, -1]
  @test coord(mesh, 2.5, 2.5, 2.5) == [0.75, 1.375, -0.625]
  @test coord(mesh, 2.5, 2.5, 2.5) == centroid(mesh, 2, 2, 2)

  xyz_coords = coords(mesh)
  centroid_coords = centroids(mesh)
  @test size(centroid_coords) == (3, 4, 8, 12)
  @test size(xyz_coords) == (3, 5, 9, 13)
end

@testitem "3D Mesh - Sphere Sector" begin
  include("common.jl")

  function sphere_grid(nr, ntheta, nphi)
    r0, r1 = (1, 3)
    (θ0, θ1) = deg2rad.((35, 180 - 35))
    (ϕ0, ϕ1) = deg2rad.((45, 360 - 45))

    r(ξ) = r0 + (r1 - r0) * ((ξ - 1) / (nr - 1))
    θ(η) = θ0 + (θ1 - θ0) * ((η - 1) / (ntheta - 1))
    ϕ(ζ) = ϕ0 + (ϕ1 - ϕ0) * ((ζ - 1) / (nphi - 1))

    x(ξ, η, ζ) = r(ξ) * sin(θ(η)) * cos(ϕ(ζ))
    y(ξ, η, ζ) = r(ξ) * sin(θ(η)) * sin(ϕ(ζ))
    z(ξ, η, ζ) = r(ξ) * cos(θ(η))

    x((ξ, η, ζ)) = x(ξ, η, ζ)
    y((ξ, η, ζ)) = y(ξ, η, ζ)
    z((ξ, η, ζ)) = z(ξ, η, ζ)

    return (x, y, z)
  end

  ni, nj, nk = (5, 9, 11)
  nhalo = 0
  x, y, z = sphere_grid(ni, nj, nk)
  @test_nowarn CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)
end