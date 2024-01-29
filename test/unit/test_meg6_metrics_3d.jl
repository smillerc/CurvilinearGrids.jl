# @testset "3D Mesh - Sphere Sector" begin
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

  return (x, y, z)
end

ni, nj, nk = (5, 9, 11)
nhalo = 0
x, y, z = sphere_grid(ni, nj, nk)
mesh = CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)

#   @test_nowarn CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)
# end

metrics = MEG6Scheme(cellsize_withhalo(mesh));
domain = mesh.iterators.cell.domain
update_metrics!(metrics, mesh.x_coord, mesh.y_coord, mesh.z_coord, domain)