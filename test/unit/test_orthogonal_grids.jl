@testset "Orthogonal reduced grids" begin
  backend = CPU()

  @testset "Cartesian 1D" begin
    x = [0.0, 1.0, 3.0]
    grid = CartesianOrthogonalGrid1D(x, 0, backend)

    domain = grid.iterators.cell.domain
    @test grid.centroid_coordinates[1][domain][1] == 0.5
    @test grid.centroid_coordinates[1][domain][2] == 2.0
    @test grid.cell_volumes[domain] == [1.0, 2.0]
    @test grid.face_areas[1][domain] == [1.0, 1.0]
  end

  @testset "Cylindrical 1D" begin
    r = [1.0, 2.0, 4.0]
    grid = CylindricalOrthogonalGrid1D(r, 0, backend)

    domain = grid.iterators.cell.domain
    @test isapprox(grid.centroid_coordinates[1][domain][1], 14 / 9; atol=1e-12)
    @test isapprox(grid.cell_volumes[domain], [3π, 12π])
    @test isapprox(grid.face_areas[1][domain], [4π, 8π])
  end

  @testset "Spherical 1D" begin
    r = [1.0, 3.0]
    grid = SphericalOrthogonalGrid1D(r, 0, backend)

    domain = grid.iterators.cell.domain
    @test isapprox(grid.centroid_coordinates[1][domain][1], 30 / 13; atol=1e-12)
    @test isapprox(grid.cell_volumes[domain][1], (4 / 3) * π * (27 - 1))
    @test isapprox(grid.face_areas[1][domain][1], 36π)
  end

  @testset "Axisymmetric RZ 2D" begin
    r = [1.0, 2.0]
    z = [0.0, 1.0, 2.0]
    grid = AxisymmetricOrthogonalGrid2D(r, z, 0, backend)

    domain = grid.iterators.cell.domain
    @test isapprox(grid.centroid_coordinates[1][domain.indices[1][1]], 14 / 9; atol=1e-12)
    @test grid.centroid_coordinates[2][domain.indices[2]] == [0.5, 1.5]
    @test isapprox(grid.cell_volumes[domain], [3π 3π])
    @test isapprox(grid.face_areas[1][domain], [4π 4π])
    @test isapprox(grid.face_areas[2][domain], [3π 3π])
  end

  @testset "Cartesian Coordinate Conversion" begin
    g_cart = CartesianOrthogonalGrid1D([0.0, 1.0, 2.0], 0, backend)
    @test cartesian_coordinates(g_cart) == coords(g_cart)

    g_axi = AxisymmetricOrthogonalGrid2D([1.0, 2.0], [0.0, 1.0], 0, backend)
    r_axi, z_axi = coords(g_axi)
    x_axi, y_axi = cartesian_coordinates(g_axi)
    @test x_axi == r_axi
    @test y_axi == z_axi

    r = [1.0, 2.0]
    θ = [π / 2, π / 2]
    ϕ = [0.0, π / 2]
    g_sph = SphericalGrid3D(r, θ, ϕ, 0, backend)
    x, y, z = cartesian_coordinates(g_sph)

    @test size(x) == (2, 2, 2)
    @test size(y) == (2, 2, 2)
    @test size(z) == (2, 2, 2)
    @test isapprox(x[1, 1, 1], 1.0; atol=1e-12)
    @test isapprox(y[1, 1, 2], 1.0; atol=1e-12)
    @test isapprox(x[2, 1, 2], 0.0; atol=1e-12)
    @test isapprox(y[2, 1, 1], 0.0; atol=1e-12)
    @test isapprox(z[1, 1, 1], 0.0; atol=1e-12)
  end
end
