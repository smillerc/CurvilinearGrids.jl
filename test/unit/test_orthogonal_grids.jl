@testset "Orthogonal reduced grids" begin
  backend = CPU()

  @testset "Cartesian 1D" begin
    x = [0.0, 1.0, 3.0]
    grid = CartesianOrthogonalGrid1D(x, 0, backend)

    domain = grid.iterators.cell.domain
    @test grid.centroid_coordinates.x[domain][1] == 0.5
    @test grid.centroid_coordinates.x[domain][2] == 2.0
    @test grid.cell_volumes[domain] == [1.0, 2.0]
    @test grid.face_areas.i₊½[domain] == [1.0, 1.0]
  end

  @testset "Cylindrical 1D" begin
    r = [1.0, 2.0, 4.0]
    grid = CylindricalOrthogonalGrid1D(r, 0, backend)

    domain = grid.iterators.cell.domain
    @test isapprox(grid.centroid_coordinates.r[domain][1], 14 / 9; atol=1e-12)
    @test isapprox(grid.cell_volumes[domain], [3π, 12π])
    @test isapprox(grid.face_areas.i₊½[domain], [4π, 8π])
  end

  @testset "Spherical 1D" begin
    r = [1.0, 3.0]
    grid = SphericalOrthogonalGrid1D(r, 0, backend)

    domain = grid.iterators.cell.domain
    @test isapprox(grid.centroid_coordinates.r[domain][1], 30 / 13; atol=1e-12)
    @test isapprox(grid.cell_volumes[domain][1], (4 / 3) * π * (27 - 1))
    @test isapprox(grid.face_areas.i₊½[domain][1], 36π)
  end

  @testset "Axisymmetric RZ 2D" begin
    r = [1.0, 2.0]
    z = [0.0, 1.0, 2.0]
    grid = AxisymmetricOrthogonalGrid2D(r, z, 0, backend)

    domain = grid.iterators.cell.domain
    @test isapprox(
      grid.centroid_coordinates.r[domain.indices[1][1]], 14 / 9; atol=1e-12
    )
    @test grid.centroid_coordinates.z[domain.indices[2]] == [0.5, 1.5]
    @test isapprox(grid.cell_volumes[domain], [3π 3π])
    @test isapprox(grid.face_areas.i₊½[domain], [4π 4π])
    @test isapprox(grid.face_areas.j₊½[domain], [3π 3π])
  end
end
