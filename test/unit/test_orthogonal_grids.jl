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

  @testset "Halo coordinate geometry" begin
    nhalo = 2

    cartesian = CartesianOrthogonalGrid1D([0.0, 1.0, 3.0], nhalo, backend)
    xnodes = Array(cartesian.node_coordinates[1])
    xcenters = Array(cartesian.centroid_coordinates[1])
    @test xnodes == [-2.0, -1.0, 0.0, 1.0, 3.0, 5.0, 7.0]
    @test all(isfinite, xcenters)
    @test xcenters == (xnodes[1:(end - 1)] .+ xnodes[2:end]) ./ 2

    radial_nodes = [0.0, 0.1, 0.23, 0.5]
    spherical = SphericalOrthogonalGrid1D(radial_nodes, nhalo, backend)
    rnodes = Array(spherical.node_coordinates[1])
    rcenters = Array(spherical.centroid_coordinates[1])
    @test rnodes[1:2] == [-0.23, -0.1]
    @test rcenters[1] == -rcenters[4]
    @test rcenters[2] == -rcenters[3]
    @test all(>(0), diff(rcenters))

    spherical_2d = SphericalOrthogonalGrid2D(
      radial_nodes, [0.0, π / 4, π / 2, 3π / 4, π], nhalo, backend
    )
    @test Array(spherical_2d.node_coordinates[1])[1:2] == [-0.23, -0.1]
    rcenters_2d = Array(spherical_2d.centroid_coordinates[1])
    @test rcenters_2d[1] == -rcenters_2d[4]
    @test rcenters_2d[2] == -rcenters_2d[3]

    spherical_3d = SphericalGrid3D(
      radial_nodes, [0.0, π / 4, π / 2, 3π / 4, π], [0.0, 0.5, 1.0, 1.5], nhalo, backend
    )
    @test Array(spherical_3d.node_coordinates[1])[1:2] == [-0.23, -0.1]
    rcenters_3d = Array(spherical_3d.centroid_coordinates[1])
    @test rcenters_3d[1] == -rcenters_3d[4]
    @test rcenters_3d[2] == -rcenters_3d[3]

    away_from_origin = SphericalOrthogonalGrid1D([1.0, 1.2, 1.7], nhalo, backend)
    @test Array(away_from_origin.node_coordinates[1])[1:2] ≈ [0.6, 0.8]

    supplied_nodes = [-0.7, -0.2, 0.0, 0.3, 0.9, 1.8, 2.8]
    supplied = CartesianOrthogonalGrid1D(
      supplied_nodes, nhalo, backend; halo_coords_included=true
    )
    @test Array(supplied.node_coordinates[1]) == supplied_nodes
    supplied_centers = Array(supplied.centroid_coordinates[1])
    @test supplied_centers == (supplied_nodes[1:(end - 1)] .+ supplied_nodes[2:end]) ./ 2

    @test_throws ArgumentError SphericalOrthogonalGrid1D([0.0, 0.1], nhalo, backend)
  end

  @testset "Unwrapped spherical polar centroids" begin
    grid = SphericalOrthogonalGrid2D([1.0, 2.0, 3.0], [0.0, 0.4, 1.5, π], 1, backend)
    θcenters = Array(grid.centroid_coordinates[2])
    @test θcenters[1] < 0
    @test θcenters[end] > π
    @test all(>(0), diff(θcenters))

    grid_3d = SphericalGrid3D(
      [1.0, 2.0, 3.0], [0.0, 0.4, 1.5, π], [0.0, 0.5, 1.0], 1, backend
    )
    θcenters_3d = Array(grid_3d.centroid_coordinates[2])
    @test θcenters_3d[1] < 0
    @test θcenters_3d[end] > π
    @test all(>(0), diff(θcenters_3d))
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

  @testset "Spherical (r, θ) 2D" begin
    r = [1.0, 2.0]
    θ = [0.0, π / 3, 2π / 3]
    grid = SphericalOrthogonalGrid2D(r, θ, 0, backend)

    domain = grid.iterators.cell.domain
    @test isapprox(grid.centroid_coordinates[1][domain.indices[1][1]], 45 / 28; atol=1e-12)
    @test grid.centroid_coordinates[2][domain.indices[2]] ≈ [acos(3 / 4), π / 2]
    @test isapprox(grid.cell_volumes[domain], [7π / 3 14π / 3])
    @test isapprox(grid.face_areas[1][domain], [4π 8π])
    @test isapprox(grid.face_areas[2][domain], [3π * sqrt(3) / 2 3π * sqrt(3) / 2])
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

    g_sph2 = SphericalOrthogonalGrid2D([1.0, 2.0], [0.0, π / 2], 0, backend)
    x2, z2 = cartesian_coordinates(g_sph2)
    @test size(x2) == (2, 2)
    @test size(z2) == (2, 2)
    @test isapprox(x2[2, 2], 2.0; atol=1e-12)
    @test isapprox(z2[2, 1], 2.0; atol=1e-12)
  end
end
