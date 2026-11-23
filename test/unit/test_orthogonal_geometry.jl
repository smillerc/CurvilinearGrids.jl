using Adapt: adapt

_requires_orthogonal_geometry(::OrthogonalGeometry) = true

@testset "Orthogonal geometry dispatch" begin
    no_cache = (; compute_metrics = false, cache_mode = :off)
    cart = (; coordinate_system = CartesianCS())
    cyl = (; coordinate_system = CylindricalCS())
    sph = (; coordinate_system = SphericalCS(), basis = SphericalBasis())
    axi_x = (; coordinate_system = AxisymmetricCS{:x}())
    axi_y = (; coordinate_system = AxisymmetricCS{:y}())

    for (dims, options) in (
            ((4,), cart), ((4, 5), cart), ((4, 5, 6), cart),
            ((4,), cyl), ((4, 5), cyl), ((4, 5, 6), cyl),
            ((4, 5), axi_x), ((4, 5), axi_y),
            ((4,), sph), ((4, 5), sph), ((4, 5, 6), sph),
        )
        @testset "$(options.coordinate_system), $(length(dims))D" begin
            n = length(dims)
            axis_maps = ntuple(d -> (t, ξ, p) -> p.origin[d] + (ξ - 1) / 16, n)
            full_maps = ntuple(d -> (t, args...) -> axis_maps[d](t, args[d], last(args)), n)
            p = (; origin = ntuple(_ -> 1.0, n))
            nodes = ntuple(d -> [axis_maps[d](0.0, i, p) for i in 1:(dims[d] + 1)], n)
            arrays = ntuple(d -> [nodes[d][I[d]] for I in CartesianIndices(dims .+ 1)], n)

            mapped = MappedGrid(OrthogonalMapping(axis_maps...), p, dims, 1; options..., no_cache...)
            full = MappedGrid(full_maps..., p, dims, 1; options..., no_cache...)
            discrete = DiscreteGrid(nodes..., 1; options..., no_cache...)
            dense = DiscreteGrid(arrays..., 1; options..., no_cache...)

            for grid in (mapped, discrete)
                @test coordinate_system(grid) === options.coordinate_system
                @test (@inferred orthogonal_geometry(grid)) === OrthogonalGeometry()
                @test orthogonal_geometry(typeof(grid)) === OrthogonalGeometry()
                @test _requires_orthogonal_geometry(orthogonal_geometry(grid))
            end
            expected = n == 1 ? OrthogonalGeometry() : GeneralGeometry()
            @test orthogonal_geometry(full) === expected
            @test orthogonal_geometry(dense) === expected
            @test mapped.node_coordinates == full.node_coordinates == discrete.node_coordinates == dense.node_coordinates
            at = ntuple(_ -> 2.25, n)
            @test coord(mapped, at) ≈ coord(discrete, at)
            @test jacobian_matrix(mapped, at) ≈ jacobian_matrix(full, at)
            @test jacobian_matrix(discrete, at) ≈ jacobian_matrix(dense, at)
        end
    end

    axis(t, ξ, p) = p.a * (ξ - 1) + t
    coupled(t, ξ, η, p) = ξ + p.a * η
    vertical(t, ξ, η, p) = η
    p = (; a = 0.25)
    grid = MappedGrid(OrthogonalMapping(axis, axis), p, (6, 7), 1)
    warped = MappedGrid(coupled, vertical, p, (6, 7), 1; no_cache...)
    @test coordinate_system(grid) === coordinate_system(warped) === CartesianCS()
    @test_throws MethodError _requires_orthogonal_geometry(orthogonal_geometry(warped))
    @test_throws ArgumentError OrthogonalMapping()

    update!(grid, 0.5, (; a = 0.5))
    @test orthogonal_geometry(grid) === OrthogonalGeometry()
    block = local_grid(grid, (3:5, 2:5))
    @test orthogonal_geometry(block) === OrthogonalGeometry()
    @test coord(block, (3.0, 2.0)) == SVector(1.5, 1.0)
    @test orthogonal_geometry(adapt(Array, block)) === OrthogonalGeometry()
    callbacks = reparameterize_mapping(grid.mapping_functions, (3.0, 2.0), (2.0, 1.0))
    rebuilt = MappedGrid(values(callbacks)..., (; a = 0.5), (3, 4), 1; t = 0.5, no_cache...)
    @test orthogonal_geometry(rebuilt) === OrthogonalGeometry()
    @test coord(rebuilt, (1.0, 1.0)) == coord(block, (3.0, 2.0))

    # Supplied halos use the input samples, including nonuniform halo spacings.
    x = [0.0, 1.0, 1.25, 2.0, 3.0, 6.0]
    y = [1.0, 2.0, 2.5, 3.5, 5.0, 8.0, 9.0]
    discrete = DiscreteGrid(x, y, 1; halo_coords_included = true)
    dense = DiscreteGrid([a for a in x, _ in y], [b for _ in x, b in y], 1; halo_coords_included = true)
    @test coordinate_system(discrete) === CartesianCS()
    @test cellsize(discrete) == (3, 4)
    @test discrete.node_coordinates == dense.node_coordinates
    @test discrete.node_coordinates[1][:, 1] == x
    @test discrete.node_coordinates[2][1, :] == y
    I = first(discrete.iterators.cell.domain)
    @test cellvolume(discrete, I) ≈ cellvolume(dense, I)
    @test face_flux_geometry(discrete, I, :ihi).metric_vector ≈ face_flux_geometry(dense, I, :ihi).metric_vector
    update!(discrete)
    @test orthogonal_geometry(discrete) === OrthogonalGeometry()
    @test orthogonal_geometry(adapt(Array, discrete)) === OrthogonalGeometry()
    @test orthogonal_geometry(dense) === GeneralGeometry()
    @test_throws ArgumentError DiscreteGrid(x[1:3], y, 1; halo_coords_included = true)

    backend = CPU()
    r = [1.0, 1.25, 2.0]
    θ = [0.5, 0.75, 1.25]
    for grid in (
            CartesianOrthogonalGrid1D(r, 1, backend),
            CartesianOrthogonalGrid2D(r, r, 1, backend),
            CartesianOrthogonalGrid3D(r, r, r, 1, backend),
            CylindricalOrthogonalGrid1D(r, 1, backend),
            AxisymmetricOrthogonalGrid2D(r, r, 1, backend; rotational_axis = :x),
            AxisymmetricOrthogonalGrid2D(r, r, 1, backend; rotational_axis = :y),
            SphericalOrthogonalGrid1D(r, 1, backend),
            SphericalOrthogonalGrid2D(r, θ, 1, backend),
            SphericalGrid3D(r, θ, θ, 1, backend),
        )
        @test orthogonal_geometry(grid) === OrthogonalGeometry()
    end
    @test !isdefined(CurvilinearGrids, :CurvilinearCS)
end
