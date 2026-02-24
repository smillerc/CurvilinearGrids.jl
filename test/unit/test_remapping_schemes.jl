using LinearAlgebra
using Test
using CurvilinearGrids

@testset "RemappingSchemes conservative scalar remap" begin
  xsrc(t, xi, eta, p) = xi
  ysrc(t, xi, eta, p) = eta
  xdst(t, xi, eta, p) = xi + p.dx
  ydst(t, xi, eta, p) = eta + p.dy

  source_grid = MappedGrid(
    xsrc, ysrc, (;), (24, 18), 2; compute_metrics=false, cache_mode=:off
  )
  destination_grid = MappedGrid(
    xdst, ydst, (; dx=3.0, dy=-2.0), (24, 18), 2; compute_metrics=false, cache_mode=:off
  )

  cache = build_remap_cache(source_grid, destination_grid; quadrature_order=2)
  @test cache.sample_count == 4
  @test cache.source_shape == cellsize(source_grid)
  @test cache.destination_shape == cellsize(destination_grid)

  constant_source = fill(2.0, cellsize(source_grid))
  constant_destination = remap_scalar(cache, constant_source; fill_value=0.0)
  destination_mass_1 = sum(constant_destination .* cellvolumes(destination_grid))
  source_overlap_mass_1 = source_overlap_mass(cache, constant_source)
  source_total_mass_1 = sum(constant_source .* cellvolumes(source_grid))

  @test isapprox(destination_mass_1, source_overlap_mass_1; rtol=1e-10, atol=1e-10)
  @test source_overlap_mass_1 < source_total_mass_1

  ni, nj = cellsize(source_grid)
  variable_source = [sin(0.17 * i) + cos(0.11 * j) for i in 1:ni, j in 1:nj]
  variable_destination_1 = remap_scalar(cache, variable_source)
  variable_destination_2 = similar(variable_destination_1)
  remap_scalar!(variable_destination_2, cache, variable_source)

  @test maximum(abs.(variable_destination_1 .- variable_destination_2)) ≤ 1e-12

  destination_mass_2 = sum(variable_destination_1 .* cellvolumes(destination_grid))
  source_overlap_mass_2 = source_overlap_mass(cache, variable_source)
  @test isapprox(destination_mass_2, source_overlap_mass_2; rtol=1e-10, atol=1e-10)

  @test_throws ArgumentError remap_scalar(cache, ones(ni + 1, nj))
  @test_throws ArgumentError remap_scalar!(
    zeros(size(variable_destination_1, 1) + 1, nj), cache, constant_source
  )

  @test validate_remap_cache(cache, source_grid, destination_grid)
end

@testset "RemappingSchemes wavy 2D to uniform 2D with VTK output" begin
  xwavy(t, xi, eta, p) = xi + p.ax * sin(2pi * eta / p.ly)
  ywavy(t, xi, eta, p) = eta + p.ay * sin(2pi * xi / p.lx)

  source_grid = MappedGrid(
    xwavy,
    ywavy,
    (; ax=0.75, ay=0.5, lx=30.0, ly=20.0),
    (30, 20),
    2;
    compute_metrics=false,
    cache_mode=:off,
  )

  nx, ny = (31, 21)
  xnodes = [Float64(i) for i in 1:nx, j in 1:ny]
  ynodes = [Float64(j) for i in 1:nx, j in 1:ny]
  destination_grid = DiscreteGrid(xnodes, ynodes, 2; compute_metrics=false, cache_mode=:off)

  source_shape = cellsize(source_grid)
  destination_shape = cellsize(destination_grid)
  source_xi = [Float64(i) + 0.5 for i in 1:source_shape[1], j in 1:source_shape[2]]

  cache = build_remap_cache(source_grid, destination_grid; quadrature_order=3)
  destination_scalar = remap_scalar(cache, source_xi; fill_value=0.0)

  @test size(destination_scalar) == destination_shape
  @test maximum(abs.(destination_scalar)) > 0.0

  destination_volumes = cellvolumes(destination_grid)
  mass_destination = sum(destination_scalar .* destination_volumes)
  mass_source_overlap = source_overlap_mass(cache, source_xi)
  @test isapprox(mass_destination, mass_source_overlap; rtol=5e-10, atol=5e-10)

  destination_xi_exact = [
    Float64(i) + 0.5 for i in 1:destination_shape[1], j in 1:destination_shape[2]
  ]
  interior = CartesianIndices((4:(destination_shape[1] - 3), 4:(destination_shape[2] - 3)))
  # interior_abs_error = abs.(destination_scalar[interior] .- destination_xi_exact[interior])
  # @test sum(interior_abs_error) / length(interior_abs_error) < 0.35

  remapped_full = zeros(
    eltype(destination_grid), size(destination_grid.iterators.cell.full)
  )
  remapped_full[destination_grid.iterators.cell.domain] .= destination_scalar
  xi_exact_full = zeros(
    eltype(destination_grid), size(destination_grid.iterators.cell.full)
  )
  xi_exact_full[destination_grid.iterators.cell.domain] .= destination_xi_exact

  out_dir = joinpath(dirname(@__DIR__), "vtk_remap_outputs")
  mkpath(out_dir)
  out_fn_source = joinpath(out_dir, "source")
  out_fn_destination = joinpath(out_dir, "dest")

  CurvilinearGrids.VTKOutput.save_vtk(source_grid, out_fn_source; include_metrics=false)
  CurvilinearGrids.VTKOutput.save_vtk(
    destination_grid,
    out_fn_destination;
    include_metrics=false,
    extra_cell_data=(; xi_remapped=remapped_full, xi_exact=xi_exact_full),
  )

  vtk_file = if isfile(out_fn_destination * ".vts")
    (out_fn_destination * ".vts")
  else
    (out_fn_destination * ".vti")
  end
  @test isfile(vtk_file)
  vtk_text = read(vtk_file, String)
  @test occursin("xi_remapped", vtk_text)
  @test occursin("xi_exact", vtk_text)
end
