using Test
using CurvilinearGrids

@testset "UnifiedGrid traits and adapters" begin
  x = collect(range(0.0, 1.0; length=8))
  legacy = CurvilinearGrid1D(x, :meg6)
  dgrid = DiscreteGrid(legacy; interpolation=:linear, cache_mode=:eager)

  @test dgrid isa DiscreteGrid
  @test dgrid.interpolation === :linear
  @test coordinate_system(dgrid) isa CurvilinearCS
  @test basis_trait(dgrid) isa ContravariantBasis

  @test dgrid.metric_caches.cell.valid
  @test dgrid.metric_caches.face.valid

  invalidate_cell_metrics!(dgrid)
  @test !dgrid.metric_caches.cell.valid
  @test dgrid.metric_caches.face.valid

  refresh_cell_metrics!(dgrid)
  @test dgrid.metric_caches.cell.valid
  @test dgrid.metric_caches.face.valid

  @test_throws ArgumentError DiscreteGrid(legacy; interpolation=:cubic)
end

@testset "MappedGrid and independent cache refresh" begin
  xmap(t, i, p) = p.amp * i + t
  params = (amp=1.0,)

  mgrid = MappedGrid(xmap, params, (8,), :meg6; cache_mode=:eager)
  @test mgrid isa MappedGrid
  @test coordinate_system(mgrid) isa CurvilinearCS
  @test basis_trait(mgrid) isa ContravariantBasis
  @test mgrid.metric_caches.cell.valid
  @test mgrid.metric_caches.face.valid

  update!(mgrid, 1.0, params)
  @test !mgrid.metric_caches.cell.valid
  @test !mgrid.metric_caches.face.valid

  refresh_face_metrics!(mgrid)
  @test !mgrid.metric_caches.cell.valid
  @test mgrid.metric_caches.face.valid
end

@testset "Face metrics include conservative hatted terms" begin
  nx, ny = (6, 7)
  x = zeros(Float64, nx, ny)
  y = zeros(Float64, nx, ny)
  for j in 1:ny
    for i in 1:nx
      x[i, j] = i
      y[i, j] = j
    end
  end

  dgrid2d = DiscreteGrid(x, y, :meg6; cache_mode=:eager)
  fm = face_metrics(dgrid2d)

  edge_name = first(keys(fm))
  edge_metrics = getproperty(fm, edge_name)

  @test Symbol("ξ̂") in keys(edge_metrics)
  @test Symbol("η̂") in keys(edge_metrics)
end

@testset "OrthogonalGrid trait behavior" begin
  r = collect(range(1.0, 2.0; length=6))
  theta = collect(range(0.2, 1.2; length=6))
  phi = collect(range(0.1, 1.0; length=6))
  legacy = SphericalGrid3D(r, theta, phi, 1)
  ogrid = OrthogonalGrid(legacy)

  @test ogrid isa OrthogonalGrid
  @test coordinate_system(ogrid) isa SphericalCS
  @test_throws ArgumentError basis_trait(ogrid)
  @test_throws ArgumentError cell_metrics(ogrid)
  @test_throws ArgumentError face_metrics(ogrid)
end
