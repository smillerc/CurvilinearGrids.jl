using Test
using CurvilinearGrids

@testset "UnifiedGrid traits and adapters" begin
  x = collect(range(0.0, 1.0; length=8))
  dgrid = DiscreteGrid(x, :meg6; interpolation=:linear, cache_mode=:eager)

  @test dgrid isa DiscreteGrid
  @test dgrid isa DiscreteGrid{1,Float64}
  @test isconcretetype(typeof(dgrid))
  @test !ismutabletype(typeof(dgrid))
  @test !hasproperty(dgrid, :core)
  @test !hasproperty(dgrid, :legacy)
  @test dgrid.interpolation === :linear
  @test coordinate_system(dgrid) isa CurvilinearCS
  @test basis_trait(dgrid) isa CartesianBasis
  @test fieldtype(typeof(dgrid), :mapping_functions) !== Any
  @test fieldtype(typeof(dgrid), :metric_functions_cache) !== Any
  @test fieldtype(typeof(dgrid), :metric_caches) !== Any
  @test fieldtype(typeof(dgrid.metric_caches.cell), :data) !== Any
  @test fieldtype(typeof(dgrid.metric_caches.face), :data) !== Any

  @test dgrid.metric_caches.cell.valid
  @test dgrid.metric_caches.face.valid

  dgrid_spherical = DiscreteGrid(x, :meg6; basis=SphericalBasis())
  @test basis_trait(dgrid_spherical) isa SphericalBasis

  invalidate_cell_metrics!(dgrid)
  @test !dgrid.metric_caches.cell.valid
  @test dgrid.metric_caches.face.valid

  refresh_cell_metrics!(dgrid)
  @test dgrid.metric_caches.cell.valid
  @test dgrid.metric_caches.face.valid

  @test_throws ArgumentError DiscreteGrid(x, :meg6; interpolation=:cubic)
end

@testset "Conserved face interpolation scheme selection" begin
  xmap(t, ξ, η, p) = ξ + 0.15 * sin(0.3 * ξ) * cos(0.2 * η)
  ymap(t, ξ, η, p) = η + 0.1 * cos(0.2 * ξ) * sin(0.3 * η)

  g_o3 = MappedGrid(
    xmap,
    ymap,
    (;),
    (12, 12),
    :meg6;
    cache_mode=:eager,
    conserved_metric_scheme=EdgeInterpolationOrder3(),
  )
  g_o2 = MappedGrid(
    xmap,
    ymap,
    (;),
    (12, 12),
    :meg6;
    cache_mode=:eager,
    conserved_metric_scheme=EdgeInterpolationOrder2(),
  )
  g_o1 = MappedGrid(
    xmap,
    ymap,
    (;),
    (12, 12),
    :meg6;
    cache_mode=:eager,
    conserved_metric_scheme=EdgeInterpolationOrder1(),
  )

  I = first(g_o3.iterators.cell.domain)
  c_o3 = face_metrics(g_o3)[1].conserved[I].jacobian_matrix
  c_o2 = face_metrics(g_o2)[1].conserved[I].jacobian_matrix
  c_o1 = face_metrics(g_o1)[1].conserved[I].jacobian_matrix
  @test c_o3 != c_o2
  @test c_o3 != c_o1

  nx, ny = (13, 13)
  x = [xmap(0.0, i, j, (;)) for i in 1:nx, j in 1:ny]
  y = [ymap(0.0, i, j, (;)) for i in 1:nx, j in 1:ny]
  d_o2 = DiscreteGrid(
    x,
    y,
    :meg6;
    cache_mode=:eager,
    conserved_metric_scheme=EdgeInterpolationOrder2(),
    interpolation=:linear,
  )

  @test face_metrics(d_o2)[1].conserved[I] isa ConservedMetric{2,Float64}
  @test_throws TypeError MappedGrid(
    xmap, ymap, (;), (8, 8), :meg6; cache_mode=:eager, conserved_metric_scheme=:not_a_scheme
  )
  @test_throws TypeError MappedGrid(
    xmap, ymap, (;), (8, 8), :meg6; cache_mode=:eager, conserved_metric_scheme=4
  )
end

@testset "MappedGrid and independent cache refresh" begin
  xmap(t, i, p) = p.amp * i + t
  params = (amp=1.0,)

  mgrid = MappedGrid(xmap, params, (8,), :meg6; cache_mode=:eager)
  @test mgrid isa MappedGrid
  @test mgrid isa MappedGrid{1,Float64}
  @test isconcretetype(typeof(mgrid))
  @test !ismutabletype(typeof(mgrid))
  @test !hasproperty(mgrid, :core)
  @test !hasproperty(mgrid, :legacy)
  @test coordinate_system(mgrid) isa CurvilinearCS
  @test basis_trait(mgrid) isa CartesianBasis
  @test fieldtype(typeof(mgrid), :mapping_functions) !== Any
  @test fieldtype(typeof(mgrid), :metric_functions_cache) !== Any
  @test fieldtype(typeof(mgrid), :metric_caches) !== Any
  @test fieldtype(typeof(mgrid.metric_caches.cell), :data) !== Any
  @test fieldtype(typeof(mgrid.metric_caches.face), :data) !== Any
  @test mgrid.metric_caches.cell.valid
  @test mgrid.metric_caches.face.valid

  mgrid_spherical = MappedGrid(xmap, params, (8,), :meg6; basis=SphericalBasis())
  @test basis_trait(mgrid_spherical) isa SphericalBasis

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

  I = first(dgrid2d.iterators.cell.domain)
  # face_metrics[edge_dim].{forward,inverse,conserved}[grid_index]
  i_face_metric = fm[1].conserved[I]

  @test i_face_metric isa ConservedMetric{2,Float64}
  @test isfinite(i_face_metric.jacobian_matrix[1, 1])
  @test isfinite(i_face_metric.jacobian_matrix[2, 2])
end

@testset "OrthogonalGrid trait behavior" begin
  r = collect(range(1.0, 2.0; length=6))
  theta = collect(range(0.2, 1.2; length=6))
  phi = collect(range(0.1, 1.0; length=6))
  legacy = SphericalGrid3D(r, theta, phi, 1)
  ogrid = OrthogonalGrid(legacy)

  @test ogrid isa OrthogonalGrid
  @test ogrid isa OrthogonalGrid{3,Float64}
  @test !ismutabletype(typeof(ogrid))
  @test coordinate_system(ogrid) isa SphericalCS
  @test_throws ArgumentError basis_trait(ogrid)
  @test_throws ArgumentError cell_metrics(ogrid)
  @test_throws ArgumentError face_metrics(ogrid)
end

@testset "Unified cellvolume trait dispatch" begin
  x1d = collect(range(1.0, 3.0; length=8))
  cyl1d = DiscreteGrid(
    x1d, :meg6; coordinate_system=CylindricalCS(), basis=CartesianBasis(), cache_mode=:eager
  )
  I1 = first(cyl1d.iterators.cell.domain)
  J1 = cell_metrics(cyl1d).forward[I1].J
  r1 = centroid(cyl1d, I1)[1]
  @test cellvolume(cyl1d, I1) ≈ (2 * π * r1) * J1

  nx, ny = (6, 7)
  x2 = zeros(Float64, nx, ny)
  y2 = zeros(Float64, nx, ny)
  for j in 1:ny
    for i in 1:nx
      x2[i, j] = i
      y2[i, j] = j
    end
  end

  cyl2d = DiscreteGrid(
    x2,
    y2,
    :MEG6;
    coordinate_system=CylindricalCS(),
    basis=CartesianBasis(),
    cache_mode=:eager,
  )
  I2 = first(cyl2d.iterators.cell.domain)
  J2 = cell_metrics(cyl2d).forward[I2].J
  r2 = centroid(cyl2d, I2)[1]
  @test cellvolume(cyl2d, I2) ≈ (2 * π * r2) * J2

  axi_y = DiscreteGrid(
    x2,
    y2,
    :MEG6;
    coordinate_system=AxisymmetricCS{:y}(),
    basis=CartesianBasis(),
    cache_mode=:eager,
  )
  Iay = first(axi_y.iterators.cell.domain)
  Jay = cell_metrics(axi_y).forward[Iay].J
  ray = centroid(axi_y, Iay)[1]
  @test cellvolume(axi_y, Iay) ≈ (2 * π * ray) * Jay

  axi_x = DiscreteGrid(
    x2,
    y2,
    :MEG6;
    coordinate_system=AxisymmetricCS{:x}(),
    basis=CartesianBasis(),
    cache_mode=:eager,
  )
  Iax = first(axi_x.iterators.cell.domain)
  Jax = cell_metrics(axi_x).forward[Iax].J
  rax = centroid(axi_x, Iax)[2]
  @test cellvolume(axi_x, Iax) ≈ (2 * π * rax) * Jax

  params = (r0=1.0, θ0=0.2, ϕ0=0.1, Δr=0.1, Δθ=0.05, Δϕ=0.07)
  rmap(t, ξ, η, ζ, p) = p.r0 + (ξ - 1) * p.Δr
  θmap(t, ξ, η, ζ, p) = p.θ0 + (η - 1) * p.Δθ
  ϕmap(t, ξ, η, ζ, p) = p.ϕ0 + (ζ - 1) * p.Δϕ
  sph3d = MappedGrid(
    rmap,
    θmap,
    ϕmap,
    params,
    (6, 6, 6),
    :meg6;
    coordinate_system=SphericalCS(),
    basis=SphericalBasis(),
    cache_mode=:eager,
  )
  I3 = first(sph3d.iterators.cell.domain)
  J3 = cell_metrics(sph3d).forward[I3].J
  c3 = centroid(sph3d, I3)
  @test cellvolume(sph3d, I3) ≈ (c3[1]^2 * sin(c3[2])) * J3

  sph2d = DiscreteGrid(
    x2,
    y2,
    :MEG6;
    coordinate_system=SphericalCS(),
    basis=SphericalBasis(),
    cache_mode=:eager,
  )
  Is2 = first(sph2d.iterators.cell.domain)
  @test_throws ArgumentError cellvolume(sph2d, Is2)
end

@testset "Mixed Real/Int tuple indexing" begin
  idx = (1, 2.3, 4)

  xmap(t, ξ, η, ζ, p) = ξ
  ymap(t, ξ, η, ζ, p) = 2 * η
  zmap(t, ξ, η, ζ, p) = 3 * ζ

  mgrid = MappedGrid(
    xmap, ymap, zmap, (;), (7, 7, 7), :meg6; basis=CartesianBasis(), cache_mode=:eager
  )

  cm = coord(mgrid, idx)
  Jm = jacobian_matrix(mgrid, idx)
  Fm_discrete = forward_cell_metrics(mgrid, (1, 2, 4))
  Fm_continuous = forward_cell_metrics(mgrid, idx)
  Gm_discrete = inverse_cell_metrics(mgrid, (1, 2, 4))
  Gm_continuous = inverse_cell_metrics(mgrid, idx)
  Vm = cellvolume(mgrid, idx)

  @test cm ≈ [1.0, 4.6, 12.0]
  @test_throws MethodError centroid(mgrid, idx)
  @test Jm[1, 1] ≈ 1.0
  @test Jm[2, 2] ≈ 2.0
  @test Jm[3, 3] ≈ 3.0
  @test Jm[1, 2] ≈ 0.0
  @test Vm ≈ 6.0
  @test Fm_discrete isa Metric{3,Float64}
  @test Fm_continuous isa Metric{3,Float64}
  @test Gm_discrete isa Metric{3,Float64}
  @test Gm_continuous isa Metric{3,Float64}
  @test inv(Fm_discrete).jacobian_matrix ≈ Gm_discrete.jacobian_matrix
  @test inv(Fm_discrete).J ≈ Gm_discrete.J
  @test inv(Fm_continuous).jacobian_matrix ≈ Gm_continuous.jacobian_matrix
  @test inv(Fm_continuous).J ≈ Gm_continuous.J

  nx, ny, nz = (8, 8, 8)
  x = [Float64(i) for i in 1:nx, j in 1:ny, k in 1:nz]
  y = [2.0 * j for i in 1:nx, j in 1:ny, k in 1:nz]
  z = [3.0 * k for i in 1:nx, j in 1:ny, k in 1:nz]

  dgrid = DiscreteGrid(
    x, y, z, :meg6; basis=CartesianBasis(), interpolation=:linear, cache_mode=:eager
  )

  cd = coord(dgrid, idx)
  Jd = jacobian_matrix(dgrid, idx)
  Fd_discrete = forward_cell_metrics(dgrid, (1, 2, 4))
  Fd_continuous = forward_cell_metrics(dgrid, idx)
  Gd_discrete = inverse_cell_metrics(dgrid, (1, 2, 4))
  Gd_continuous = inverse_cell_metrics(dgrid, idx)
  Vd = cellvolume(dgrid, idx)

  @test cd ≈ [1.0, 4.6, 12.0]
  @test_throws MethodError centroid(dgrid, idx)
  @test Jd[1, 1] ≈ 1.0
  @test Jd[2, 2] ≈ 2.0
  @test Jd[3, 3] ≈ 3.0
  @test Jd[2, 1] ≈ 0.0
  @test Vd ≈ 6.0
  @test Fd_discrete isa Metric{3,Float64}
  @test Fd_continuous isa Metric{3,Float64}
  @test Gd_discrete isa Metric{3,Float64}
  @test Gd_continuous isa Metric{3,Float64}
  @test inv(Fd_discrete).jacobian_matrix ≈ Gd_discrete.jacobian_matrix
  @test inv(Fd_discrete).J ≈ Gd_discrete.J
  @test inv(Fd_continuous).jacobian_matrix ≈ Gd_continuous.jacobian_matrix
  @test inv(Fd_continuous).J ≈ Gd_continuous.J
end
