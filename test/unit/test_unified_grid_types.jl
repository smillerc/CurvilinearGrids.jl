using Test
using CurvilinearGrids

@testset "UnifiedGrid traits and adapters" begin
  x = collect(range(0.0, 1.0; length=8))
  dgrid = DiscreteGrid(x, 5; interpolation=:linear, cache_mode=:eager)

  @test dgrid isa DiscreteGrid
  @test dgrid isa DiscreteGrid{1,Float64}
  @test isconcretetype(typeof(dgrid))
  @test !ismutabletype(typeof(dgrid))
  @test !hasproperty(dgrid, :core)
  @test !hasproperty(dgrid, :legacy)
  @test !hasproperty(dgrid, :discretization_scheme)
  @test !hasproperty(dgrid, :discretization_scheme_name)
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

  dgrid_spherical = DiscreteGrid(x, 5; basis=SphericalBasis())
  @test basis_trait(dgrid_spherical) isa SphericalBasis
  dgrid_nhalo3 = DiscreteGrid(x, 3)
  @test dgrid_nhalo3.nhalo == 3
  @test_throws MethodError DiscreteGrid(x; interpolation=:linear)

  invalidate_cell_metrics!(dgrid)
  @test !dgrid.metric_caches.cell.valid
  @test dgrid.metric_caches.face.valid

  refresh_cell_metrics!(dgrid)
  @test dgrid.metric_caches.cell.valid
  @test dgrid.metric_caches.face.valid

  @test_throws ArgumentError DiscreteGrid(x, 5; interpolation=:cubic)
end

@testset "Conserved face interpolation scheme selection" begin
  xmap(t, ξ, η, p) = ξ + 0.15 * sin(0.3 * ξ) * cos(0.2 * η)
  ymap(t, ξ, η, p) = η + 0.1 * cos(0.2 * ξ) * sin(0.3 * η)

  g_o3 = MappedGrid(
    xmap,
    ymap,
    (;),
    (12, 12),
    5;
    cache_mode=:eager,
    conserved_metric_scheme=EdgeInterpolationOrder3(),
  )
  g_o2 = MappedGrid(
    xmap,
    ymap,
    (;),
    (12, 12),
    5;
    cache_mode=:eager,
    conserved_metric_scheme=EdgeInterpolationOrder2(),
  )
  g_o1 = MappedGrid(
    xmap,
    ymap,
    (;),
    (12, 12),
    5;
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
    5;
    cache_mode=:eager,
    conserved_metric_scheme=EdgeInterpolationOrder2(),
    interpolation=:linear,
  )

  @test face_metrics(d_o2)[1].conserved[I] isa ConservedMetric{2,Float64}
  @test_throws TypeError MappedGrid(
    xmap, ymap, (;), (8, 8), 5; cache_mode=:eager, conserved_metric_scheme=:not_a_scheme
  )
  @test_throws TypeError MappedGrid(
    xmap, ymap, (;), (8, 8), 5; cache_mode=:eager, conserved_metric_scheme=4
  )
end

@testset "MappedGrid and independent cache refresh" begin
  xmap(t, i, p) = p.amp * i + t
  params = (amp=1.0,)

  mgrid = MappedGrid(xmap, params, (8,), 5; cache_mode=:eager)
  @test mgrid isa MappedGrid
  @test mgrid isa MappedGrid{1,Float64}
  @test isconcretetype(typeof(mgrid))
  @test !ismutabletype(typeof(mgrid))
  @test !hasproperty(mgrid, :core)
  @test !hasproperty(mgrid, :legacy)
  @test !hasproperty(mgrid, :discretization_scheme)
  @test !hasproperty(mgrid, :discretization_scheme_name)
  @test coordinate_system(mgrid) isa CurvilinearCS
  @test basis_trait(mgrid) isa CartesianBasis
  @test fieldtype(typeof(mgrid), :mapping_functions) !== Any
  @test fieldtype(typeof(mgrid), :metric_functions_cache) !== Any
  @test fieldtype(typeof(mgrid), :metric_caches) !== Any
  @test fieldtype(typeof(mgrid.metric_caches.cell), :data) !== Any
  @test fieldtype(typeof(mgrid.metric_caches.face), :data) !== Any
  @test mgrid.metric_caches.cell.valid
  @test mgrid.metric_caches.face.valid

  mgrid_spherical = MappedGrid(xmap, params, (8,), 5; basis=SphericalBasis())
  @test basis_trait(mgrid_spherical) isa SphericalBasis
  mgrid_nhalo3 = MappedGrid(xmap, params, (8,), 3)
  @test mgrid_nhalo3.nhalo == 3
  @test_throws MethodError MappedGrid(xmap, params, (8,); cache_mode=:eager)

  update!(mgrid, 1.0, params)
  @test !mgrid.metric_caches.cell.valid
  @test !mgrid.metric_caches.face.valid

  refresh_face_metrics!(mgrid)
  @test !mgrid.metric_caches.cell.valid
  @test mgrid.metric_caches.face.valid
end

@testset "Metrics can be fully disabled" begin
  xmap(t, ξ, p) = ξ + p.shift
  mgrid = MappedGrid(xmap, (; shift=1.0), (8,), 5; compute_metrics=false, cache_mode=:off)

  @test mgrid.metric_functions_cache !== nothing
  @test mgrid.metric_caches === nothing
  c0 = coord(mgrid, (1,))
  @test isfinite(c0[1])
  update!(mgrid, 1.0, (; shift=2.0))
  c1 = coord(mgrid, (1,))
  @test c1[1] ≈ c0[1] + 1.0

  @test_throws ArgumentError cell_metrics(mgrid)
  @test_throws ArgumentError face_metrics(mgrid)
  @test forward_cell_metrics(mgrid, (1.2,)) isa Metric{1,Float64}
  @test inverse_cell_metrics(mgrid, (1.2,)) isa Metric{1,Float64}
  @test isfinite(cellvolume(mgrid, (1.2,)))
  @test forward_cell_metrics(mgrid, (1,)) isa Metric{1,Float64}
  @test inverse_cell_metrics(mgrid, (1,)) isa Metric{1,Float64}
  @test isfinite(cellvolume(mgrid, (1,)))

  x = collect(range(0.0, 1.0; length=8))
  dgrid = DiscreteGrid(x, 5; compute_metrics=false, cache_mode=:off)
  @test dgrid.metric_functions_cache !== nothing
  @test dgrid.metric_caches === nothing
  @test isfinite(coord(dgrid, (1,))[1])

  @test_throws ArgumentError cell_metrics(dgrid)
  @test_throws ArgumentError face_metrics(dgrid)
  @test forward_cell_metrics(dgrid, (1.2,)) isa Metric{1,Float64}
  @test inverse_cell_metrics(dgrid, (1.2,)) isa Metric{1,Float64}
  @test isfinite(cellvolume(dgrid, (1.2,)))
  @test forward_cell_metrics(dgrid, (1,)) isa Metric{1,Float64}
  @test inverse_cell_metrics(dgrid, (1,)) isa Metric{1,Float64}
  @test isfinite(cellvolume(dgrid, (1,)))
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

  dgrid2d = DiscreteGrid(x, y, 5; cache_mode=:eager)
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
    x1d, 5; coordinate_system=CylindricalCS(), basis=CartesianBasis(), cache_mode=:eager
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
    x2, y2, 5; coordinate_system=CylindricalCS(), basis=CartesianBasis(), cache_mode=:eager
  )
  I2 = first(cyl2d.iterators.cell.domain)
  J2 = cell_metrics(cyl2d).forward[I2].J
  r2 = centroid(cyl2d, I2)[1]
  @test cellvolume(cyl2d, I2) ≈ (2 * π * r2) * J2

  axi_y = DiscreteGrid(
    x2,
    y2,
    5;
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
    5;
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
    5;
    coordinate_system=SphericalCS(),
    basis=SphericalBasis(),
    cache_mode=:eager,
  )
  I3 = first(sph3d.iterators.cell.domain)
  J3 = cell_metrics(sph3d).forward[I3].J
  c3 = centroid(sph3d, I3)
  @test cellvolume(sph3d, I3) ≈ (c3[1]^2 * sin(c3[2])) * J3

  sph2d = DiscreteGrid(
    x2, y2, 5; coordinate_system=SphericalCS(), basis=SphericalBasis(), cache_mode=:eager
  )
  Is2 = first(sph2d.iterators.cell.domain)
  Js2 = cell_metrics(sph2d).forward[Is2].J
  cs2 = centroid(sph2d, Is2)
  @test cellvolume(sph2d, Is2) ≈ (cs2[1]^2 * sin(cs2[2])) * Js2
end

@testset "Mixed Real/Int tuple indexing" begin
  idx = (1, 2.3, 4)

  xmap(t, ξ, η, ζ, p) = ξ
  ymap(t, ξ, η, ζ, p) = 2 * η
  zmap(t, ξ, η, ζ, p) = 3 * ζ

  mgrid = MappedGrid(
    xmap, ymap, zmap, (;), (7, 7, 7), 5; basis=CartesianBasis(), cache_mode=:eager
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
  @test inv(Fm_discrete).J ≈ inv(Gm_discrete.J)
  @test inv(Fm_continuous).jacobian_matrix ≈ Gm_continuous.jacobian_matrix
  @test inv(Fm_continuous).J ≈ inv(Gm_continuous.J)

  nx, ny, nz = (8, 8, 8)
  x = [Float64(i) for i in 1:nx, j in 1:ny, k in 1:nz]
  y = [2.0 * j for i in 1:nx, j in 1:ny, k in 1:nz]
  z = [3.0 * k for i in 1:nx, j in 1:ny, k in 1:nz]

  dgrid = DiscreteGrid(
    x, y, z, 5; basis=CartesianBasis(), interpolation=:linear, cache_mode=:eager
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
  @test inv(Fd_discrete).J ≈ inv(Gd_discrete.J)
  @test inv(Fd_continuous).jacobian_matrix ≈ Gd_continuous.jacobian_matrix
  @test inv(Fd_continuous).J ≈ inv(Gd_continuous.J)
end

@testset "Inverse mapping: physical -> computational coordinate" begin
  params = (a11=1.3, a12=0.2, a21=-0.4, a22=1.1, b1=0.5, b2=-0.2)
  x1(t, ξ, η, p) = p.a11 * ξ + p.a12 * η + p.b1
  x2(t, ξ, η, p) = p.a21 * ξ + p.a22 * η + p.b2

  mgrid = MappedGrid(x1, x2, params, (9, 10), 5; cache_mode=:off, compute_metrics=false)

  ξtrue = (3.2, 5.7)
  xphys = Tuple(coord(mgrid, ξtrue))
  ξmapped = computational_coordinate(mgrid, xphys)
  @test all(abs.(Tuple(ξmapped) .- ξtrue) .<= 1e-10)

  result = computational_coordinate(
    mgrid, xphys; guess=(1.0, 1.0), return_result=true, throw_on_failure=false
  )
  @test result isa InverseCoordinateResult{2,Float64}
  @test result.converged
  @test all(abs.(Tuple(result.coordinate) .- ξtrue) .<= 1e-10)

  x = [x1(0.0, i, j, params) for i in 1:10, j in 1:11]
  y = [x2(0.0, i, j, params) for i in 1:10, j in 1:11]
  dgrid = DiscreteGrid(
    x, y, 5; interpolation=:linear, cache_mode=:off, compute_metrics=false
  )
  ξdiscrete = computational_coordinate(dgrid, xphys)
  @test all(abs.(Tuple(ξdiscrete) .- ξtrue) .<= 1e-10)

  @test_throws ArgumentError computational_coordinate(mgrid, (-100.0, -100.0))

  x1n(t, ξ, η, p) = ξ + 0.35 * sin(0.8 * ξ) * cos(0.4 * η)
  x2n(t, ξ, η, p) = η + 0.3 * cos(0.6 * ξ) * sin(0.7 * η)
  nonlinear_grid = MappedGrid(
    x1n, x2n, (;), (9, 10), 5; cache_mode=:off, compute_metrics=false
  )
  ξnonlinear = (7.9, 8.2)
  xnonlinear = Tuple(coord(nonlinear_grid, ξnonlinear))
  @test_throws ErrorException computational_coordinate(
    nonlinear_grid, xnonlinear; guess=(1.0, 1.0), maxiters=1, tol=1e-14
  )
  failure = computational_coordinate(
    nonlinear_grid,
    xnonlinear;
    guess=(1.0, 1.0),
    maxiters=1,
    tol=1e-14,
    throw_on_failure=false,
    return_result=true,
  )
  @test !failure.converged
end

@testset "VTK output for unified grids" begin
  xmap(t, ξ, η, p) = ξ + 0.1 * sin(0.25 * η)
  ymap(t, ξ, η, p) = η + 0.05 * cos(0.2 * ξ)

  mapped2d = MappedGrid(
    xmap,
    ymap,
    (;),
    (8, 7),
    2;
    cache_mode=:eager,
    conserved_metric_scheme=EdgeInterpolationOrder2(),
  )

  nx, ny = (9, 8)
  x2 = [Float64(i) for i in 1:nx, j in 1:ny]
  y2 = [Float64(j) for i in 1:nx, j in 1:ny]
  discrete2d = DiscreteGrid(x2, y2, 2; compute_metrics=false, cache_mode=:off)

  nx3, ny3, nz3 = (5, 4, 3)
  x3 = [Float64(i) for i in 1:nx3, j in 1:ny3, k in 1:nz3]
  y3 = [Float64(j) for i in 1:nx3, j in 1:ny3, k in 1:nz3]
  z3 = [Float64(k) for i in 1:nx3, j in 1:ny3, k in 1:nz3]
  discrete3d = DiscreteGrid(x3, y3, z3, 1; compute_metrics=false, cache_mode=:off)

  mktempdir() do dir
    mapped_fn = joinpath(dir, "mapped2d")
    save_vtk(mapped2d, mapped_fn)
    @test any(startswith("mapped2d"), readdir(dir))

    discrete_fn = joinpath(dir, "discrete2d")
    payload2d = zeros(eltype(discrete2d), size(discrete2d.iterators.cell.full))
    payload2d[discrete2d.iterators.cell.domain] .= 1.0
    save_vtk(
      discrete2d, discrete_fn; include_metrics=false, extra_cell_data=(; marker=payload2d)
    )
    @test any(startswith("discrete2d"), readdir(dir))
    @test_throws ArgumentError save_vtk(discrete2d, joinpath(dir, "discrete2d_metrics"))

    discrete3d_fn = joinpath(dir, "discrete3d")
    save_vtk(discrete3d, discrete3d_fn; include_metrics=false)
    @test any(startswith("discrete3d"), readdir(dir))
  end
end
