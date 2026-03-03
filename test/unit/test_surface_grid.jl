using LinearAlgebra
using Test
using CurvilinearGrids

@testset "SurfaceGrid boundary parsing and validation" begin
  x2(t, ξ, η, p) = ξ
  y2(t, ξ, η, p) = η
  mapped2 = MappedGrid(x2, y2, (;), (8, 9), 2; compute_metrics=false, cache_mode=:off)

  s1 = SurfaceGrid(mapped2, :ILO)
  @test s1.axis == 1
  @test s1.side == :lo

  s2 = SurfaceGrid(mapped2, :xmax)
  @test s2.axis == 1
  @test s2.side == :hi

  s3 = SurfaceGrid(mapped2, :YLo)
  @test s3.axis == 2
  @test s3.side == :lo

  @test_throws ArgumentError SurfaceGrid(mapped2, :bad_boundary)
  @test_throws ArgumentError SurfaceGrid(mapped2, :klo)
  @test_throws ArgumentError CurvilinearGrids._build_surface_index((1,), 1, 1, Val(3))

  x1(t, ξ, p) = ξ
  mapped1 = MappedGrid(x1, (;), (8,), 2; compute_metrics=false, cache_mode=:off)
  @test_throws ArgumentError SurfaceGrid(mapped1, :ilo)
end

@testset "SurfaceGrid mapped 2D geometry and orientation" begin
  x2(t, ξ, η, p) = ξ
  y2(t, ξ, η, p) = η
  mapped2 = MappedGrid(x2, y2, (;), (8, 9), 2; compute_metrics=false, cache_mode=:off)

  out = SurfaceGrid(mapped2, :ilo, true)
  inn = SurfaceGrid(mapped2, :ilo, false)

  @test out isa SurfaceGrid{2,Float64}
  @test !hasproperty(out, :parent)
  @test !hasproperty(out, :face_parent_indices)
  @test out.boundary == :ilo
  @test length(out.node_coordinates[1]) == length(mapped2.iterators.node.domain.indices[2])
  @test length(out.face_areas) == length(mapped2.iterators.cell.domain.indices[2])

  I = first(CartesianIndices(out.face_areas))
  @test out.face_areas[I] ≈ 1.0 atol = 1e-12
  @test out.face_normals[1][I] ≈ -1.0 atol = 1e-12
  @test out.face_normals[2][I] ≈ 0.0 atol = 1e-12
  @test hypot(out.face_normals[1][I], out.face_normals[2][I]) ≈ 1.0 atol = 1e-12

  @test inn.face_areas[I] ≈ out.face_areas[I] atol = 1e-12
  @test inn.face_normals[1][I] ≈ -out.face_normals[1][I] atol = 1e-12
  @test inn.face_normals[2][I] ≈ -out.face_normals[2][I] atol = 1e-12

  jhi = SurfaceGrid(mapped2, :jhi, true)
  Ij = first(CartesianIndices(jhi.face_areas))
  @test jhi.face_normals[1][Ij] ≈ 0.0 atol = 1e-12
  @test jhi.face_normals[2][Ij] ≈ 1.0 atol = 1e-12

  out2 = extract_surface_mesh(mapped2, :ilo)
  @test out2 isa SurfaceGrid{2,Float64}
end

@testset "SurfaceGrid mapped 3D geometry and zero-area fallback" begin
  x3(t, ξ, η, ζ, p) = ξ
  y3(t, ξ, η, ζ, p) = η
  z3(t, ξ, η, ζ, p) = ζ
  mapped3 = MappedGrid(
    x3, y3, z3, (;), (6, 7, 8), 2; compute_metrics=false, cache_mode=:off
  )

  khi_out = SurfaceGrid(mapped3, :khi, true)
  khi_in = SurfaceGrid(mapped3, :khi, false)
  klo_out = SurfaceGrid(mapped3, :klo, true)
  zhi_out = SurfaceGrid(mapped3, :ZHI, true)
  @test zhi_out.side == :hi
  @test zhi_out.axis == 3

  I = first(CartesianIndices(khi_out.face_areas))
  @test khi_out.face_areas[I] ≈ 1.0 atol = 1e-12
  @test khi_out.face_normals[1][I] ≈ 0.0 atol = 1e-12
  @test khi_out.face_normals[2][I] ≈ 0.0 atol = 1e-12
  @test khi_out.face_normals[3][I] ≈ 1.0 atol = 1e-12
  @test khi_in.face_normals[3][I] ≈ -1.0 atol = 1e-12
  @test klo_out.face_normals[3][I] ≈ -1.0 atol = 1e-12

  x0(t, ξ, η, ζ, p) = zero(ξ)
  y0(t, ξ, η, ζ, p) = zero(η)
  z0(t, ξ, η, ζ, p) = zero(ζ)
  degenerate3 = MappedGrid(
    x0, y0, z0, (;), (4, 4, 4), 1; compute_metrics=false, cache_mode=:off
  )
  sdeg3 = SurfaceGrid(degenerate3, :ilo, true)
  Id = first(CartesianIndices(sdeg3.face_areas))
  @test sdeg3.face_areas[Id] == 0.0
  @test sdeg3.face_normals[1][Id] == 0.0
  @test sdeg3.face_normals[2][Id] == 0.0
  @test sdeg3.face_normals[3][Id] == 0.0
end

@testset "Unified face area and outward normal API" begin
  x3(t, ξ, η, ζ, p) = ξ
  y3(t, ξ, η, ζ, p) = η
  z3(t, ξ, η, ζ, p) = ζ
  mapped3 = MappedGrid(
    x3, y3, z3, (;), (6, 7, 8), 2; compute_metrics=false, cache_mode=:off
  )
  cell_ranges = Tuple(mapped3.iterators.cell.domain.indices)
  i0 = first(cell_ranges[1])
  j0 = first(cell_ranges[2])
  k0 = first(cell_ranges[3])
  k1 = last(cell_ranges[3])
  idx_lo = (i0, j0, k0)
  idx_hi = (i0, j0, k1)

  area_lo = face_area(mapped3, idx_lo, :klo)
  area_hi = face_area(mapped3, idx_hi, :khi)
  normal_lo = outward_face_normal(mapped3, idx_lo, :klo)
  normal_hi = outward_face_normal(mapped3, idx_hi, :khi)

  @test area_lo ≈ 1.0 atol = 1e-12
  @test area_hi ≈ 1.0 atol = 1e-12
  @test normal_lo[1] ≈ 0.0 atol = 1e-12
  @test normal_lo[2] ≈ 0.0 atol = 1e-12
  @test normal_lo[3] ≈ -1.0 atol = 1e-12
  @test normal_hi[1] ≈ 0.0 atol = 1e-12
  @test normal_hi[2] ≈ 0.0 atol = 1e-12
  @test normal_hi[3] ≈ 1.0 atol = 1e-12
end

@testset "SurfaceGrid axisymmetric rotated area" begin
  params = (; r0=2.0, dr=0.2, z0=-1.0, dz=0.4)
  rmap(t, ξ, η, p) = p.r0 + p.dr * ξ
  zmap(t, ξ, η, p) = p.z0 + p.dz * η

  axisym = MappedGrid(
    rmap,
    zmap,
    params,
    (6, 7),
    2;
    coordinate_system=AxisymmetricCS{:y}(),
    basis=CartesianBasis(),
    cache_mode=:eager,
  )

  map_surface = SurfaceGrid(axisym, :ihi, true)

  @test all(
    I -> begin
      r_face = map_surface.face_centers[1][I]
      expected_area = 2π * abs(r_face) * abs(params.dz)
      isapprox(map_surface.face_areas[I], expected_area; atol=1e-12, rtol=1e-12)
    end,
    CartesianIndices(map_surface.face_areas),
  )

  I = first(CartesianIndices(map_surface.face_areas))
  @test map_surface.face_normals[1][I] ≈ 1.0 atol = 1e-12
  @test map_surface.face_normals[2][I] ≈ 0.0 atol = 1e-12

  cell_ranges = Tuple(axisym.iterators.cell.domain.indices)
  idx_api = (last(cell_ranges[1]), first(cell_ranges[2]))
  area_api = face_area(axisym, idx_api, :ihi)
  normal_api = outward_face_normal(axisym, idx_api, :ihi)
  @test area_api ≈ map_surface.face_areas[I] atol = 1e-12
  @test normal_api[1] ≈ map_surface.face_normals[1][I] atol = 1e-12
  @test normal_api[2] ≈ map_surface.face_normals[2][I] atol = 1e-12
end

@testset "SurfaceGrid spherical-basis sector geometry" begin
  params = (; r0=1.0, dr=0.05, theta0=0.35, dtheta=0.03, phi0=0.20, dphi=0.04)
  r(t, ξ, η, ζ, p) = p.r0 + p.dr * ξ
  theta(t, ξ, η, ζ, p) = p.theta0 + p.dtheta * η
  phi(t, ξ, η, ζ, p) = p.phi0 + p.dphi * ζ

  spherical = MappedGrid(
    r,
    theta,
    phi,
    params,
    (6, 7, 8),
    2;
    coordinate_system=SphericalCS(),
    basis=SphericalBasis(),
    cache_mode=:eager,
  )

  map_surface = SurfaceGrid(spherical, :ihi, true)

  @test all(
    I -> begin
      cx = map_surface.face_centers[1][I]
      cy = map_surface.face_centers[2][I]
      cz = map_surface.face_centers[3][I]
      rmag = sqrt(cx^2 + cy^2 + cz^2)
      nx_expected = cx / rmag
      ny_expected = cy / rmag
      nz_expected = cz / rmag
      isapprox(map_surface.face_normals[1][I], nx_expected; atol=1e-12) &&
        isapprox(map_surface.face_normals[2][I], ny_expected; atol=1e-12) &&
        isapprox(map_surface.face_normals[3][I], nz_expected; atol=1e-12)
    end,
    CartesianIndices(map_surface.face_areas),
  )

  @test all(
    I -> begin
      cx = map_surface.face_centers[1][I]
      cy = map_surface.face_centers[2][I]
      cz = map_surface.face_centers[3][I]
      r_face = sqrt(cx^2 + cy^2 + cz^2)
      theta_face = acos(clamp(cz / r_face, -1.0, 1.0))
      expected_area = r_face^2 * sin(theta_face) * params.dtheta * params.dphi
      isapprox(map_surface.face_areas[I], expected_area; rtol=1e-12, atol=1e-12)
    end,
    CartesianIndices(map_surface.face_areas),
  )

  cell_ranges = Tuple(spherical.iterators.cell.domain.indices)
  idx_api = (last(cell_ranges[1]), first(cell_ranges[2]), first(cell_ranges[3]))
  Iapi = first(CartesianIndices(map_surface.face_areas))
  area_api = face_area(spherical, idx_api, :ihi)
  normal_api = outward_face_normal(spherical, idx_api, :ihi)
  @test area_api ≈ map_surface.face_areas[Iapi] atol = 1e-12
  @test normal_api[1] ≈ map_surface.face_normals[1][Iapi] atol = 1e-12
  @test normal_api[2] ≈ map_surface.face_normals[2][Iapi] atol = 1e-12
  @test normal_api[3] ≈ map_surface.face_normals[3][Iapi] atol = 1e-12

  cart_pts = CurvilinearGrids._surface_vtk_points(map_surface)
  Inode = first(CartesianIndices(cart_pts))
  node_ranges = Tuple(spherical.iterators.node.domain.indices)
  i = last(node_ranges[1])
  j = first(node_ranges[2])
  k = first(node_ranges[3])
  rnode = spherical.node_coordinates[1][i, j, k]
  thetanode = spherical.node_coordinates[2][i, j, k]
  phinode = spherical.node_coordinates[3][i, j, k]
  expected_x = rnode * sin(thetanode) * cos(phinode)
  expected_y = rnode * sin(thetanode) * sin(phinode)
  expected_z = rnode * cos(thetanode)
  @test isapprox(cart_pts[Inode][1], expected_x; atol=1e-12)
  @test isapprox(cart_pts[Inode][2], expected_y; atol=1e-12)
  @test isapprox(cart_pts[Inode][3], expected_z; atol=1e-12)

  dir = joinpath(dirname(@__DIR__), "vtk_surface_outputs")
  mkpath(dir)
  fn_sphere = joinpath(dir, "surface3d_spherical_sector")
  for ext in (".vti", ".vts")
    isfile(fn_sphere * ext) && rm(fn_sphere * ext)
  end
  save_vtk(map_surface, fn_sphere)
  vtk_sphere = isfile(fn_sphere * ".vts") ? (fn_sphere * ".vts") : (fn_sphere * ".vti")
  @test isfile(vtk_sphere)
  txt_sphere = read(vtk_sphere, String)
  @test occursin("Points", txt_sphere)
  @test occursin("face_area", txt_sphere)
  @test occursin("face_normal", txt_sphere)
end

@testset "SurfaceGrid discrete geometry" begin
  nx, ny = (10, 11)
  x2 = [1.5 * i + 0.1 * j for i in 1:nx, j in 1:ny]
  y2 = [0.2 * i + 2.0 * j for i in 1:nx, j in 1:ny]
  discrete2 = DiscreteGrid(
    x2, y2, 2; interpolation=:linear, compute_metrics=false, cache_mode=:off
  )

  s2 = SurfaceGrid(discrete2, :ymin, true)
  @test s2.axis == 2
  @test s2.side == :lo
  @test minimum(s2.face_areas) > 0.0
  @test all(
    I -> isapprox(hypot(s2.face_normals[1][I], s2.face_normals[2][I]), 1.0; atol=1e-10),
    CartesianIndices(s2.face_areas),
  )

  nx3, ny3, nz3 = (8, 9, 10)
  x3 = [i + 0.1 * j for i in 1:nx3, j in 1:ny3, k in 1:nz3]
  y3 = [j + 0.2 * k for i in 1:nx3, j in 1:ny3, k in 1:nz3]
  z3 = [k + 0.3 * i for i in 1:nx3, j in 1:ny3, k in 1:nz3]
  discrete3 = DiscreteGrid(
    x3, y3, z3, 2; interpolation=:linear, compute_metrics=false, cache_mode=:off
  )
  out = SurfaceGrid(discrete3, :ilo, true)
  inn = SurfaceGrid(discrete3, :ilo, false)
  @test minimum(out.face_areas) > 0.0
  @test all(
    I -> begin
      all(isapprox(out.face_normals[d][I], -inn.face_normals[d][I]; atol=1e-12) for d in 1:3)
    end,
    CartesianIndices(out.face_areas),
  )

  out3 = extract_surface_mesh(discrete3, :ilo)
  @test out3 isa SurfaceGrid{3,Float64}

  x0(t, ξ, η, p) = zero(ξ)
  y0(t, ξ, η, p) = zero(η)
  degenerate2 = MappedGrid(x0, y0, (;), (4, 4), 1; compute_metrics=false, cache_mode=:off)
  sdeg2 = SurfaceGrid(degenerate2, :ilo, true)
  I = first(CartesianIndices(sdeg2.face_areas))
  @test sdeg2.face_areas[I] == 0.0
  @test sdeg2.face_normals[1][I] == 0.0
  @test sdeg2.face_normals[2][I] == 0.0
end

@testset "Legacy extract_surface_mesh compatibility" begin
  x = collect(range(0.0, 1.0; length=8))
  y = collect(range(-1.0, 1.0; length=9))
  z = collect(range(2.0, 3.0; length=7))

  legacy2 = RectilinearGrid2D(x, y, :MEG6)
  sx, sy = extract_surface_mesh(legacy2, :ihi)
  @test size(sx) == (1, length(y))
  @test size(sy) == (1, length(y))
  @test_throws ArgumentError extract_surface_mesh(legacy2, :klo)

  legacy3 = RectilinearGrid3D(x, y, z, :MEG6)
  sx3, sy3, sz3 = extract_surface_mesh(legacy3, :khi)
  @test size(sx3) == (length(x), length(y), 1)
  @test size(sy3) == (length(x), length(y), 1)
  @test size(sz3) == (length(x), length(y), 1)
  @test_throws ArgumentError extract_surface_mesh(legacy3, :qhi)
end

@testset "SurfaceGrid VTK output fields" begin
  x2(t, ξ, η, p) = ξ + 0.5 * sin(0.25 * η)
  y2(t, ξ, η, p) = η + 0.10 * sin(0.20 * ξ) * cos(0.30 * η)
  x3(t, ξ, η, ζ, p) = ξ + 0.8 * sin(0.20 * η) * cos(0.25 * ζ)
  y3(t, ξ, η, ζ, p) = η + 0.7 * sin(0.30 * ζ) * cos(0.20 * ξ)
  z3(t, ξ, η, ζ, p) = ζ + 0.6 * sin(0.25 * ξ) * cos(0.20 * η)

  mapped2 = MappedGrid(x2, y2, (;), (8, 9), 2; compute_metrics=false, cache_mode=:off)
  mapped3 = MappedGrid(
    x3, y3, z3, (;), (6, 7, 8), 2; compute_metrics=false, cache_mode=:off
  )
  surf2 = SurfaceGrid(mapped2, :ilo, true)
  surf3 = SurfaceGrid(mapped3, :khi, true)

  @test !all(isapprox.(surf2.node_coordinates[1], surf2.node_coordinates[1][1]; atol=1e-12))
  @test !all(
    isapprox.(surf3.node_coordinates[3], surf3.node_coordinates[3][1, 1]; atol=1e-12)
  )

  dir = joinpath(dirname(@__DIR__), "vtk_surface_outputs")
  mkpath(dir)

  fn2 = joinpath(dir, "surface2d_curved")
  for ext in (".vti", ".vts")
    isfile(fn2 * ext) && rm(fn2 * ext)
  end
  save_vtk(surf2, fn2)
  vtk2 = isfile(fn2 * ".vts") ? (fn2 * ".vts") : (fn2 * ".vti")
  @test isfile(vtk2)
  txt2 = read(vtk2, String)
  @test occursin("Points", txt2)
  @test occursin("face_area", txt2)
  @test occursin("face_normal", txt2)
  @test !occursin("face_center", txt2)
  @test !occursin("parent_cell_index", txt2)

  fn3 = joinpath(dir, "surface3d_curved")
  for ext in (".vti", ".vts")
    isfile(fn3 * ext) && rm(fn3 * ext)
  end
  save_vtk(surf3, fn3)
  vtk3 = isfile(fn3 * ".vts") ? (fn3 * ".vts") : (fn3 * ".vti")
  @test isfile(vtk3)
  txt3 = read(vtk3, String)
  @test occursin("Points", txt3)
  @test occursin("face_area", txt3)
  @test occursin("face_normal", txt3)
  @test !occursin("face_center", txt3)
  @test !occursin("parent_cell_index", txt3)
  # end
end
