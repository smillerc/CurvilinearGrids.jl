using Test
using LinearAlgebra
using StaticArrays
using KernelAbstractions
using CurvilinearGrids

function _rect_nodes_2d(x0, x1, y0, y1, ni, nj)
  xnodes = collect(range(x0, x1; length=ni + 1))
  ynodes = collect(range(y0, y1; length=nj + 1))
  x = [xnodes[i] for i in eachindex(xnodes), j in eachindex(ynodes)]
  y = [ynodes[j] for i in eachindex(xnodes), j in eachindex(ynodes)]
  return x, y
end

function _abutting_two_block_setup(; nhalo=2, gap=0.0)
  ni, nj = (6, 5)
  x1, y1 = _rect_nodes_2d(0.0, 1.0, 0.0, 1.0, ni, nj)
  x2, y2 = _rect_nodes_2d(1.0 + gap, 2.0 + gap, 0.0, 1.0, ni, nj)
  b1 = DiscreteGrid(x1, y1, nhalo)
  b2 = DiscreteGrid(x2, y2, nhalo)
  iface = (b1, :ihi) => (b2, :ilo)
  return b1, b2, iface
end

@inline _linear_coord(s::Real, n::Int, lo::Real, hi::Real) = lo + ((s - 1) / n) * (hi - lo)

function _abutting_mapped_3d_cartesian_spherical_setup(; nhalo=2)
  ni, nj, nk = (5, 4, 4)

  params_cart = (; ni, nj, nk, xlo=1.0, xhi=2.0, ylo=0.6, yhi=1.0, zlo=0.4, zhi=0.9)
  x1c(t, ξ, η, ζ, p) = _linear_coord(ξ, p.ni, p.xlo, p.xhi)
  x2c(t, ξ, η, ζ, p) = _linear_coord(η, p.nj, p.ylo, p.yhi)
  x3c(t, ξ, η, ζ, p) = _linear_coord(ζ, p.nk, p.zlo, p.zhi)

  params_sph = (; ni, nj, nk, xlo=2.0, xhi=3.0, ylo=0.6, yhi=1.0, zlo=0.4, zhi=0.9)
  x_phys(t, ξ, η, ζ, p) = _linear_coord(ξ, p.ni, p.xlo, p.xhi)
  y_phys(t, ξ, η, ζ, p) = _linear_coord(η, p.nj, p.ylo, p.yhi)
  z_phys(t, ξ, η, ζ, p) = _linear_coord(ζ, p.nk, p.zlo, p.zhi)
  rmap(t, ξ, η, ζ, p) = sqrt(
    x_phys(t, ξ, η, ζ, p)^2 + y_phys(t, ξ, η, ζ, p)^2 + z_phys(t, ξ, η, ζ, p)^2
  )
  θmap(t, ξ, η, ζ, p) = begin
    r = rmap(t, ξ, η, ζ, p)
    acos(z_phys(t, ξ, η, ζ, p) / r)
  end
  ϕmap(t, ξ, η, ζ, p) = atan(y_phys(t, ξ, η, ζ, p), x_phys(t, ξ, η, ζ, p))

  g_cart = MappedGrid(
    x1c,
    x2c,
    x3c,
    params_cart,
    (ni, nj, nk),
    nhalo;
    coordinate_system=CartesianCS(),
    basis=CartesianBasis(),
  )
  g_sph = MappedGrid(
    rmap,
    θmap,
    ϕmap,
    params_sph,
    (ni, nj, nk),
    nhalo;
    coordinate_system=SphericalCS(),
    basis=SphericalBasis(),
  )

  iface = (g_cart, :ihi) => (g_sph, :ilo)
  return g_cart, g_sph, iface
end

function _abutting_orthogonal_discrete_setup(; nhalo=1)
  backend = CPU()
  ni, nj = (6, 5)
  left = OrthogonalGrid(
    CartesianOrthogonalGrid2D(
      collect(range(0.0, 1.0; length=ni + 1)),
      collect(range(0.0, 1.0; length=nj + 1)),
      nhalo,
      backend,
    ),
  )
  xnodes, ynodes = _rect_nodes_2d(1.0, 2.0, 0.0, 1.0, ni, nj)
  right = DiscreteGrid(xnodes, ynodes, nhalo)
  iface = (left, :ihi) => (right, :ilo)
  return left, right, iface
end

@inline function _spherical_basis_to_cartesian_matrix(q::SVector{3,<:Real})
  _, θ, ϕ = q
  sθ, cθ = sin(θ), cos(θ)
  sϕ, cϕ = sin(ϕ), cos(ϕ)
  return @SMatrix [
    sθ * cϕ cθ * cϕ -sϕ
    sθ * sϕ cθ * sϕ cϕ
    cθ -sθ 0.0
  ]
end

@testset "MultiBlock strict topology and abutting validation" begin
  b1, b2, iface = _abutting_two_block_setup()

  mb = MultiBlockMesh((b1, b2), (iface,); tolerance=1e-12)
  @test length(mb.interfaces) == 1

  # Interfaces must be provided as a tuple.
  @test_throws ArgumentError MultiBlockMesh((b1, b2), iface; tolerance=1e-12)

  # Duplicate interface face assignment must fail.
  @test_throws ArgumentError MultiBlockMesh((b1, b2), (iface, iface); tolerance=1e-12)

  # Invalid face symbol must fail.
  bad_iface = (b1, :foo) => (b2, :ilo)
  @test_throws ArgumentError MultiBlockMesh((b1, b2), (bad_iface,))

  # Non-abutting interface must fail.
  b1g, b2g, ifaceg = _abutting_two_block_setup(; gap=0.1)
  @test_throws ArgumentError MultiBlockMesh((b1g, b2g), (ifaceg,); tolerance=1e-12)
end

@testset "MultiBlock scalar/vector/tensor exchange" begin
  b1, b2, iface = _abutting_two_block_setup()
  mb = MultiBlockMesh((b1, b2), (iface,); tolerance=1e-12)

  f1 = zeros(Float64, size(b1.iterators.cell.full))
  f2 = zeros(Float64, size(b2.iterators.cell.full))

  i1_interior = last(b1.iterators.cell.domain.indices[1])
  i1_ghost = i1_interior + 1
  i2_interior = first(b2.iterators.cell.domain.indices[1])
  i2_ghost = i2_interior - 1

  for j in b1.iterators.cell.domain.indices[2]
    f1[i1_interior, j] = 3.0
    f2[i2_interior, j] = 7.0
  end

  exchange_interface!(mb, 1, [f1, f2]; field_kind=:scalar)
  for j in b1.iterators.cell.domain.indices[2]
    @test f1[i1_ghost, j] == 7.0
    @test f2[i2_ghost, j] == 3.0
  end

  f1_depth = zeros(Float64, size(b1.iterators.cell.full))
  f2_depth = zeros(Float64, size(b2.iterators.cell.full))
  i1_interior_2 = i1_interior - 1
  i1_ghost_2 = i1_interior + 2
  i2_interior_2 = i2_interior + 1
  i2_ghost_2 = i2_interior - 2
  for j in b1.iterators.cell.domain.indices[2]
    f1_depth[i1_interior, j] = 31.0
    f1_depth[i1_interior_2, j] = 32.0
    f2_depth[i2_interior, j] = 71.0
    f2_depth[i2_interior_2, j] = 72.0
  end
  exchange_interface!(mb, 1, [f1_depth, f2_depth]; field_kind=:scalar, depth=2)
  for j in b1.iterators.cell.domain.indices[2]
    @test f1_depth[i1_ghost, j] == 71.0
    @test f1_depth[i1_ghost_2, j] == 72.0
    @test f2_depth[i2_ghost, j] == 31.0
    @test f2_depth[i2_ghost_2, j] == 32.0
  end
  @test_throws ArgumentError exchange_interface!(
    mb, 1, [f1_depth, f2_depth]; field_kind=:scalar, depth=3
  )

  v1 = fill((@SVector [0.0, 0.0]), size(b1.iterators.cell.full))
  v2 = fill((@SVector [0.0, 0.0]), size(b2.iterators.cell.full))
  for j in b1.iterators.cell.domain.indices[2]
    v1[i1_interior, j] = @SVector [1.0, 2.0]
    v2[i2_interior, j] = @SVector [3.0, 4.0]
  end
  exchange_interface!(mb, 1, [v1, v2]; field_kind=:vector)
  for j in b1.iterators.cell.domain.indices[2]
    @test v1[i1_ghost, j] ≈ @SVector [3.0, 4.0]
    @test v2[i2_ghost, j] ≈ @SVector [1.0, 2.0]
  end

  t1 = fill((@SMatrix [0.0 0.0; 0.0 0.0]), size(b1.iterators.cell.full))
  t2 = fill((@SMatrix [0.0 0.0; 0.0 0.0]), size(b2.iterators.cell.full))
  for j in b1.iterators.cell.domain.indices[2]
    t1[i1_interior, j] = @SMatrix [1.0 2.0; 3.0 4.0]
    t2[i2_interior, j] = @SMatrix [5.0 6.0; 7.0 8.0]
  end
  exchange_interface!(mb, 1, [t1, t2]; field_kind=:tensor)
  for j in b1.iterators.cell.domain.indices[2]
    @test t1[i1_ghost, j] ≈ @SMatrix [5.0 6.0; 7.0 8.0]
    @test t2[i2_ghost, j] ≈ @SMatrix [1.0 2.0; 3.0 4.0]
  end
end

@testset "MultiBlock computational coordinate transfer" begin
  b1, b2, iface = _abutting_two_block_setup()
  mb = MultiBlockMesh((b1, b2), (iface,); tolerance=1e-12)

  i1_face = last(b1.iterators.cell.domain.indices[1]) + 0.5
  i2_face = first(b2.iterators.cell.domain.indices[1]) - 0.5
  ξ_left = (i1_face, 4.25)

  ξ_right = computational_coordinate(mb, 1, ξ_left; from=:left)
  @test ξ_right[1] ≈ i2_face
  @test ξ_right[2] ≈ ξ_left[2]

  back = computational_coordinate(mb, 1, ξ_right; from=:right)
  @test back[1] ≈ i1_face
  @test back[2] ≈ ξ_left[2]

  transferred = computational_coordinate(mb, 1, 1, ξ_left)
  @test transferred.block_id == 2
  @test transferred.coordinate ≈ ξ_right
end

@testset "MultiBlock orthogonal/discrete validation and exchange" begin
  left, right, iface = _abutting_orthogonal_discrete_setup()
  mb = MultiBlockMesh((left, right), (iface,); tolerance=1e-12)

  f_left = zeros(Float64, size(left.iterators.cell.full))
  f_right = zeros(Float64, size(right.iterators.cell.full))
  i_left_interior = last(left.iterators.cell.domain.indices[1])
  i_left_ghost = i_left_interior + 1
  i_right_interior = first(right.iterators.cell.domain.indices[1])
  i_right_ghost = i_right_interior - 1

  for j in left.iterators.cell.domain.indices[2]
    f_left[i_left_interior, j] = 2.5
    f_right[i_right_interior, j] = 6.5
  end

  exchange_interface!(mb, 1, [f_left, f_right]; field_kind=:scalar)
  for j in left.iterators.cell.domain.indices[2]
    @test f_left[i_left_ghost, j] == 6.5
    @test f_right[i_right_ghost, j] == 2.5
  end

  v_left = fill((@SVector [0.0, 0.0]), size(left.iterators.cell.full))
  v_right = fill((@SVector [0.0, 0.0]), size(right.iterators.cell.full))
  for j in left.iterators.cell.domain.indices[2]
    v_left[i_left_interior, j] = @SVector [1.0, -2.0]
    v_right[i_right_interior, j] = @SVector [3.0, 4.0]
  end

  exchange_interface!(mb, 1, [v_left, v_right]; field_kind=:vector)
  for j in left.iterators.cell.domain.indices[2]
    @test v_left[i_left_ghost, j] ≈ @SVector [3.0, 4.0]
    @test v_right[i_right_ghost, j] ≈ @SVector [1.0, -2.0]
  end
end

if QUICK_TESTS
  @info "Skipping 3D multiblock mapped vector transfer (--quick)"
else
  @testset "MultiBlock 3D mapped vector transfer: Cartesian <-> Spherical basis" begin
    g_cart, g_sph, iface = _abutting_mapped_3d_cartesian_spherical_setup()
    mb = MultiBlockMesh((g_cart, g_sph), (iface,); tolerance=1e-10)

    v_cart = fill((@SVector [0.0, 0.0, 0.0]), size(g_cart.iterators.cell.full))
    v_sph = fill((@SVector [0.0, 0.0, 0.0]), size(g_sph.iterators.cell.full))

    i_cart_interior = last(g_cart.iterators.cell.domain.indices[1])
    i_sph_interior = first(g_sph.iterators.cell.domain.indices[1])
    i_sph_ghost = i_sph_interior - 1

    donor_vec = @SVector [1.3, -0.7, 2.1]
    for j in g_cart.iterators.cell.domain.indices[2],
      k in g_cart.iterators.cell.domain.indices[3]

      v_cart[i_cart_interior, j, k] = donor_vec
    end

    exchange_interface!(
      mb, 1, [v_cart, v_sph]; field_kind=:vector, direction=:left_to_right
    )

    for j in g_cart.iterators.cell.domain.indices[2],
      k in g_cart.iterators.cell.domain.indices[3]

      transferred_components_sph = v_sph[i_sph_ghost, j, k]
      q_sph = centroid(g_sph, (i_sph_interior, j, k))
      Q = _spherical_basis_to_cartesian_matrix(q_sph)
      transferred_cart = Q * transferred_components_sph

      mag_donor = norm(donor_vec)
      mag_transferred = norm(transferred_cart)
      @test isapprox(mag_transferred, mag_donor; rtol=1e-12, atol=1e-12)

      cos_angle = dot(transferred_cart, donor_vec) / (mag_transferred * mag_donor)
      @test isapprox(cos_angle, 1.0; atol=1.0e-12)
    end
  end
end
