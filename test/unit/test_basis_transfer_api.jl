@inline function _spherical_basis_to_cartesian_matrix_3d(q::SVector{3,<:Real})
  _, θ, ϕ = q
  sθ, cθ = sin(θ), cos(θ)
  sϕ, cϕ = sin(ϕ), cos(ϕ)
  return @SMatrix [
    sθ * cϕ cθ * cϕ -sϕ
    sθ * sϕ cθ * sϕ cϕ
    cθ -sθ 0.0
  ]
end

@testset "Public basis-transfer API" begin
  params2 = (; ni=8, nj=7, x0=0.0, x1=1.0, y0=-0.2, y1=0.9)
  xmap2(t, ξ, η, p) = p.x0 + ((ξ - 1) / p.ni) * (p.x1 - p.x0)
  ymap2(t, ξ, η, p) = p.y0 + ((η - 1) / p.nj) * (p.y1 - p.y0)
  cart2 = MappedGrid(
    xmap2,
    ymap2,
    params2,
    (params2.ni, params2.nj),
    2;
    coordinate_system=CartesianCS(),
    basis=CartesianBasis(),
    compute_metrics=false,
    cache_mode=:off,
  )

  Rcart = basis_transfer_matrix(cart2, (2.25, 3.5), (6.75, 5.25))
  expected_cart = @SMatrix [1.0 0.0; 0.0 1.0]
  @test Rcart ≈ expected_cart

  backend = CPU()
  sph2 = SphericalOrthogonalGrid2D([1.0, 1.5, 2.0], [0.2, 0.6, 1.0], 0, backend)
  θa = 0.35
  θb = 0.9
  Δθ = θb - θa
  Rsph2 = basis_transfer_matrix(sph2, SVector(1.2, θa), SVector(1.2, θb))
  expected_sph2 = @SMatrix [cos(Δθ) sin(Δθ); -sin(Δθ) cos(Δθ)]
  @test Rsph2 ≈ expected_sph2 atol=1.0e-12

  params3 = (;
    ni=5,
    nj=4,
    nk=3,
    r0=1.4,
    θ0=0.35,
    ϕ0=0.2,
    Δr=0.1,
    Δθ=0.08,
    Δϕ=0.09,
  )
  rmap(t, ξ, η, ζ, p) = p.r0 + (ξ - 1) * p.Δr
  θmap(t, ξ, η, ζ, p) = p.θ0 + (η - 1) * p.Δθ
  ϕmap(t, ξ, η, ζ, p) = p.ϕ0 + (ζ - 1) * p.Δϕ
  xmap3(t, ξ, η, ζ, p) = p.r0 + (ξ - 1) * p.Δr
  ymap3(t, ξ, η, ζ, p) = p.θ0 + (η - 1) * p.Δθ
  zmap3(t, ξ, η, ζ, p) = p.ϕ0 + (ζ - 1) * p.Δϕ
  cart3 = MappedGrid(
    xmap3,
    ymap3,
    zmap3,
    params3,
    (params3.ni, params3.nj, params3.nk),
    2;
    coordinate_system=CartesianCS(),
    basis=CartesianBasis(),
    compute_metrics=false,
    cache_mode=:off,
  )
  sph3 = MappedGrid(
    rmap,
    θmap,
    ϕmap,
    params3,
    (params3.ni, params3.nj, params3.nk),
    2;
    coordinate_system=SphericalCS(),
    basis=SphericalBasis(),
    compute_metrics=false,
    cache_mode=:off,
  )
  qa = SVector(1.6, 0.55, 0.4)
  qb = SVector(1.8, 0.82, 0.73)
  Qa = _spherical_basis_to_cartesian_matrix_3d(SVector(qa[1], qa[2], qa[3]))
  Qb = _spherical_basis_to_cartesian_matrix_3d(SVector(qb[1], qb[2], qb[3]))
  Rsph3 = basis_transfer_matrix(sph3, qa, qb)
  @test Rsph3 ≈ transpose(Qb) * Qa atol=1.0e-12
  @test Rsph3 * transpose(Rsph3) ≈ I atol=1.0e-12

  qcart = SVector(1.6, 0.8, 0.9)
  rs = norm(qcart)
  qsph = SVector(rs, acos(qcart[3] / rs), atan(qcart[2], qcart[1]))
  Qsph = _spherical_basis_to_cartesian_matrix_3d(qsph)
  Rcart_to_sph = basis_transfer_matrix(cart3, qcart, sph3, qsph)
  Rsph_to_cart = basis_transfer_matrix(sph3, qsph, cart3, qcart)
  @test Rcart_to_sph ≈ transpose(Qsph) atol=1.0e-12
  @test Rsph_to_cart ≈ Qsph atol=1.0e-12

  bad2 = MappedGrid(
    xmap2,
    ymap2,
    params2,
    (params2.ni, params2.nj),
    2;
    coordinate_system=CylindricalCS(),
    basis=SphericalBasis(),
    compute_metrics=false,
    cache_mode=:off,
  )
  @test_throws ArgumentError basis_transfer_matrix(bad2, SVector(1.0, 0.2), SVector(1.1, 0.3))
end

@testset "Public face-coordinate and face-flux geometry API" begin
  backend = CPU()

  cart_orth = CartesianOrthogonalGrid2D([0.0, 1.0, 3.0], [0.0, 2.0, 5.0], 0, backend)
  @test centroid(cart_orth, (1, 2)) == @SVector [0.5, 3.5]
  @test face_coordinate(cart_orth, (1, 2), :ihi) == @SVector [1.0, 3.5]
  @test face_coordinate(cart_orth, (1, 2), :ilo) == @SVector [0.0, 3.5]
  @test face_coordinate(cart_orth, (1, 2), :jhi) == @SVector [0.5, 5.0]

  axi_orth = AxisymmetricOrthogonalGrid2D([1.0, 2.0, 4.0], [0.0, 1.5, 3.0], 0, backend)
  @test centroid(axi_orth, (1, 2)) ≈ (@SVector [14 / 9, 2.25]) atol=1.0e-12
  @test face_coordinate(axi_orth, (1, 2), :ihi) ≈ (@SVector [2.0, 2.25]) atol=1.0e-12
  @test face_coordinate(axi_orth, (1, 2), :jlo) ≈ (@SVector [14 / 9, 1.5]) atol=1.0e-12

  sph_orth = SphericalOrthogonalGrid2D([1.0, 2.0], [0.0, π / 3, 2π / 3], 0, backend)
  @test centroid(sph_orth, (1, 1)) ≈ (@SVector [45 / 28, acos(3 / 4)]) atol=1.0e-12
  @test face_coordinate(sph_orth, (1, 1), :ihi) ≈ (@SVector [2.0, acos(3 / 4)]) atol=1.0e-12
  @test face_coordinate(sph_orth, (1, 1), :jhi) ≈ (@SVector [45 / 28, π / 3]) atol=1.0e-12
  @test_throws ArgumentError face_flux_geometry(sph_orth, (1, 1), :ihi)

  params2 = (; ni=9, nj=8, x0=-0.5, x1=1.3, y0=0.2, y1=1.1, warp=0.04)
  xmap2(t, ξ, η, p) =
    p.x0 + ((ξ - 1) / p.ni) * (p.x1 - p.x0) + p.warp * sin(2π * (η - 1) / p.nj)
  ymap2(t, ξ, η, p) =
    p.y0 + ((η - 1) / p.nj) * (p.y1 - p.y0) + 0.3 * p.warp * cos(2π * (ξ - 1) / p.ni)
  cart2 = MappedGrid(
    xmap2,
    ymap2,
    params2,
    (params2.ni, params2.nj),
    2;
    coordinate_system=CartesianCS(),
    basis=CartesianBasis(),
    cache_mode=:eager,
  )

  I2 = (cart2.nhalo + 3, cart2.nhalo + 4)
  @test face_coordinate(cart2, I2, :ilo) ≈
    face_coordinate(cart2, (I2[1] - 1, I2[2]), :ihi) atol=1.0e-12
  @test face_coordinate(cart2, I2, :jlo) ≈
    face_coordinate(cart2, (I2[1], I2[2] - 1), :jhi) atol=1.0e-12

  geom_hi = face_flux_geometry(cart2, I2, :ihi)
  geom_lo = face_flux_geometry(cart2, I2, :ilo)
  Ghi = face_metrics(cart2)[1].conserved[I2...].jacobian_matrix
  Glo = face_metrics(cart2)[1].conserved[I2[1] - 1, I2[2]].jacobian_matrix
  @test geom_hi.coordinate ≈ face_coordinate(cart2, I2, :ihi) atol=1.0e-12
  @test geom_hi.metric_vector ≈ (@SVector [Ghi[1, 1], Ghi[1, 2]]) atol=1.0e-12
  @test geom_hi.area ≈ norm(geom_hi.metric_vector) atol=1.0e-12
  @test geom_hi.normal * geom_hi.area ≈ geom_hi.metric_vector atol=1.0e-12
  @test geom_lo.metric_vector ≈ (-@SVector [Glo[1, 1], Glo[1, 2]]) atol=1.0e-12

  params3 = (; ni=6, nj=5, nk=4, r0=1.4, θ0=0.45, ϕ0=0.3, Δr=0.08, Δθ=0.07, Δϕ=0.06)
  rmap(t, ξ, η, ζ, p) = p.r0 + (ξ - 1) * p.Δr
  θmap(t, ξ, η, ζ, p) = p.θ0 + (η - 1) * p.Δθ
  ϕmap(t, ξ, η, ζ, p) = p.ϕ0 + (ζ - 1) * p.Δϕ
  sph3 = MappedGrid(
    rmap,
    θmap,
    ϕmap,
    params3,
    (params3.ni, params3.nj, params3.nk),
    2;
    coordinate_system=SphericalCS(),
    basis=SphericalBasis(),
    cache_mode=:eager,
  )

  I3 = (sph3.nhalo + 2, sph3.nhalo + 2, sph3.nhalo + 2)
  geom3 = face_flux_geometry(sph3, I3, :ihi)
  outward3 = outward_face_normal(sph3, I3, :ihi)
  q3 = face_coordinate(sph3, I3, :ihi)
  Q3 = _spherical_basis_to_cartesian_matrix_3d(SVector(q3[1], q3[2], q3[3]))
  @test geom3.normal ≈ (@SVector [1.0, 0.0, 0.0]) atol=1.0e-12
  @test outward3 ≈ SVector(Q3[1, 1], Q3[2, 1], Q3[3, 1]) atol=1.0e-12
  @test norm(outward3 - geom3.normal) > 1.0e-3
end
