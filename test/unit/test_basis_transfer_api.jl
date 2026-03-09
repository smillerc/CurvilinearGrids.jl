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
