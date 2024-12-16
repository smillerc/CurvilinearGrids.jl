
@testset "MEG6 Gradients + Edges" begin
  using ForwardDiff
  using LinearAlgebra

  using CurvilinearGrids.MetricDiscretizationSchemes.MonotoneExplicit6thOrderScheme:
    ∂!, ∂²!, toedge!

  n_i = 20
  f(i) = sinpi(0.5((i - 1) / (n_i - 1)))
  ∂f(i) = ForwardDiff.derivative(f, i)
  ∂²f(i) = ForwardDiff.derivative(∂f, i)

  xn = -6:31
  xc = xn[2:end] .- 0.5

  fdata = f.(xc)
  ∂fdata = similar(fdata)
  ∂²fdata = similar(fdata)

  domain = CartesianIndices(fdata)
  inner_domain = expand(domain, -6)

  ∂!(∂fdata, fdata, domain, 1)
  ∂²!(∂²fdata, ∂fdata, fdata, domain, 1)

  ∂f_err = @. abs(∂f(xc) - ∂fdata)
  ∂²f_err = @. abs(∂²f(xc) - ∂²fdata)

  L2_∂f_err = mapreduce(√, +, ∂f_err)
  L2_∂²f_err = mapreduce(√, +, ∂²f_err)

  # The inner domain is more accurate
  @test norm(∂f_err[inner_domain]) < 1e-6
  @test norm(∂²f_err[inner_domain]) < 1e-6

  # Lower order boundaries
  @test norm(∂f_err) < 1e-3
  @test norm(∂²f_err) < 1e-3

  inner_L2_∂f_err = mapreduce(√, +, ∂f_err[inner_domain])
  inner_L2_∂²f_err = mapreduce(√, +, ∂²f_err[inner_domain])

  # Interpolate to the edges
  fhalf = similar(fdata)
  fill!(fhalf, NaN)
  toedge!(fhalf, ∂²fdata, ∂fdata, fdata, inner_domain, 1)

  fhalf_err = @. f(xc + 0.5) - fhalf
  @test norm(fhalf_err[inner_domain]) < 1e-2

  # Turn plots on for diagnosis if need be
  # using Plots
  # plot(∂f_err[inner_domain]; label="∂f L2 = $L2_∂f_err")

  # plot(xc, fdata; label="fdata")
  # plot!(xc .+ 0.5, fhalf; label="fᵢ₊½", marker=:circle, ms=1.5)

  # plot!(xc, ∂f.(xc); label="AD ∂f")
  # plot!(xc, ∂²f.(xc); label="AD ∂²f")
  # plot!(xc, ∂fdata; ms=2, marker=:square, label="MEG6 ∂f")
  # plot!(xc, ∂²fdata; ms=2, marker=:square, label="MEG6 ∂²f")
end

@testset "MEG6 Gradients + Edges + CurvilinearGrid3D" begin
  # include("common.jl")

  using CurvilinearGrids.MetricDiscretizationSchemes.MonotoneExplicit6thOrderScheme:
    ∂!, ∂²!, ∂x∂ξ!, toedge!, conserved_metric!

  # include("common.jl")

  x0, x1 = (0.0, 2.0)
  y0, y1 = (1, 3)
  z0, z1 = (-1, 2)
  ni, nj, nk = (4, 8, 12)
  nhalo = 4

  mesh = rectlinear_grid((x0, y0, z0), (x1, y1, z1), (ni, nj, nk), nhalo)
  full_domain = mesh.iterators.cell.domain

  meg6 = mesh.discretization_scheme

  ∂x_∂ξ = meg6.metric
  xᵢ₊½ = meg6.xᵢ₊½
  ∂x = meg6.∂x
  ∂²x = meg6.∂²x

  x = mesh.centroid_coordinates.x
  y = mesh.centroid_coordinates.y
  z = mesh.centroid_coordinates.z

  ∂!(∂x, x, full_domain, 1)
  @test all(∂x[full_domain] .≈ 0.5)
  ∂!(∂x, x, full_domain, 2)
  @test all(iszero.(∂x[full_domain]))
  ∂!(∂x, x, full_domain, 3)
  @test all(iszero.(∂x[full_domain]))

  ∂!(∂x, y, full_domain, 1)
  @test all(iszero.(∂x[full_domain]))
  ∂!(∂x, y, full_domain, 2)
  @test all(∂x[full_domain] .≈ 0.25)
  ∂!(∂x, y, full_domain, 3)
  @test all(iszero.(∂x[full_domain]))

  ∂!(∂x, z, full_domain, 1)
  @test all(iszero.(∂x[full_domain]))
  ∂!(∂x, z, full_domain, 2)
  @test all(iszero.(∂x[full_domain]))
  ∂!(∂x, z, full_domain, 3)
  @test all(∂x[full_domain] .≈ 0.25)

  for (dim, _x) in zip((:x, :y, :z), (x, y, z))
    for (n, ax) in zip((:ξ, :η, :ζ), (1, 2, 3))
      # println("∂²$(dim)/∂$(n)")
      ∂²!(∂²x, ∂x, _x, full_domain, ax)
      passes = all(iszero.(∂²x[full_domain]))
      if !passes
        @error "∂²$(dim)/∂$(n) failed"
        @test passes
      end
    end
  end

  ∂x∂ξ!(meg6, ∂x_∂ξ, x, full_domain, 1) # xξ
  @test all(∂x_∂ξ[mesh.iterators.cell.domain] .≈ 0.5)
  ∂x∂ξ!(meg6, ∂x_∂ξ, x, full_domain, 2) # xη
  @test all(∂x_∂ξ[mesh.iterators.cell.domain] .≈ 0.0)
  ∂x∂ξ!(meg6, ∂x_∂ξ, x, full_domain, 3) # xζ
  @test all(∂x_∂ξ[mesh.iterators.cell.domain] .≈ 0.0)

  ∂x∂ξ!(meg6, ∂x_∂ξ, y, full_domain, 1) # yξ
  @test all(∂x_∂ξ[mesh.iterators.cell.domain] .≈ 0.0)
  ∂x∂ξ!(meg6, ∂x_∂ξ, y, full_domain, 2) # yη
  @test all(∂x_∂ξ[mesh.iterators.cell.domain] .≈ 0.25)
  ∂x∂ξ!(meg6, ∂x_∂ξ, y, full_domain, 3) # yζ
  @test all(∂x_∂ξ[mesh.iterators.cell.domain] .≈ 0.0)

  ∂x∂ξ!(meg6, ∂x_∂ξ, z, full_domain, 1) # zξ
  @test all(∂x_∂ξ[mesh.iterators.cell.domain] .≈ 0.0)
  ∂x∂ξ!(meg6, ∂x_∂ξ, z, full_domain, 2) # zη
  @test all(∂x_∂ξ[mesh.iterators.cell.domain] .≈ 0.0)
  ∂x∂ξ!(meg6, ∂x_∂ξ, z, full_domain, 3) # zζ
  @test all(∂x_∂ξ[mesh.iterators.cell.domain] .≈ 0.25)

  ξ, η, ζ = (1, 2, 3)
  ξ̂x = meg6.metric
  conserved_metric!(meg6, ξ̂x, y, η, z, ζ, full_domain)

  @test all(ξ̂x[mesh.iterators.cell.domain] .≈ 0.0625)
end
