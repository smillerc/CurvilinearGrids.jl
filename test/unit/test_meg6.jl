
using Test

@testset "MEG6 Gradients + Edges" begin
  using ForwardDiff
  using LinearAlgebra

  # include("../../src/metric_schemes/indexing_fun.jl")
  # include("../../src/metric_schemes/6th_order_explicit/gradients.jl")
  # include("../../src/metric_schemes/6th_order_explicit/interpolation.jl")
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
  # fhalf
end

@testset "Gradients" begin
  include("common.jl")
  # include("../../src/metric_schemes/indexing_fun.jl")
  # include("../../src/metric_schemes/6th_order_explicit/gradients.jl")
  # include("../../src/metric_schemes/6th_order_explicit/interpolation.jl")

  # include("../../src/metric_schemes/6th_order_explicit/MonotoneExplicit6thOrderScheme.jl")
  using CurvilinearGrids.MetricDiscretizationSchemes.MonotoneExplicit6thOrderScheme:
    ∂!, ∂²!, ∂x∂ξ!, toedge!, conserved_metric!

  function rect_grid(nx, ny, nz)
    x0, x1 = (0.0, 2.0)
    y0, y1 = (1, 3)
    z0, z1 = (-1, 2)

    x(ξ, η, ζ) = @. x0 + (x1 - x0) * ((ξ - 1) / (nx - 1))
    y(ξ, η, ζ) = @. y0 + (y1 - y0) * ((η - 1) / (ny - 1))
    z(ξ, η, ζ) = @. z0 + (z1 - z0) * ((ζ - 1) / (nz - 1))

    return (x, y, z)
  end

  ni, nj, nk = (5, 9, 13)
  nhalo = 6
  x, y, z = rect_grid(ni, nj, nk)

  mesh = CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)
  full_domain = mesh.iterators.cell.domain

  # @show mesh.iterators.cell.full
  meg6 = mesh.discretization_scheme

  ∂x_∂ξ = meg6.cache.metric
  xᵢ₊½ = meg6.cache.xᵢ₊½
  ∂x = meg6.cache.∂x
  ∂²x = meg6.cache.∂²x

  x = mesh.centroids.x
  y = mesh.centroids.y
  z = mesh.centroids.z

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
        #   @show extrema(∂²x[full_domain])
        @error "∂²$(dim)/∂$(n) failed"
        @test passes
        #   # break
      end
    end
  end

  # toedge!(xᵢ₊½, ∂²x, ∂x, x, full_domain, 3)
  # display(xᵢ₊½)

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
  ξ̂x = meg6.cache.metric
  conserved_metric!(meg6, ξ̂x, y, η, z, ζ, full_domain)

  # expanded_dom = expand(mesh.iterators.cell.domain, +1)
  # @show ξ̂x[expanded_dom]
  @test all(ξ̂x[mesh.iterators.cell.domain] .≈ 0.0625)

  # toedge!(xᵢ₊½, ∂²x, ∂x, ξ̂x, expand(full_domain, ξ, -1), ξ)
  # @test all(xᵢ₊½[mesh.iterators.cell.domain] .≈ 0.0625)

  # toedge!(xᵢ₊½, ∂²x, ∂x, ξ̂x, expand(full_domain, η, -1), η)
  # @test all(xᵢ₊½[mesh.iterators.cell.domain] .≈ 0.0625)

  # toedge!(xᵢ₊½, ∂²x, ∂x, ξ̂x, expand(full_domain, ζ, -1), ζ)
  # @test all(xᵢ₊½[mesh.iterators.cell.domain] .≈ 0.0625)

  # # @show xᵢ₊½
  # # toedge!(xᵢ₊½, ∂²x, ∂x, ξ̂x, full_domain, η)
  # # toedge!(xᵢ₊½, ∂²x, ∂x, ξ̂x, full_domain, ζ)

  # η̂y = meg6.cache.metric
  # conserved_metric!(meg6, η̂y, z, ζ, x, ξ, full_domain)
  # @test all(η̂y[mesh.iterators.cell.domain] .≈ 0.125)

  # ζ̂z = meg6.cache.metric
  # conserved_metric!(meg6, ζ̂z, x, ξ, y, η, full_domain)
  # @test all(ζ̂z[mesh.iterators.cell.domain] .≈ 0.125)

  # @show extrema(ξ̂x[mesh.iterators.cell.domain])
  # ∂!(∂x, x, domain, 3)
end
