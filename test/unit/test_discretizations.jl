using Test
using CartesianDomains, KernelAbstractions
using CurvilinearGrids
using StructArrays
using CurvilinearGrids.DiscretizationSchemes
using CurvilinearGrids.DiscretizationSchemes: central_derivative
using CurvilinearGrids.MetricDiscretizationSchemes

@testset "MonotoneExplicitGradientScheme basic constructors" begin
  second_order = SecondOrder()
  fourth_order = FourthOrder()
  sixth_order = SixthOrder()
  eighth_order = EighthOrder()

  celldims = (30, 40, 50)
  @test_nowarn begin
    meg2 = MonotoneExplicitGradientScheme(2; use_cache=false)
  end

  @test_nowarn begin
    meg4 = MonotoneExplicitGradientScheme(4; use_cache=false)
  end

  @test_nowarn begin
    meg6 = MonotoneExplicitGradientScheme(6; use_cache=false)
  end

  @test_nowarn begin
    meg2 = MonotoneExplicitGradientScheme(2; celldims=celldims, use_cache=true)
  end

  @test_nowarn begin
    meg4 = MonotoneExplicitGradientScheme(4; celldims=celldims, use_cache=true)
  end

  @test_nowarn begin
    meg6 = MonotoneExplicitGradientScheme(6; celldims=celldims, use_cache=true)
  end
end

@testset "MEG{2,4,6} 1st and 2nd derivatives" begin
  celldims = (30, 40, 50)

  function init_xyz(cell_dims)
    x = zeros(cell_dims...)
    y = zeros(cell_dims...)
    z = zeros(cell_dims...)

    for I in CartesianIndices(x)
      i, j, k = I.I
      x[i, j, k] = i
      y[i, j, k] = 2j
      z[i, j, k] = 3k
    end

    ∂x_∂ξ = zeros(cell_dims...) .* NaN
    ∂y_∂η = zeros(cell_dims...) .* NaN
    ∂z_∂ζ = zeros(cell_dims...) .* NaN

    ∂²x_∂ξ² = zeros(cell_dims...) .* NaN
    ∂²y_∂η² = zeros(cell_dims...) .* NaN
    ∂²z_∂ζ² = zeros(cell_dims...) .* NaN

    return (x, y, z, ∂x_∂ξ, ∂y_∂η, ∂z_∂ζ, ∂²x_∂ξ², ∂²y_∂η², ∂²z_∂ζ²)
  end

  x, y, z, ∂x_∂ξ, ∂y_∂η, ∂z_∂ζ, ∂²x_∂ξ², ∂²y_∂η², ∂²z_∂ζ² = init_xyz(celldims)
  iaxis, jaxis, kaxis = (1, 2, 3)
  full_domain = CartesianIndices(x)

  for order in (2, 4, 6)
    meg = MonotoneExplicitGradientScheme(order; celldims=celldims, use_cache=true)

    inner_domain = expand(full_domain, -meg.nhalo)

    ∂halo = (meg.nhalo + 1) ÷ 2
    ∂domain = expand(full_domain, -∂halo)
    ∂²domain = expand(full_domain, -∂halo - 1)
    # ∂domain = expand(full_domain, -3)
    # ∂²domain = expand(full_domain, -4)

    # @test DiscretizationSchemes.first_derivative_domain(meg, full_domain) == ∂domain
    # @test DiscretizationSchemes.second_derivative_domain(meg, ∂domain) == ∂²domain

    backend = get_backend(x)

    compute_first_derivatives!(
      meg, ∂x_∂ξ, x, iaxis, ∂domain, backend; use_one_sided_on_edges=false
    )
    compute_first_derivatives!(
      meg, ∂y_∂η, y, jaxis, ∂domain, backend; use_one_sided_on_edges=false
    )
    compute_first_derivatives!(
      meg, ∂z_∂ζ, z, kaxis, ∂domain, backend; use_one_sided_on_edges=false
    )

    compute_second_derivatives!(
      meg, ∂²x_∂ξ², ∂x_∂ξ, x, iaxis, ∂²domain, backend; use_one_sided_on_edges=false
    )
    compute_second_derivatives!(
      meg, ∂²y_∂η², ∂y_∂η, y, jaxis, ∂²domain, backend; use_one_sided_on_edges=false
    )
    compute_second_derivatives!(
      meg, ∂²z_∂ζ², ∂z_∂ζ, z, kaxis, ∂²domain, backend; use_one_sided_on_edges=false
    )

    @test all(∂x_∂ξ[∂domain] .≈ 1.0)
    @test all(∂y_∂η[∂domain] .≈ 2.0)
    @test all(∂z_∂ζ[∂domain] .≈ 3.0)
    @test all(∂²x_∂ξ²[∂²domain] .≈ 0.0)
    @test all(∂²y_∂η²[∂²domain] .≈ 0.0)
    @test all(∂²z_∂ζ²[∂²domain] .≈ 0.0)

    # #   Now use one-sided edge derivatives

    if order == 6
      x, y, z, ∂x_∂ξ, ∂y_∂η, ∂z_∂ζ, ∂²x_∂ξ², ∂²y_∂η², ∂²z_∂ζ² = init_xyz(celldims)
      compute_first_derivatives!(
        meg, ∂x_∂ξ, x, iaxis, ∂domain, backend; use_one_sided_on_edges=true
      )
      compute_first_derivatives!(
        meg, ∂y_∂η, y, jaxis, ∂domain, backend; use_one_sided_on_edges=true
      )
      compute_first_derivatives!(
        meg, ∂z_∂ζ, z, kaxis, ∂domain, backend; use_one_sided_on_edges=true
      )

      floor = 1e-14
      compute_second_derivatives!(
        meg,
        ∂²x_∂ξ²,
        ∂x_∂ξ,
        x,
        iaxis,
        ∂²domain,
        backend;
        use_one_sided_on_edges=true,
        floor=floor,
      )
      compute_second_derivatives!(
        meg,
        ∂²y_∂η²,
        ∂y_∂η,
        y,
        jaxis,
        ∂²domain,
        backend;
        use_one_sided_on_edges=true,
        floor=floor,
      )
      compute_second_derivatives!(
        meg,
        ∂²z_∂ζ²,
        ∂z_∂ζ,
        z,
        kaxis,
        ∂²domain,
        backend;
        use_one_sided_on_edges=true,
        floor=floor,
      )

      @test all(∂x_∂ξ[inner_domain] .≈ 1.0)
      @test all(∂y_∂η[inner_domain] .≈ 2.0)
      @test all(∂z_∂ζ[inner_domain] .≈ 3.0)

      @test all(∂²x_∂ξ²[inner_domain] .≈ 0.0)
      @test all(∂²y_∂η²[inner_domain] .≈ 0.0)
      @test all(∂²z_∂ζ²[inner_domain] .≈ 0.0)
    end
  end
end

@testset "MEG6 metrics" begin
  celldims = (30, 40)

  function init_xyz(cell_dims)
    x = zeros(cell_dims...)
    y = zeros(cell_dims...)

    for I in CartesianIndices(x)
      i, j = I.I
      x[I] = i
      y[I] = 2j
      #   z[i, j, k] = 3k
    end

    return (x, y)
  end

  x, y = init_xyz(celldims)
  iaxis, jaxis, kaxis = (1, 2, 3)
  full_domain = CartesianIndices(x)

  backend = get_backend(x)
  cell_center_metrics, edge_metrics = CurvilinearGrids.GridTypes.get_metric_soa(
    celldims, backend, Float64
  )

  for order in (6,)
    meg = MonotoneExplicitGradientScheme(order; celldims=celldims, use_cache=true)
    inner_domain = expand(full_domain, -meg.nhalo)

    CurvilinearGrids.MetricDiscretizationSchemes.forward_metrics!(
      meg, cell_center_metrics, x, y, inner_domain
    )
    CurvilinearGrids.MetricDiscretizationSchemes.conservative_metrics!(
      meg, cell_center_metrics, x, y, inner_domain
    )

    @test all(cell_center_metrics.x₁.ξ[inner_domain] .≈ 1.0)
    @test all(cell_center_metrics.x₂.η[inner_domain] .≈ 2.0)

    @test all(cell_center_metrics.ξ.x₁[inner_domain] .≈ 1.0)
    @test all(cell_center_metrics.η.x₂[inner_domain] .≈ 0.5)
  end
end

# begin
#   ni, nj = (40, 80)
#   x0, x1 = (0, 2)
#   y0, y1 = (1, 3)

#   T = Float64

#   x = collect(range(x0, x1; length=ni + 1))
#   y = collect(range(y0, y1; length=nj + 1))
#   x2d = zeros(T, length(x), length(y))
#   y2d = zeros(T, length(x), length(y))

#   @inbounds for j in eachindex(y)
#     for i in eachindex(x)
#       x2d[i, j] = x[i]
#       y2d[i, j] = y[j]
#     end
#   end

#   scheme_name = :MEG6
#   use_symmetric_conservative_metric_scheme = false
#   if scheme_name === :MEG6
#     order = 6
#     MetricDiscretizationScheme = MonotoneExplicitGradientScheme
#   elseif scheme_name === :MEG6_SYMMETRIC
#     order = 6
#     MetricDiscretizationScheme = MonotoneExplicitGradientScheme
#     use_symmetric_conservative_metric_scheme = true
#   else
#     error("Only $(:MEG6) is supported for now")
#   end

#   nhalo = DiscretizationSchemes.nhalo_lookup[scheme_name]

#   limits, iterators = get_iterators(size(x2d), halo_coords_included, nhalo)
#   celldims = size(iterators.cell.full)

#   backend = CPU()
#   meg = MetricDiscretizationScheme(
#     order;
#     use_cache=true,
#     celldims=celldims,
#     backend=backend,
#     T=T,
#     use_symmetric_conservative_metric_scheme=use_symmetric_conservative_metric_scheme,
#   )

#   node_coordinates, centroid_coordinates, node_velocities, nnodes = CurvilinearGrids.GridTypes._grid_constructor(
#     x2d, y2d, iterators; backend=backend
#   )

#   cell_center_metrics, edge_metrics = CurvilinearGrids.GridTypes.get_metric_soa(
#     celldims, backend, T
#   )
#   # cell_center_metrics, edge_metrics = CurvilinearGrids.GridTypes.get_metric_soa_rectilinear2d(
#   #   celldims, backend, T
#   # )

#   inner_domain = iterators.cell.domain

#   CurvilinearGrids.MetricDiscretizationSchemes.forward_metrics!(
#     meg, cell_center_metrics, StructArrays.components(centroid_coordinates)..., inner_domain
#   )

#   CurvilinearGrids.MetricDiscretizationSchemes.conservative_metrics!(
#     meg, cell_center_metrics, StructArrays.components(centroid_coordinates)..., inner_domain
#   )
#   cell_center_metrics.x₁.ξ[inner_domain]
# end

@testset "Stencil Floating Point Error" begin
  function central_derivative_naive(fm3, fm2, fm1, fp1, fp2, fp3)
    return (
      -1 / 60 * fm3 + 3 / 20 * fm2 - 3 / 4 * fm1 + 3 / 4 * fp1 - 3 / 20 * fp2 + 1 / 60 * fp3
    )
  end

  vals = 1e4ones(8)
  for i in eachindex(vals)
    pert = 100eps() * rand() # add perturbations > than epsilon (~2e-16)
    # @show pert
    vals[i] += pert
  end

  # The naive version reports a derivative of ~ 1e-13, when it _should_ be 0.0! The `central_derivative`
  # function uses the Kahan compensated sum to counteract floating-point roundoff error
  # @show abs(central_derivative_naive(vals[1:6]...))
  # @show central_derivative(vals[1:2]...)
  # @show central_derivative(vals[1:4]...)
  # @show central_derivative(vals[1:6]...)
  # @show central_derivative(vals[1:8]...)

  @test iszero(central_derivative(vals[1:2]...))
  @test iszero(central_derivative(vals[1:4]...))
  @test iszero(central_derivative(vals[1:6]...))
  @test iszero(central_derivative(vals[1:8]...))

  @test abs(central_derivative_naive(vals[1:6]...)) > 0 &&
    iszero(central_derivative(vals[1:6]...))
end