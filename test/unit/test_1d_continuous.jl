
@testset "Uniform ContinuousCurvilinearGrid1D" begin
  function get_uniform_mapping()
    function x(t, i, p)
      @unpack xmin, Δx = p
      return xmin + (i - 1) * Δx
    end

    return x
  end

  function get_uniform_params()
    celldims = (40,)
    xmin, xmax = (0.0, 2.0)

    ni, = celldims

    Δx = (xmax - xmin) / ni

    params = (; xmin, Δx)
    return params, celldims
  end

  params, celldims = get_uniform_params()
  x = get_uniform_mapping()

  backend = AutoForwardDiff()
  mesh = ContinuousCurvilinearGrid1D(x, params, celldims, :meg6, CPU(), backend)

  I1 = CurvilinearGrids.GridTypes.gcl(mesh.edge_metrics, mesh.iterators.cell.domain)
  @test all(abs.(extrema(I1)) .< eps())

  cell_volume = 0.05

  @test all(mesh.cell_center_metrics.J .≈ cell_volume)
  @test all(mesh.cell_center_metrics.x₁.ξ .≈ 0.05)
  @test all(mesh.cell_center_metrics.ξ̂.x₁ .≈ 1.0)
  @test all(mesh.cell_center_metrics.ξ.x₁ .≈ 20.0)
  @test all(mesh.cell_center_metrics.ξ.t .≈ 0.0)

  iaxis, jaxis, kaxis = (1, 2, 3)
  domain = mesh.iterators.cell.domain
  i₊½_domain = expand(domain, iaxis, -1)

  @test all(mesh.edge_metrics.i₊½.ξ̂.x₁[i₊½_domain] .≈ 1.0)
  @test all(mesh.edge_metrics.i₊½.ξ̂.t[i₊½_domain] .≈ 0.0)
end
