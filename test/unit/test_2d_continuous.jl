
function uniform_mapping(xmin, xmax, ymin, ymax, ncells::NTuple{2,Int})
  ni, nj = ncells

  Δx = (xmax - xmin) / ni
  Δy = (ymax - ymin) / nj

  params = (; Δx, Δy, xmin, xmax, ymin, ymax)

  function x(t, i, j, p)
    @unpack Δx = p
    return xmin + (i - 1) * Δx
  end

  function y(t, i, j, p)
    @unpack Δy = p
    return ymin + (j - 1) * Δy
  end

  return (x, y, params)
end

function cylindrical_sector_mapping(rmin, rmax, θmin, θmax, ncells::NTuple{2,Int})
  ni, nj = ncells

  # Uniform spacings
  Δr = (rmax - rmin) / ni
  Δθ = (θmax - θmin) / nj

  params = (; Δr, Δθ, rmax, rmin, θmax, θmin)

  function x(t, i, j, p)
    @unpack rmin, θmin, Δr, Δθ = p
    r(i) = (rmin + (i - 1) * Δr)
    θ(j) = (θmin + (j - 1) * Δθ)

    return r(i) * cos(θ(j))
  end

  function y(t, i, j, p)
    @unpack rmin, θmin, Δr, Δθ = p

    r(i) = (rmin + (i - 1) * Δr)
    θ(j) = (θmin + (j - 1) * Δθ)

    return r(i) * sin(θ(j))
  end

  return (x, y, params)
end

function wavy_mapping(ncells::NTuple{2,Int})
  ni, nj = ncells
  Lx = Ly = 12

  xmin = -Lx / 2
  ymin = -Ly / 2

  Δx0 = Lx / ni
  Δy0 = Ly / nj

  Ax = 0.25 / Δx0
  Ay = 0.5 / Δy0

  n = 0.5

  params = (; Lx, Ly, xmin, ymin, Δx0, Δy0, Ax, Ay, n)

  function x(t, i, j, p)
    @unpack xmin, Ax, Δy0, Δx0 = p

    return xmin + Δx0 * ((i - 1) + Ax * sinpi(n * (j - 1) * Δy0))
  end

  function y(t, i, j, p)
    @unpack ymin, Ay, Δy0, Δx0 = p

    return ymin + Δy0 * ((j - 1) + Ay * sinpi(n * (i - 1) * Δx0))
  end

  return (x, y, params)
end

@testset "Wavy ContinuousCurvilinearGrid2D" begin
  ncells = (401, 401)

  (x, y, params) = wavy_mapping(ncells)

  @unpack Lx, Ly, xmin, ymin, Δx0, Δy0, Ax, Ay, n = params
  params_2 = (; Lx, Ly, xmin, ymin, Δx0, Δy0, Ax=2Ax, Ay=0.5Ay, n)

  mesh = ContinuousCurvilinearGrid2D(x, y, params, ncells, :meg6, CPU())

  t, i, j = (0, 10, 20)
  c1 = coord(mesh, t, i, j)
  # @show c1
  # save_vtk(mesh, "wavy_2d_ad_t1")
  I1, I2 = CurvilinearGrids.GridTypes.gcl(mesh.edge_metrics, mesh.iterators.cell.domain)

  # @show extrema(I1)
  # @show extrema(I2)
  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)

  CurvilinearGrids.GridTypes.update_mapping_functions!(mesh, t, params_2)

  c2 = coord(mesh, t, i, j)

  @test !isapprox(c1, c2)
  # @show c2
  # save_vtk(mesh, "wavy_2d_ad_t2")

  I1, I2 = CurvilinearGrids.GridTypes.gcl(mesh.edge_metrics, mesh.iterators.cell.domain)

  # @show extrema(I1)
  # @show extrema(I2)
  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)

  nothing
end

@testset "Cylindrical Sector ContinuousCurvilinearGrid2D" begin
  nhalo = 5
  nr, nθ = 41, 41
  celldims = (nr, nθ)
  rmin, rmax = 1.0, 4.0
  θmin, θmax = deg2rad(5), deg2rad(25)   # narrow polar band around equator

  (x, y, params) = cylindrical_sector_mapping(rmin, rmax, θmin, θmax, celldims)

  mesh = ContinuousCurvilinearGrid2D(x, y, params, celldims, :meg6, CPU())
  I1, I2 = CurvilinearGrids.GridTypes.gcl(mesh.edge_metrics, mesh.iterators.cell.domain)
  # @show extrema(I1)
  # @show extrema(I2)
  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)
end

@testset "Uniform ContinuousCurvilinearGrid2D" begin
  x0, x1 = (0.0, 2.0)
  y0, y1 = (1, 3)
  celldims = (40, 80)
  (x, y, params) = uniform_mapping(x0, x1, y0, y1, celldims)

  mesh = ContinuousCurvilinearGrid2D(x, y, params, celldims, :meg6, CPU())
  I1, I2 = CurvilinearGrids.GridTypes.gcl(mesh.edge_metrics, mesh.iterators.cell.domain)
  @test all(abs.(extrema(I1)) .< eps())
  @test all(abs.(extrema(I2)) .< eps())

  cell_volume = 0.05 * 0.025

  @test all(mesh.cell_center_metrics.J .≈ cell_volume)

  @test all(mesh.cell_center_metrics.x₁.ξ .≈ 0.05)
  @test all(mesh.cell_center_metrics.x₂.ξ .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₁.η .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₂.η .≈ 0.025)

  @test all(mesh.cell_center_metrics.ξ̂.x₁ .≈ 0.025)
  @test all(mesh.cell_center_metrics.ξ̂.x₂ .≈ 0.0)
  @test all(mesh.cell_center_metrics.η̂.x₁ .≈ 0.0)
  @test all(mesh.cell_center_metrics.η̂.x₂ .≈ 0.05)

  @test all(mesh.cell_center_metrics.ξ.x₁ .≈ 20.0)
  @test all(mesh.cell_center_metrics.ξ.x₂ .≈ 0.0)
  @test all(mesh.cell_center_metrics.η.x₁ .≈ 0.0)
  @test all(mesh.cell_center_metrics.η.x₂ .≈ 40.0)

  @test all(mesh.cell_center_metrics.ξ.t .≈ 0.0)
  @test all(mesh.cell_center_metrics.η.t .≈ 0.0)

  iaxis, jaxis, kaxis = (1, 2, 3)
  domain = mesh.iterators.cell.domain
  i₊½_domain = expand(domain, iaxis, -1)
  j₊½_domain = expand(domain, jaxis, -1)

  @test all(mesh.edge_metrics.i₊½.ξ̂.x₁[i₊½_domain] .≈ 0.025)
  @test all(mesh.edge_metrics.j₊½.ξ̂.x₁[j₊½_domain] .≈ 0.025)
  @test all(mesh.edge_metrics.i₊½.ξ̂.x₂[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ξ̂.x₂[j₊½_domain] .≈ 0.0)

  @test all(mesh.edge_metrics.i₊½.η̂.x₁[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.η̂.x₁[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.η̂.x₂[i₊½_domain] .≈ 0.05)
  @test all(mesh.edge_metrics.j₊½.η̂.x₂[j₊½_domain] .≈ 0.05)

  @test all(mesh.edge_metrics.i₊½.ξ̂.t[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ξ̂.t[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.η̂.t[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.η̂.t[j₊½_domain] .≈ 0.0)
end

@testset "ContinuousCurvilinearGrid2D vs CurvilinearGrid2D" begin
  nhalo = 5
  nr, nθ = 50, 51
  celldims = (nr, nθ)

  rmin, rmax = (0.0400, 0.05500)

  θmin, θmax = ((5, 25) .|> deg2rad)

  (x, y, params) = cylindrical_sector_mapping(rmin, rmax, θmin, θmax, celldims)

  @info "ContinuousCurvilinearGrid2D"
  cm = ContinuousCurvilinearGrid2D(x, y, params, celldims, :meg6, CPU())

  I1, I2 = CurvilinearGrids.GridTypes.gcl(cm.edge_metrics, cm.iterators.cell.domain)
  # @show extrema(I1)
  # @show extrema(I2)
  @test all(abs.(extrema(I1)) .< eps())
  @test all(abs.(extrema(I2)) .< eps())
  # CurvilinearGrids.save_vtk(cm, "cylindrical_sector_ad")

  @info "CurvilinearGrid2D"
  xdom = cm.node_coordinates.x[cm.iterators.node.full]
  ydom = cm.node_coordinates.y[cm.iterators.node.full]
  dm = CurvilinearGrid2D(xdom, ydom, :meg6; halo_coords_included=true)
  I1, I2 = CurvilinearGrids.GridTypes.gcl(dm.edge_metrics, dm.iterators.cell.domain)
  # @show extrema(I1)
  # @show extrema(I2)
  @test all(abs.(extrema(I1)) .< eps())
  @test all(abs.(extrema(I2)) .< eps())
  # CurvilinearGrids.save_vtk(dm, "cylindrical_sector_meg")

  dom = cm.iterators.cell.domain

  # for dim in (:ξ, :η, :ξ̂, :η̂, :x₁, :x₂)
  #   for ((dm_name, dm_component), (cm_name, cm_component)) in zip(
  #     StructArrays.components(dm.cell_center_metrics[dim]) |> pairs,
  #     StructArrays.components(cm.cell_center_metrics[dim]) |> pairs,
  #   )

  #     # we have to use a relatively coarse tolerance for this mesh! This is why AD is better :)
  #     @test all(isapprox.(dm_component[dom], cm_component[dom], rtol=1e-5))
  #     #   passes = all(isapprox.(dm_component[dom], cm_component[dom], rtol=1e-5))
  #     #   @info "Dim: $dim, $dm_name, passes? $passes"
  #     #   if !passes
  #     #     @show extrema(dm_component[dom])
  #     #     @show extrema(cm_component[dom])
  #     #     println()
  #     #   end
  #   end
  #   # println()
  # end

  nothing
end
