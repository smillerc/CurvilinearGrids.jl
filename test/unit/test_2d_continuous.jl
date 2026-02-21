function wavy_params()
  celldims = (401, 401)
  ni, nj = celldims
  Lx = Ly = 12

  xmin = -Lx / 2
  ymin = -Ly / 2

  Δx0 = Lx / ni
  Δy0 = Ly / nj

  Ax = 0.25 / Δx0
  Ay = 0.5 / Δy0

  n = 0.5

  params = (; Lx, Ly, xmin, ymin, Δx0, Δy0, Ax, Ay, n)

  return params, celldims
end

function wavy_mapping()
  function x(t, i, j, p)
    @unpack xmin, Ax, Δy0, Δx0, n = p
    return xmin + Δx0 * ((i - 1) + Ax * sinpi(n * (j - 1) * Δy0))
  end

  function y(t, i, j, p)
    @unpack ymin, Ay, Δy0, Δx0, n = p
    return ymin + Δy0 * ((j - 1) + Ay * sinpi(n * (i - 1) * Δx0))
  end

  return x, y
end

function cylindrical_params()
  celldims = (41, 41)
  ni, nj = celldims

  rmin, rmax = 1.0, 4.0
  θmin, θmax = deg2rad(5), deg2rad(25)

  Δr = (rmax - rmin) / ni
  Δθ = (θmax - θmin) / nj

  params = (; Δr, Δθ, rmax, rmin, θmax, θmin)
  return params, celldims
end

function cylindrical_mapping()
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

  return x, y
end

@testset "Wavy MappedGrid2D" begin
  params, celldims = wavy_params()
  x, y = wavy_mapping()

  @unpack Lx, Ly, xmin, ymin, Δx0, Δy0, Ax, Ay, n = params
  params_2 = (; Lx, Ly, xmin, ymin, Δx0, Δy0, Ax=2Ax, Ay=0.5Ay, n)

  mesh = MappedGrid(x, y, params, celldims, 5; backend=CPU())

  i, j = (10, 20)
  c1 = coord(mesh, (i, j))

  I1, I2 = CurvilinearGrids.GridTypes.gcl(face_metrics(mesh), mesh.iterators.cell.domain)
  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)

  CurvilinearGrids.GridTypes.update!(mesh, 0.0, params_2)
  c2 = coord(mesh, (i, j))

  @test !isapprox(c1, c2)

  I1, I2 = CurvilinearGrids.GridTypes.gcl(face_metrics(mesh), mesh.iterators.cell.domain)
  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)
end

@testset "Cylindrical Sector MappedGrid2D" begin
  params, celldims = cylindrical_params()
  x, y = cylindrical_mapping()

  mesh = MappedGrid(x, y, params, celldims, 5; backend=CPU())
  I1, I2 = CurvilinearGrids.GridTypes.gcl(face_metrics(mesh), mesh.iterators.cell.domain)
  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)
end

@testset "Uniform MappedGrid2D" begin
  function get_uniform_mapping()
    function x(t, i, j, p)
      @unpack xmin, Δx = p
      return xmin + (i - 1) * Δx
    end

    function y(t, i, j, p)
      @unpack ymin, Δy = p
      return ymin + (j - 1) * Δy
    end

    return x, y
  end

  function get_uniform_params()
    celldims = (40, 80)
    xmin, xmax = (0.0, 2.0)
    ymin, ymax = (1, 3)

    ni, nj = celldims

    Δx = (xmax - xmin) / ni
    Δy = (ymax - ymin) / nj

    params = (; xmin, ymin, Δx, Δy)
    return params, celldims
  end

  params, celldims = get_uniform_params()
  x, y = get_uniform_mapping()

  backend = AutoForwardDiff()
  mesh = MappedGrid(x, y, params, celldims, 5; backend=CPU(), diff_backend=backend)

  I1, I2 = CurvilinearGrids.GridTypes.gcl(face_metrics(mesh), mesh.iterators.cell.domain)
  @test all(abs.(extrema(I1)) .< eps())
  @test all(abs.(extrema(I2)) .< eps())

  cm = cell_metrics(mesh)
  cell_volume = 0.05 * 0.025

  @test all([m.J for m in cm.forward] .≈ cell_volume)
  @test all([m[1, 1] for m in cm.forward] .≈ 0.05)
  @test all([m[1, 2] for m in cm.forward] .≈ 0.0)
  @test all([m[2, 1] for m in cm.forward] .≈ 0.0)
  @test all([m[2, 2] for m in cm.forward] .≈ 0.025)

  @test all([m[1, 1] for m in cm.inverse] .≈ 20.0)
  @test all([m[1, 2] for m in cm.inverse] .≈ 0.0)
  @test all([m[2, 1] for m in cm.inverse] .≈ 0.0)
  @test all([m[2, 2] for m in cm.inverse] .≈ 40.0)

  @test all([(m[1, 1] * m.J) for m in cm.inverse] .≈ 0.025)
  @test all([(m[1, 2] * m.J) for m in cm.inverse] .≈ 0.0)
  @test all([(m[2, 1] * m.J) for m in cm.inverse] .≈ 0.0)
  @test all([(m[2, 2] * m.J) for m in cm.inverse] .≈ 0.05)

  iaxis, jaxis = (1, 2)
  domain = mesh.iterators.cell.domain
  i₊½_domain = expand(domain, iaxis, -1)
  j₊½_domain = expand(domain, jaxis, -1)
  fm = face_metrics(mesh)

  @test all([m[1, 1] for m in fm[1].conserved[i₊½_domain]] .≈ 0.025)
  @test all([m[1, 1] for m in fm[2].conserved[j₊½_domain]] .≈ 0.025)
  @test all([m[1, 2] for m in fm[1].conserved[i₊½_domain]] .≈ 0.0)
  @test all([m[1, 2] for m in fm[2].conserved[j₊½_domain]] .≈ 0.0)
  @test all([m[2, 1] for m in fm[1].conserved[i₊½_domain]] .≈ 0.0)
  @test all([m[2, 1] for m in fm[2].conserved[j₊½_domain]] .≈ 0.0)
  @test all([m[2, 2] for m in fm[1].conserved[i₊½_domain]] .≈ 0.05)
  @test all([m[2, 2] for m in fm[2].conserved[j₊½_domain]] .≈ 0.05)
end

@testset "MappedGrid2D vs DiscreteGrid2D" begin
  params, celldims = cylindrical_params()
  x, y = cylindrical_mapping()

  mm = MappedGrid(x, y, params, celldims, 5; backend=CPU())

  I1, I2 = CurvilinearGrids.GridTypes.gcl(face_metrics(mm), mm.iterators.cell.domain)
  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)

  @info "DiscreteGrid2D"
  xdom = mm.node_coordinates[1][mm.iterators.node.full]
  ydom = mm.node_coordinates[2][mm.iterators.node.full]
  dm = DiscreteGrid(xdom, ydom, 5; halo_coords_included=true)

  I1, I2 = CurvilinearGrids.GridTypes.gcl(face_metrics(dm), dm.iterators.cell.domain)
  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)

  nothing
end
