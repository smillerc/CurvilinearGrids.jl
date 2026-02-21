using CurvilinearGrids,
  Test,
  KernelAbstractions,
  DifferentiationInterface,
  BenchmarkTools,
  UnPack,
  CartesianDomains

@testset "Wavy MappedGrid3D" begin
  function wavy_mapping()
    function x(t, i, j, k, p)
      @unpack Ax, Ay, Az, n, Δx0, Δy0, Δz0, xmin = p
      xmin +
      Δx0 * ((i - 1) + Ax * sin(pi * n * (j - 1) * Δy0) * sin(pi * n * (k - 1) * Δz0))
    end

    function y(t, i, j, k, p)
      @unpack Ax, Ay, Az, n, Δx0, Δy0, Δz0, ymin = p
      ymin +
      Δy0 * ((j - 1) + Ay * sin(pi * n * (k - 1) * Δz0) * sin(pi * n * (i - 1) * Δx0))
    end

    function z(t, i, j, k, p)
      @unpack Ax, Ay, Az, n, Δx0, Δy0, Δz0, zmin = p
      zmin +
      Δz0 * ((k - 1) + Az * sin(pi * n * (i - 1) * Δx0) * sin(pi * n * (j - 1) * Δy0))
    end

    return (x, y, z)
  end

  function wavy_params()
    celldims = (41, 41, 41)
    ni, nj, nk = celldims
    Lx = Ly = Lz = 12

    xmin = -Lx / 2
    ymin = -Ly / 2
    zmin = -Lz / 2

    Δx0 = Lx / ni
    Δy0 = Ly / nj
    Δz0 = Lz / nk

    Ax = 0.2 / Δx0
    Ay = 0.4 / Δy0
    Az = 0.6 / Δz0

    n = 0.5

    params = (; Ax, Ay, Az, n, Δx0, Δy0, Δz0, xmin, ymin, zmin)
    return params, celldims
  end

  x, y, z = wavy_mapping()
  params, celldims = wavy_params()
  mesh = MappedGrid(x, y, z, params, celldims, :meg6; backend=CPU())
  I1, I2, I3 = CurvilinearGrids.GridTypes.gcl(
    face_metrics(mesh), mesh.iterators.cell.domain
  )
  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)
  @test all(abs.(extrema(I3)) .< 1e-14)
end

@testset "Wavy MappedGrid3D Edge Scheme Runtime Comparison" begin
  function wavy_mapping()
    function x(t, i, j, k, p)
      @unpack Ax, Ay, Az, n, Δx0, Δy0, Δz0, xmin = p
      xmin +
      Δx0 * ((i - 1) + Ax * sin(pi * n * (j - 1) * Δy0) * sin(pi * n * (k - 1) * Δz0))
    end

    function y(t, i, j, k, p)
      @unpack Ax, Ay, Az, n, Δx0, Δy0, Δz0, ymin = p
      ymin +
      Δy0 * ((j - 1) + Ay * sin(pi * n * (k - 1) * Δz0) * sin(pi * n * (i - 1) * Δx0))
    end

    function z(t, i, j, k, p)
      @unpack Ax, Ay, Az, n, Δx0, Δy0, Δz0, zmin = p
      zmin +
      Δz0 * ((k - 1) + Az * sin(pi * n * (i - 1) * Δx0) * sin(pi * n * (j - 1) * Δy0))
    end

    return (x, y, z)
  end

  function wavy_params()
    celldims = (41, 41, 41)
    ni, nj, nk = celldims
    Lx = Ly = Lz = 12

    xmin = -Lx / 2
    ymin = -Ly / 2
    zmin = -Lz / 2

    Δx0 = Lx / ni
    Δy0 = Ly / nj
    Δz0 = Lz / nk

    Ax = 0.2 / Δx0
    Ay = 0.4 / Δy0
    Az = 0.6 / Δz0
    n = 0.5

    params = (; Ax, Ay, Az, n, Δx0, Δy0, Δz0, xmin, ymin, zmin)
    return params, celldims
  end

  function run_wavy_case(x, y, z, params, celldims, scheme)
    mesh = MappedGrid(
      x, y, z, params, celldims, :meg6; backend=CPU(), conserved_metric_scheme=scheme
    )
    I1, I2, I3 = CurvilinearGrids.GridTypes.gcl(
      face_metrics(mesh), mesh.iterators.cell.domain
    )
    m1 = maximum(abs, I1)
    m2 = maximum(abs, I2)
    m3 = maximum(abs, I3)
    return (m1=m1, m2=m2, m3=m3, m=max(m1, m2, m3))
  end

  x, y, z = wavy_mapping()
  params, celldims = wavy_params()
  schemes = (
    ("EdgeInterpolationOrder1", CurvilinearGrids.GridTypes.EdgeInterpolationOrder1()),
    ("EdgeInterpolationOrder2", CurvilinearGrids.GridTypes.EdgeInterpolationOrder2()),
    ("EdgeInterpolationOrder3", CurvilinearGrids.GridTypes.EdgeInterpolationOrder3()),
  )

  # Warm up each scheme to avoid JIT compilation noise in runtime comparison.
  for (_, scheme) in schemes
    run_wavy_case(x, y, z, params, celldims, scheme)
  end

  runtimes = Float64[]
  for (name, scheme) in schemes
    t0 = time_ns()
    errors = run_wavy_case(x, y, z, params, celldims, scheme)
    elapsed_seconds = (time_ns() - t0) * 1e-9
    push!(runtimes, elapsed_seconds)

    @info "Wavy MappedGrid3D edge interpolation comparison" scheme = name max_I1 =
      errors.m1 max_I2 = errors.m2 max_I3 = errors.m3 max_error = errors.m runtime_seconds =
      elapsed_seconds

    @test errors.m1 < 1e-14
    @test errors.m2 < 1e-14
    @test errors.m3 < 1e-14
    @test elapsed_seconds > 0
  end

  @test runtimes[1] < runtimes[2]
  @test runtimes[2] < runtimes[3]
end

@testset "Sphere Sector MappedGrid3D" begin
  function get_sector_parameters()
    celldims = (40, 40, 40)
    ni, nj, nk = celldims
    rmin, rmax = 1.0, 4.0
    θmin, θmax = π / 2 - deg2rad(5), π / 2 + deg2rad(5)
    ϕmin, ϕmax = -deg2rad(10), deg2rad(10)

    Δr = (rmax - rmin) / ni
    Δθ = (θmax - θmin) / nj
    Δϕ = (ϕmax - ϕmin) / nk

    params = (; Δr, rmax, rmin, Δθ, θmax, θmin, Δϕ, ϕmax, ϕmin)
    return params, celldims
  end

  function get_sector_mapping()
    function x(t, i, j, k, p)
      @unpack Δr, rmax, rmin, Δθ, θmax, θmin, Δϕ, ϕmax, ϕmin = p
      r(i) = (rmin + (i - 1) * Δr)
      θ(j) = (θmin + (j - 1) * Δθ)
      ϕ(k) = (ϕmin + (k - 1) * Δϕ)
      return r(i) * sin(θ(j)) * cos(ϕ(k))
    end

    function y(t, i, j, k, p)
      @unpack Δr, rmax, rmin, Δθ, θmax, θmin, Δϕ, ϕmax, ϕmin = p
      r(i) = (rmin + (i - 1) * Δr)
      θ(j) = (θmin + (j - 1) * Δθ)
      ϕ(k) = (ϕmin + (k - 1) * Δϕ)
      return r(i) * sin(θ(j)) * sin(ϕ(k))
    end

    function z(t, i, j, k, p)
      @unpack Δr, rmax, rmin, Δθ, θmax, θmin, Δϕ, ϕmax, ϕmin = p
      r(i) = (rmin + (i - 1) * Δr)
      θ(j) = (θmin + (j - 1) * Δθ)
      ϕ(k) = (ϕmin + (k - 1) * Δϕ)
      return r(i) * cos(θ(j))
    end

    return (x, y, z)
  end

  sector_params, celldims = get_sector_parameters()
  x, y, z = get_sector_mapping()

  backend = AutoForwardDiff()
  mesh = MappedGrid(
    x, y, z, sector_params, celldims, :meg6; backend=CPU(), diff_backend=backend
  )
  I1, I2, I3 = CurvilinearGrids.GridTypes.gcl(
    face_metrics(mesh), mesh.iterators.cell.domain
  )

  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)
  @test all(abs.(extrema(I3)) .< 1e-14)
end

@testset "Uniform MappedGrid3D" begin
  function get_uniform_mapping()
    function x(t, i, j, k, p)
      @unpack xmin, Δx = p
      return xmin + (i - 1) * Δx
    end

    function y(t, i, j, k, p)
      @unpack ymin, Δy = p
      return ymin + (j - 1) * Δy
    end

    function z(t, i, j, k, p)
      @unpack zmin, Δz = p
      return zmin + (k - 1) * Δz
    end

    return x, y, z
  end

  function get_uniform_params()
    celldims = (40, 80, 120)
    xmin, xmax = (0.0, 2.0)
    ymin, ymax = (1, 3)
    zmin, zmax = (-1, 2)

    ni, nj, nk = celldims

    Δx = (xmax - xmin) / ni
    Δy = (ymax - ymin) / nj
    Δz = (zmax - zmin) / nk

    params = (; xmin, ymin, zmin, Δx, Δy, Δz)
    return params, celldims
  end

  params, celldims = get_uniform_params()
  x, y, z = get_uniform_mapping()

  backend = AutoForwardDiff()
  mesh = MappedGrid(x, y, z, params, celldims, :meg6; backend=CPU(), diff_backend=backend)

  I1, I2, I3 = CurvilinearGrids.GridTypes.gcl(
    face_metrics(mesh), mesh.iterators.cell.domain
  )
  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)
  @test all(abs.(extrema(I3)) .< 1e-14)

  cm = cell_metrics(mesh)
  cell_volume = 0.05 * 0.025 * 0.025

  @test all([m.J for m in cm.forward] .≈ cell_volume)

  @test all([m[1, 1] for m in cm.forward] .≈ 0.05)
  @test all([m[2, 2] for m in cm.forward] .≈ 0.025)
  @test all([m[3, 3] for m in cm.forward] .≈ 0.025)
  @test all([m[1, 2] for m in cm.forward] .≈ 0.0)
  @test all([m[1, 3] for m in cm.forward] .≈ 0.0)
  @test all([m[2, 1] for m in cm.forward] .≈ 0.0)
  @test all([m[2, 3] for m in cm.forward] .≈ 0.0)
  @test all([m[3, 1] for m in cm.forward] .≈ 0.0)
  @test all([m[3, 2] for m in cm.forward] .≈ 0.0)

  @test all([m[1, 1] for m in cm.inverse] .≈ 20.0)
  @test all([m[2, 2] for m in cm.inverse] .≈ 40.0)
  @test all([m[3, 3] for m in cm.inverse] .≈ 40.0)
  @test all([m[1, 2] for m in cm.inverse] .≈ 0.0)
  @test all([m[1, 3] for m in cm.inverse] .≈ 0.0)
  @test all([m[2, 1] for m in cm.inverse] .≈ 0.0)
  @test all([m[2, 3] for m in cm.inverse] .≈ 0.0)
  @test all([m[3, 1] for m in cm.inverse] .≈ 0.0)
  @test all([m[3, 2] for m in cm.inverse] .≈ 0.0)

  @test all([(m[1, 1] * m.J) for m in cm.inverse] .≈ 0.000625)
  @test all([(m[2, 2] * m.J) for m in cm.inverse] .≈ 0.00125)
  @test all([(m[3, 3] * m.J) for m in cm.inverse] .≈ 0.00125)

  iaxis, jaxis, kaxis = (1, 2, 3)
  domain = mesh.iterators.cell.domain
  i₊½_domain = expand(domain, iaxis, -1)
  j₊½_domain = expand(domain, jaxis, -1)
  k₊½_domain = expand(domain, kaxis, -1)

  fm = face_metrics(mesh)
  for (axis, dom) in ((1, i₊½_domain), (2, j₊½_domain), (3, k₊½_domain))
    cons = fm[axis].conserved[dom]
    @test all([m[1, 1] for m in cons] .≈ 0.000625)
    @test all([m[1, 2] for m in cons] .≈ 0.0)
    @test all([m[1, 3] for m in cons] .≈ 0.0)
    @test all([m[2, 1] for m in cons] .≈ 0.0)
    @test all([m[2, 2] for m in cons] .≈ 0.00125)
    @test all([m[2, 3] for m in cons] .≈ 0.0)
    @test all([m[3, 1] for m in cons] .≈ 0.0)
    @test all([m[3, 2] for m in cons] .≈ 0.0)
    @test all([m[3, 3] for m in cons] .≈ 0.00125)
  end
end
