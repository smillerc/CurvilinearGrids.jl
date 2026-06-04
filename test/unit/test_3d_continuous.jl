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
  mesh = MappedGrid(x, y, z, params, celldims, 5; backend=CPU())
  I1, I2, I3 = CurvilinearGrids.GridTypes.gcl(
    face_metrics(mesh), mesh.iterators.cell.domain
  )
  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)
  @test all(abs.(extrema(I3)) .< 1e-14)
end

@testset "Wavy MappedGrid3D Face Reconstruction Runtime Comparison" begin
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
      x, y, z, params, celldims, 5; backend=CPU(), conserved_metric_scheme=scheme
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
    (
      "EndpointAverageReconstruction",
      CurvilinearGrids.GridTypes.EndpointAverageReconstruction(),
    ),
    (
      "GradientCorrectedReconstruction",
      CurvilinearGrids.GridTypes.GradientCorrectedReconstruction(),
    ),
    (
      "CurvatureCorrectedReconstruction",
      CurvilinearGrids.GridTypes.CurvatureCorrectedReconstruction(),
    ),
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

    @info "Wavy MappedGrid3D face reconstruction comparison" scheme = name max_I1 = errors.m1 max_I2 =
      errors.m2 max_I3 = errors.m3 max_error = errors.m runtime_seconds = elapsed_seconds

    @test errors.m1 < 1e-14
    @test errors.m2 < 1e-14
    @test errors.m3 < 1e-14
    @test elapsed_seconds > 0
  end

  @test runtimes[1] < runtimes[2]
  @test runtimes[2] < runtimes[3]
end

@testset "AD Thomas-Lombard MappedGrid3D" begin
  function max_gcl_residual(mesh)
    I1, I2, I3 = CurvilinearGrids.GridTypes.gcl(
      face_metrics(mesh), mesh.iterators.cell.domain
    )
    domain = mesh.iterators.cell.domain
    return maximum((
      maximum(abs, I1[domain]), maximum(abs, I2[domain]), maximum(abs, I3[domain])
    ))
  end

  function max_active_row_difference(reference, prototype)
    reference_face_metrics = face_metrics(reference)
    prototype_face_metrics = face_metrics(prototype)
    domain = reference.iterators.cell.domain
    max_error = 0.0
    for axis in 1:3, I in domain, component in 1:3
      reference_value =
        reference_face_metrics[axis].conserved[I].jacobian_matrix[axis, component]
      prototype_value =
        prototype_face_metrics[axis].conserved[I].jacobian_matrix[axis, component]
      max_error = max(max_error, abs(reference_value - prototype_value))
    end
    return max_error
  end

  function warped_mapping()
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

    return x, y, z
  end

  function warped_params(celldims)
    ni, nj, nk = celldims
    Lx = Ly = Lz = 12
    Δx0 = Lx / ni
    Δy0 = Ly / nj
    Δz0 = Lz / nk
    return (;
      Ax=0.2 / Δx0,
      Ay=0.4 / Δy0,
      Az=0.6 / Δz0,
      n=0.5,
      Δx0,
      Δy0,
      Δz0,
      xmin=-Lx / 2,
      ymin=-Ly / 2,
      zmin=-Lz / 2,
    )
  end

  function spherical_wedge_mapping()
    function x(t, i, j, k, p)
      @unpack rmin, θmin, ϕmin, Δr, Δθ, Δϕ = p
      r = rmin + (i - 1) * Δr
      θ = θmin + (j - 1) * Δθ
      ϕ = ϕmin + (k - 1) * Δϕ
      return r * sin(θ) * cos(ϕ)
    end

    function y(t, i, j, k, p)
      @unpack rmin, θmin, ϕmin, Δr, Δθ, Δϕ = p
      r = rmin + (i - 1) * Δr
      θ = θmin + (j - 1) * Δθ
      ϕ = ϕmin + (k - 1) * Δϕ
      return r * sin(θ) * sin(ϕ)
    end

    function z(t, i, j, k, p)
      @unpack rmin, θmin, Δr, Δθ = p
      r = rmin + (i - 1) * Δr
      θ = θmin + (j - 1) * Δθ
      return r * cos(θ)
    end

    return x, y, z
  end

  function spherical_wedge_params(celldims)
    ni, nj, nk = celldims
    rmin, rmax = 1.0, 4.0
    θmin, θmax = π / 2 - deg2rad(5), π / 2 + deg2rad(5)
    ϕmin, ϕmax = -deg2rad(10), deg2rad(10)
    return (;
      rmin,
      θmin,
      ϕmin,
      Δr=(rmax - rmin) / ni,
      Δθ=(θmax - θmin) / nj,
      Δϕ=(ϕmax - ϕmin) / nk,
    )
  end


  @testset "AD Thomas-Lombard mapped" begin
    ad_tl_cases = (
      ("warped", warped_mapping(), warped_params((4, 4, 4)), (4, 4, 4)),
      (
        "spherical wedge",
        spherical_wedge_mapping(),
        spherical_wedge_params((4, 4, 4)),
        (4, 4, 4),
      ),
    )

    for (name, mapping, params, celldims) in ad_tl_cases
      @testset "$name" begin
        x, y, z = mapping
        reference = MappedGrid(
          x,
          y,
          z,
          params,
          celldims,
          3;
          backend=CPU(),
          diff_backend=AutoForwardDiff(),
          conserved_metric_scheme=CurvatureCorrectedReconstruction(),
        )
        ad_tl = MappedGrid(
          x,
          y,
          z,
          params,
          celldims,
          3;
          backend=CPU(),
          diff_backend=AutoForwardDiff(),
          conserved_metric_scheme=ADThomasLombardMetric(),
        )

        @test max_gcl_residual(ad_tl) < 1e-12
        @test max_active_row_difference(reference, ad_tl) < 1e-12
      end
    end
  end

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
    x,
    y,
    z,
    sector_params,
    celldims,
    5;
    backend=CPU(),
    diff_backend=backend,
    conserved_metric_scheme=CurvilinearGrids.GridTypes.EndpointAverageReconstruction(),
  )
  I1, I2, I3 = CurvilinearGrids.GridTypes.gcl(
    face_metrics(mesh), mesh.iterators.cell.domain
  )

  @show extrema(I1)
  @show extrema(I2)
  @show extrema(I3)

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
  mesh = MappedGrid(x, y, z, params, celldims, 5; backend=CPU(), diff_backend=backend)

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
