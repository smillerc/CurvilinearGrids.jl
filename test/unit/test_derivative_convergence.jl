
using Test
using Plots
using LinearAlgebra
using StaticArrays
using KernelAbstractions
using CartesianDomains
using CurvilinearGrids
using Polynomials

@testset "SphericalGrid3D Derivative Convergence (r, θ, ϕ)" begin

  # -------------------------------------------------------------
  # Analytic fields
  # -------------------------------------------------------------

  # r-only
  T_r(r) = r^3
  dTdr_analytic(r) = 3r^2

  # theta-only
  T_θ(θ) = cos(3θ)
  dTdθ_analytic(θ) = -3sin(3θ)

  # phi-only
  T_ϕ(ϕ) = sin(4ϕ)
  dTdϕ_analytic(ϕ) = 4cos(4ϕ)

  # -------------------------------------------------------------
  # Helpers
  # -------------------------------------------------------------
  function convergence_rate(hs, errs)
    # Fit: log(err) = p log(h) + C
    p = Polynomials.fit(log.(hs), log.(errs), 1)[1]
    return p
  end

  # -------------------------------------------------------------
  # 1. Radial derivative convergence test
  # -------------------------------------------------------------
  function radial_derivative_error(n_r)
    backend = KernelAbstractions.CPU()
    nhalo = 2

    θmid = π / 2
    ϕmid = 0.0
    nr, ntheta, nphi = (n_r, n_r, n_r)

    r_nodes = range(1.0, 2.0; length=nr) |> collect
    θ_nodes = range(θmid - 0.1, θmid + 0.1; length=ntheta) |> collect
    ϕ_nodes = range(ϕmid - 0.1, ϕmid + 0.1; length=nphi) |> collect

    grid = SphericalGrid3D(
      r_nodes, θ_nodes, ϕ_nodes, nhalo, backend; halo_coords_included=true
    )

    rc = grid.centroid_coordinates.r

    T = similar(grid.cell_volumes)

    # fill T(r) = r^3
    for I in grid.iterators.cell.full
      i, j, k = I.I
      T[I] = T_r(rc[i])
    end

    maxerr = 0.0

    for I in grid.iterators.cell.domain
      i, j, k = I.I
      r = rc[i]
      dT_num = cell_center_derivative(grid, T, I, 1)
      dT_ex = dTdr_analytic(r)
      maxerr = max(maxerr, abs(dT_num - dT_ex))
    end

    return maxerr
  end

  # -------------------------------------------------------------
  # 2. Theta derivative convergence test
  # -------------------------------------------------------------
  function theta_derivative_error(nθ)
    backend = KernelAbstractions.CPU()
    nhalo = 2

    θmid = π / 2
    ϕmid = 0.0
    nr, ntheta, nphi = (nθ, nθ, nθ)

    r_nodes = range(1.0, 2.0; length=nr) |> collect
    θ_nodes = range(θmid - 0.1, θmid + 0.1; length=ntheta) |> collect
    ϕ_nodes = range(ϕmid - 0.1, ϕmid + 0.1; length=nphi) |> collect

    grid = SphericalGrid3D(
      r_nodes, θ_nodes, ϕ_nodes, nhalo, backend; halo_coords_included=true
    )

    rc = grid.centroid_coordinates.r
    θc = grid.centroid_coordinates.θ

    T = similar(grid.cell_volumes)

    # fill T(θ) = cos(3θ)
    for I in grid.iterators.cell.full
      i, j, k = I.I
      T[I] = T_θ(θc[j])
    end

    maxerr = 0.0

    for I in grid.iterators.cell.domain
      i, j, k = I.I
      θ = θc[j]

      dT_dθ_coord = cell_center_derivative(grid, T, I, 2)
      dT_exact = (1 / rc[i]) * dTdθ_analytic(θ)

      maxerr = max(maxerr, abs(dT_dθ_coord - dT_exact))
    end

    return maxerr
  end

  # -------------------------------------------------------------
  # 3. Phi derivative convergence test
  # -------------------------------------------------------------
  function phi_derivative_error(nϕ)
    backend = KernelAbstractions.CPU()
    nhalo = 2

    θmid = π / 2
    ϕmid = 0.0
    nr, ntheta, nphi = (nϕ, nϕ, nϕ)

    r_nodes = range(1.0, 2.0; length=nr) |> collect
    θ_nodes = range(θmid - 0.1, θmid + 0.1; length=ntheta) |> collect
    ϕ_nodes = range(ϕmid - 0.1, ϕmid + 0.1; length=nphi) |> collect

    grid = SphericalGrid3D(
      r_nodes, θ_nodes, ϕ_nodes, nhalo, backend; halo_coords_included=true
    )

    rc = grid.centroid_coordinates.r
    θc = grid.centroid_coordinates.θ
    ϕc = grid.centroid_coordinates.ϕ

    T = similar(grid.cell_volumes)

    for I in grid.iterators.cell.full
      i, j, k = I.I
      T[I] = T_ϕ(ϕc[k])
    end

    maxerr = 0.0

    for I in grid.iterators.cell.domain
      i, j, k = I.I
      ϕ = ϕc[k]

      dT_dϕ = cell_center_derivative(grid, T, I, 3)
      dT_exact = (1 / (rc[i] * sin(θc[j]))) * dTdϕ_analytic(ϕ)

      maxerr = max(maxerr, abs(dT_dϕ - dT_exact))
    end

    return maxerr
  end

  # -------------------------------------------------------------
  # Run convergence tests for each direction
  # -------------------------------------------------------------

  Ns = [32, 64, 128, 256]
  hs = 1 ./ Ns

  # --- radial test ---
  errs_r = [radial_derivative_error(N) for N in Ns]
  p_r = convergence_rate(hs, errs_r)
  @info "Radial derivative convergence rate p ≈ $p_r"
  @test abs(p_r - 2) < 0.3

  # --- theta test ---
  errs_θ = [theta_derivative_error(N) for N in Ns]
  p_θ = convergence_rate(hs, errs_θ)
  @info "Theta derivative convergence rate p ≈ $p_θ"
  @test abs(p_θ - 2) < 0.3

  # --- phi test ---
  errs_ϕ = [phi_derivative_error(N) for N in Ns]
  p_ϕ = convergence_rate(hs, errs_ϕ)
  @info "Phi derivative convergence rate p ≈ $p_ϕ"
  @test abs(p_ϕ - 2) < 0.3

  # --- Plot ---
  plot(
    hs,
    errs_r;
    xscale=:log10,
    yscale=:log10,
    marker=:circle,
    label="radial error",
    lw=2,
    xlabel="h",
    ylabel="L∞ error",
  )
  plot!(hs, errs_θ; marker=:diamond, label="theta error", lw=2)
  plot!(hs, errs_ϕ; marker=:square, label="phi error", lw=2)

  # reference h²
  plot!(hs, hs .^ 2 * errs_r[end] / hs[end]^2; lw=2, ls=:dash, label="h²")

  savefig("spherical_derivative_convergence_r_theta_phi.png")
end
