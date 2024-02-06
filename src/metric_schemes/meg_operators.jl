using CurvilinearGrids
using UnPack
using BenchmarkTools

include("indexing_fun.jl")

@inline function ∂(ϕ, idx, dim)
  @inbounds begin
    f1 = 3 / 4
    f2 = 3 / 20
    f3 = 1 / 60

    ᵢ₊₁ = up(idx, dim, 1)
    ᵢ₊₂ = up(idx, dim, 2)
    ᵢ₊₃ = up(idx, dim, 3)
    ᵢ₋₁ = down(idx, dim, 1)
    ᵢ₋₂ = down(idx, dim, 2)
    ᵢ₋₃ = down(idx, dim, 3)

    ϕᵢ₊₁ = ϕ(ᵢ₊₁)
    ϕᵢ₊₂ = ϕ(ᵢ₊₂)
    ϕᵢ₊₃ = ϕ(ᵢ₊₃)
    ϕᵢ₋₁ = ϕ(ᵢ₋₁)
    ϕᵢ₋₂ = ϕ(ᵢ₋₂)
    ϕᵢ₋₃ = ϕ(ᵢ₋₃)

    ∂ϕ = (f1 * (ϕᵢ₊₁ - ϕᵢ₋₁) - f2 * (ϕᵢ₊₂ - ϕᵢ₋₂) + f3 * (ϕᵢ₊₃ - ϕᵢ₋₃))
  end
  return ∂ϕ
end

@inline function ∂²(ϕ, idx, dim)
  @inbounds begin
    f1 = 1 / 2

    ᵢ₊₁ = up(idx, dim, 1)
    ᵢ₋₁ = down(idx, dim, 1)
    ϕᵢ₊₁ = ϕ(ᵢ₊₁)
    ϕᵢ = ϕ(idx)
    ϕᵢ₋₁ = ϕ(ᵢ₋₁)

    ∂ϕᵢ₊₁ = ∂(ϕ, ᵢ₊₁, dim)
    ∂ϕᵢ₋₁ = ∂(ϕ, ᵢ₋₁, dim)

    ∂²ϕ = (2(ϕᵢ₊₁ - 2ϕᵢ + ϕᵢ₋₁) - f1 * (∂ϕᵢ₊₁ - ∂ϕᵢ₋₁))
  end
  return ∂²ϕ
end

function ∂ϕ(ϕ, idx, dim)
  @inbounds begin
    ᵢ₊₁ = up(idx, dim, 1)
    ᵢ₋₁ = down(idx, dim, 1)

    xᵢ = ϕ(idx)
    xᵢ₊₁ = ϕ(ᵢ₊₁)
    xᵢ₋₁ = ϕ(ᵢ₋₁)

    xᴸᵢ₊½ = xᵢ + 0.5∂(ϕ, idx, dim) + ∂²(ϕ, idx, dim) / 12
    xᴿᵢ₊½ = xᵢ₊₁ - 0.5∂(ϕ, ᵢ₊₁, dim) + ∂²(ϕ, ᵢ₊₁, dim) / 12
    xᵢ₊½ = (xᴿᵢ₊½ + xᴸᵢ₊½) / 2

    xᴸᵢ₋½ = xᵢ₋₁ + 0.5∂(ϕ, ᵢ₋₁, dim) + ∂²(ϕ, ᵢ₋₁, dim) / 12
    xᴿᵢ₋½ = xᵢ - 0.5∂(ϕ, idx, dim) + ∂²(ϕ, idx, dim) / 12

    xᵢ₋½ = (xᴿᵢ₋½ + xᴸᵢ₋½) / 2
  end
  return xᵢ₊½ - xᵢ₋½
end

function ϕᵢ₊½(ϕ, idx, dim)
  @inbounds begin
    ᵢ₊₁ = up(idx, dim, 1)
    xᵢ = ϕ(idx)
    xᵢ₊₁ = ϕ(ᵢ₊₁)
    xᴸᵢ₊½ = xᵢ + 0.5∂(ϕ, idx, dim) + ∂²(ϕ, idx, dim) / 12
    xᴿᵢ₊½ = xᵢ₊₁ - 0.5∂(ϕ, ᵢ₊₁, dim) + ∂²(ϕ, ᵢ₊₁, dim) / 12
    xᵢ₊½ = (xᴿᵢ₊½ + xᴸᵢ₊½) / 2
  end
  return xᵢ₊½
end

i, j, k = (10, 20, 15)
idx = CartesianIndex(i, j, k)

function metrics()
  function wavy_grid(ni, nj, nk)
    Lx = Ly = Lz = 4.0

    xmin = -Lx / 2
    ymin = -Ly / 2
    zmin = -Lz / 2

    Δx0 = Lx / ni
    Δy0 = Ly / nj
    Δz0 = Lz / nk

    x(i, j, k) = xmin + Δx0 * ((i - 1) + sinpi((j - 1) * Δy0) * sinpi((k - 1) * Δz0))
    y(i, j, k) = ymin + Δy0 * ((j - 1) + sinpi((k - 1) * Δz0) * sinpi((i - 1) * Δx0))
    z(i, j, k) = zmin + Δz0 * ((k - 1) + sinpi((i - 1) * Δx0) * sinpi((j - 1) * Δy0))

    return (x, y, z)
  end

  ni = nj = nk = 20
  nhalo = 5
  x, y, z = wavy_grid(ni, nj, nk)
  mesh = CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)
  ξ, η, ζ = (1, 2, 3)

  # @inline xc(idx) = mesh.x_coord[idx]
  # @inline yc(idx) = mesh.y_coord[idx]
  # @inline zc(idx) = mesh.z_coord[idx]
  xc(idx) = mesh.coord_funcs.x((idx.I .+ 0.5)...)
  yc(idx) = mesh.coord_funcs.y((idx.I .+ 0.5)...)
  zc(idx) = mesh.coord_funcs.z((idx.I .+ 0.5)...)

  xζy(idx) = ∂ϕ(xc, idx, ζ) * yc(idx)
  xηy(idx) = ∂ϕ(xc, idx, η) * yc(idx)
  xξy(idx) = ∂ϕ(xc, idx, ξ) * yc(idx)
  yζz(idx) = ∂ϕ(yc, idx, ζ) * zc(idx)
  yηz(idx) = ∂ϕ(yc, idx, η) * zc(idx)
  yξz(idx) = ∂ϕ(yc, idx, ξ) * zc(idx)
  zζx(idx) = ∂ϕ(zc, idx, ζ) * xc(idx)
  zηx(idx) = ∂ϕ(zc, idx, η) * xc(idx)
  zξx(idx) = ∂ϕ(zc, idx, ξ) * xc(idx)

  xζy_η(idx) = ∂ϕ(xζy, idx, η)
  xζy_ξ(idx) = ∂ϕ(xζy, idx, ξ)
  xηy_ζ(idx) = ∂ϕ(xηy, idx, ζ)
  xηy_ξ(idx) = ∂ϕ(xηy, idx, ξ)
  xξy_ζ(idx) = ∂ϕ(xξy, idx, ζ)
  xξy_η(idx) = ∂ϕ(xξy, idx, η)
  yζz_η(idx) = ∂ϕ(yζz, idx, η)
  yζz_ξ(idx) = ∂ϕ(yζz, idx, ξ)
  yηz_ζ(idx) = ∂ϕ(yηz, idx, ζ)
  yηz_ξ(idx) = ∂ϕ(yηz, idx, ξ)
  yξz_ζ(idx) = ∂ϕ(yξz, idx, ζ)
  yξz_η(idx) = ∂ϕ(yξz, idx, η)
  zζx_η(idx) = ∂ϕ(zζx, idx, η)
  zζx_ξ(idx) = ∂ϕ(zζx, idx, ξ)
  zηx_ζ(idx) = ∂ϕ(zηx, idx, ζ)
  zηx_ξ(idx) = ∂ϕ(zηx, idx, ξ)
  zξx_ζ(idx) = ∂ϕ(zξx, idx, ζ)
  zξx_η(idx) = ∂ϕ(zξx, idx, η)

  ξ̂x(idx) = yηz_ζ(idx) - yζz_η(idx)
  ξ̂y(idx) = zηx_ζ(idx) - zζx_η(idx)
  ξ̂z(idx) = xηy_ζ(idx) - xζy_η(idx)

  η̂x(idx) = yζz_ξ(idx) - yξz_ζ(idx)
  η̂y(idx) = zζx_ξ(idx) - zξx_ζ(idx)
  η̂z(idx) = xζy_ξ(idx) - xξy_ζ(idx)

  ζ̂x(idx) = yξz_η(idx) - yηz_ξ(idx)
  ζ̂y(idx) = zξx_η(idx) - zηx_ξ(idx)
  ζ̂z(idx) = xξy_η(idx) - xηy_ξ(idx)

  ξ̂xᵢ₊½(idx) = ϕᵢ₊½(ξ̂x, idx, ξ)
  η̂xᵢ₊½(idx) = ϕᵢ₊½(η̂x, idx, ξ)
  ζ̂xᵢ₊½(idx) = ϕᵢ₊½(ζ̂x, idx, ξ)

  ξ̂yⱼ₊½(idx) = ϕᵢ₊½(ξ̂y, idx, η)
  η̂yⱼ₊½(idx) = ϕᵢ₊½(η̂y, idx, η)
  ζ̂yⱼ₊½(idx) = ϕᵢ₊½(ζ̂y, idx, η)

  ξ̂zₖ₊½(idx) = ϕᵢ₊½(ξ̂z, idx, ζ)
  η̂zₖ₊½(idx) = ϕᵢ₊½(η̂z, idx, ζ)
  ζ̂zₖ₊½(idx) = ϕᵢ₊½(ζ̂z, idx, ζ)

  ∂ξ̂x∂ξ(idx) = ∂ϕ(ξ̂x, idx, ξ)
  ∂η̂x∂η(idx) = ∂ϕ(η̂x, idx, η)
  ∂ζ̂x∂ζ(idx) = ∂ϕ(ζ̂x, idx, ζ)
  ∂ξ̂y∂ξ(idx) = ∂ϕ(ξ̂y, idx, ξ)
  ∂η̂y∂η(idx) = ∂ϕ(η̂y, idx, η)
  ∂ζ̂y∂ζ(idx) = ∂ϕ(ζ̂y, idx, ζ)
  ∂ξ̂z∂ξ(idx) = ∂ϕ(ξ̂z, idx, ξ)
  ∂η̂z∂η(idx) = ∂ϕ(η̂z, idx, η)
  ∂ζ̂z∂ζ(idx) = ∂ϕ(ζ̂z, idx, ζ)

  # I₁ = ∂ξ̂x∂ξ(idx) + ∂η̂x∂η(idx) + ∂ζ̂x∂ζ(idx)
  # I₂ = ∂ξ̂y∂ξ(idx) + ∂η̂y∂η(idx) + ∂ζ̂y∂ζ(idx)
  # I₃ = ∂ξ̂z∂ξ(idx) + ∂η̂z∂η(idx) + ∂ζ̂z∂ζ(idx)
  # return I₁, I₂, I₃

  return (;
    xc,
    yc,
    zc,
    xζy,
    ∂ξ̂x∂ξ,
    ∂η̂x∂η,
    ∂ζ̂x∂ζ,
    ∂ξ̂y∂ξ,
    ∂η̂y∂η,
    ∂ζ̂y∂ζ,
    ∂ξ̂z∂ξ,
    ∂η̂z∂η,
    ∂ζ̂z∂ζ,
    ξ̂xᵢ₊½,
    η̂xᵢ₊½,
    ζ̂xᵢ₊½,
    ξ̂yⱼ₊½,
    η̂yⱼ₊½,
    ζ̂yⱼ₊½,
    ξ̂zₖ₊½,
    η̂zₖ₊½,
    ζ̂zₖ₊½,
  )
end

# domain = mesh.iterators.cell.domain
m = metrics();

small = CartesianIndices((6:10, 6:10, 6:10))

_xζy = zeros(20, 20, 20);
@time begin
  for idx in small
    # m.ξ̂xᵢ₊½(idx)
    _xζy[idx] = m.xζy(idx)
  end
end

# ic = i + 0.5
# jc = j + 0.5
# kc = k + 0.5

# @benchmark m.xc($idx)
# @code_warntype m.xc(idx)
# @benchmark $x($ic, $jc, $kc)

@benchmark $m.ξ̂yⱼ₊½($idx)

# @time begin
for idx in small
  I₁ = m.∂ξ̂x∂ξ(idx) + m.∂η̂x∂η(idx) + m.∂ζ̂x∂ζ(idx)
  I₂ = m.∂ξ̂y∂ξ(idx) + m.∂η̂y∂η(idx) + m.∂ζ̂y∂ζ(idx)
  I₃ = m.∂ξ̂z∂ξ(idx) + m.∂η̂z∂η(idx) + m.∂ζ̂z∂ζ(idx)
  @show I₁, I₂, I₃
end
# end
# @benchmark ∂ξ̂x∂ξ($idx)

# map(∂ξ̂x∂ξ, domain)

# function ξ̂x(idx)
#   ξ = mesh.cell_center_metrics.ξ[idx]
#   J = mesh.cell_center_metrics.J[idx]
#   return ξ.x * J
# end

# function ξ̂xAD(idx)
#   return conservative_metrics(mesh, idx.I).ξ̂.x
# end

# function η̂x(idx)
#   η = mesh.cell_center_metrics.η[idx]
#   J = mesh.cell_center_metrics.J[idx]
#   return η.x * J
# end

# function ζ̂x(idx)
#   ζ = mesh.cell_center_metrics.ζ[idx]
#   J = mesh.cell_center_metrics.J[idx]
#   return ζ.x * J
# end

# ξx(idx) = mesh.cell_center_metrics.ξ[idx].x
# ηx(idx) = mesh.cell_center_metrics.η[idx].x
# ζx(idx) = mesh.cell_center_metrics.ζ[idx].x
# J(idx) = mesh.cell_center_metrics.J[idx]

# ξ̂xᵢ₊½, ξ̂xᵢ₋½ = ϕ_pm_half(ξ̂x, idx, 1)

# ξ̂xᵢ₋½_AD = conservative_metrics(mesh, (i, j + 0.5, k + 0.5)).ξ̂.x
# ξ̂xᵢ₊½_AD = conservative_metrics(mesh, (i + 1, j + 0.5, k + 0.5)).ξ̂.x
# η̂xⱼ₋½_AD = conservative_metrics(mesh, (i + 0.5, j, k + 0.5)).η̂.x
# η̂xⱼ₊½_AD = conservative_metrics(mesh, (i + 0.5, j + 1, k + 0.5)).η̂.x
# ζ̂ₖ₋½_AD = conservative_metrics(mesh, (i + 0.5, j + 0.5, k)).ζ̂.x
# ζ̂ₖ₊½_AD = conservative_metrics(mesh, (i + 0.5, j + 0.5, k + 1)).ζ̂.x

# (ξ̂xᵢ₋½_AD - ξ̂xᵢ₊½_AD) + (η̂xⱼ₋½_AD - η̂xⱼ₊½_AD) + (ζ̂ₖ₋½_AD - ζ̂ₖ₊½_AD)

# ξ̂x(idx)

# ϕ_pm_half(ξ̂x, idx, 1)
# ϕ_pm_half(J, idx, 1)
# ϕ_pm_half(ξx, idx, 1)
# ϕ_pm_half(ηx, idx, 1)

# I1 = similar(mesh.cell_center_metrics.J)

# for idx in mesh.iterators.cell.domain
#   ξ̂xᵢ₊½, ξ̂xᵢ₋½ = ϕ_pm_half(ξ̂x, idx, 1)
#   η̂xⱼ₊½, η̂xⱼ₋½ = ϕ_pm_half(η̂x, idx, 2)
#   ζxₖ₊½, ζxₖ₋½ = ϕ_pm_half(ζ̂x, idx, 3)

#   I1[idx] = (ξ̂xᵢ₊½ - ξ̂xᵢ₋½) + (η̂xⱼ₊½ - η̂xⱼ₋½) + (ζxₖ₊½ - ζxₖ₋½)
#   # @show I1
# end

# using Plots

# I1d = @views I1[mesh.iterators.cell.domain]
# # ∂.(ξx, CartesianIndices((5:10, 5:5, 10:10)), 1)

# extrema(I1d)
# heatmap(I1d[:, 1, :])