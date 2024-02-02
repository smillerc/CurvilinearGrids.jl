
using BenchmarkTools
using CurvilinearGrids
using MappedArrays
using Polyester
using KernelAbstractions

include("indexing_fun.jl")

function ∂!(∂ϕ, ϕ, dim)
  f1 = 3 / 4
  f2 = 3 / 20
  f3 = 1 / 60

  domain = CartesianIndices(ϕ)
  inner_domain = expand(domain, dim, -4)
  # the ∂ϕ kernel requires i±3, and ultimately we need ∂ϕ
  # at ilo-1 (in order to get i+1/2 at ilo-1), where ilo is the start of the domain. 
  # All this to say, we need 5 halo zones, and start the ∂ϕ
  # at ilo - 2 (where the stencil ultimately extends to ilo - 5).
  # The same is applies to the ihi side
  @batch for i in inner_domain
    ᵢ₊₁ = up(i, dim, 1)
    ᵢ₊₂ = up(i, dim, 2)
    ᵢ₊₃ = up(i, dim, 3)
    ᵢ₋₁ = down(i, dim, 1)
    ᵢ₋₂ = down(i, dim, 2)
    ᵢ₋₃ = down(i, dim, 3)

    ∂ϕ[i] = (f1 * (ϕ[ᵢ₊₁] - ϕ[ᵢ₋₁]) - f2 * (ϕ[ᵢ₊₂] - ϕ[ᵢ₋₂]) + f3 * (ϕ[ᵢ₊₃] - ϕ[ᵢ₋₃]))
  end
  return nothing
end

function ∂²!(∂²ϕ, ∂ϕ, ϕ, dim)
  domain = CartesianIndices(ϕ)
  inner_domain = expand(domain, dim, -4)

  f1 = 1 / 2
  # Similar to the domain extents of the ∂ϕ kernel, ∂²ϕ
  # needs to be computed at ilo-1, where ilo is the start of the domain
  # indices, e.g. ilo:ihi is the domain. This method requires a halo
  # region of 5 cells, so the inner_domain is 4 cells smaller on each
  # dimension.
  @batch for i in inner_domain
    ᵢ₊₁ = up(i, dim, 1)
    ᵢ₋₁ = down(i, dim, 1)

    ∂²ϕ[i] = (2(ϕ[ᵢ₊₁] - 2ϕ[i] + ϕ[ᵢ₋₁]) - f1 * (∂ϕ[ᵢ₊₁] - ∂ϕ[ᵢ₋₁]))
  end
  return nothing
end

@kernel function ∂²_kernel!(
  ∂²ϕ::AbstractArray{T,N},
  @Const(∂ϕ::AbstractArray{T,N}),
  @Const(ϕ::AbstractArray{T,N}),
  @Const(dim::Int),
  @Const(offset::Int)
) where {T,N}
  # domain = CartesianIndices(ϕ)
  # inner_domain = expand(domain, dim, -4)

  # get the global index and adjust for the halo size
  idx = @index(Global, NTuple) .+ offset
  i = CartesianIndex(idx)
  f1 = 1 / 2
  # Similar to the domain extents of the ∂ϕ kernel, ∂²ϕ
  # needs to be computed at ilo-1, where ilo is the start of the domain
  # indices, e.g. ilo:ihi is the domain. This method requires a halo
  # region of 5 cells, so the inner_domain is 4 cells smaller on each
  # dimension.
  ᵢ₊₁ = up(i, dim, 1)
  ᵢ₋₁ = down(i, dim, 1)

  ∂²ϕ[i] = (2(ϕ[ᵢ₊₁] - 2ϕ[i] + ϕ[ᵢ₋₁]) - f1 * (∂ϕ[ᵢ₊₁] - ∂ϕ[ᵢ₋₁]))
  return nothing
end

function ∂x∂ξ!(∂ϕ∂ξ, ∂²ϕ, ∂ϕ, ϕ, dim)
  domain = CartesianIndices(ϕ)
  inner_domain = expand(domain, dim, -1)

  @batch for i in inner_domain
    ᵢ₊₁ = up(i, dim, 1)
    ᵢ₋₁ = down(i, dim, 1)

    xᴸᵢ₊½ = ϕ[i] + 0.5∂ϕ[i] + ∂²ϕ[i] / 12
    xᴿᵢ₊½ = ϕ[ᵢ₊₁] - 0.5∂ϕ[ᵢ₊₁] + ∂²ϕ[ᵢ₊₁] / 12

    xᴸᵢ₋½ = ϕ[ᵢ₋₁] + 0.5∂ϕ[ᵢ₋₁] + ∂²ϕ[ᵢ₋₁] / 12
    xᴿᵢ₋½ = ϕ[i] - 0.5∂ϕ[i] + ∂²ϕ[i] / 12

    xᵢ₊½ = (xᴿᵢ₊½ + xᴸᵢ₊½) / 2
    xᵢ₋½ = (xᴿᵢ₋½ + xᴸᵢ₋½) / 2
    ∂ϕ∂ξ[i] = xᵢ₊½ - xᵢ₋½
  end

  return nothing
end

# --------------------------------------------------------
# --------------------------------------------------------

# Inner partial product, i.e. (∂x/∂ξ)z
function inner∂prod!(
  yξ_z::AbstractArray{T,N}, y::AbstractArray{T,N}, z::AbstractArray{T,N}, ξ::Int
) where {T,N}
  f1 = 2282 / 2880
  f2 = -546 / 2880
  f3 = 89 / 2880
  f4 = -3 / 2880
  f5 = -1 / 2880

  domain = CartesianIndices(y)
  inner_domain = expand(domain, -5)

  @batch for i in inner_domain
    ᵢ₊₁ = up(i, ξ, 1)
    ᵢ₊₂ = up(i, ξ, 2)
    ᵢ₊₃ = up(i, ξ, 3)
    ᵢ₊₄ = up(i, ξ, 4)
    ᵢ₊₅ = up(i, ξ, 5)
    ᵢ₋₁ = down(i, ξ, 1)
    ᵢ₋₂ = down(i, ξ, 2)
    ᵢ₋₃ = down(i, ξ, 3)
    ᵢ₋₄ = down(i, ξ, 4)
    ᵢ₋₅ = down(i, ξ, 5)

    yξ_z[i] =
      (
        f1 * (y[ᵢ₊₁] - y[ᵢ₋₁]) + #
        f2 * (y[ᵢ₊₂] - y[ᵢ₋₂]) + #
        f3 * (y[ᵢ₊₃] - y[ᵢ₋₃]) + #
        f4 * (y[ᵢ₊₄] - y[ᵢ₋₄]) + #
        f5 * (y[ᵢ₊₅] - y[ᵢ₋₅])
      ) * z[i]
  end
end

function conserved!(
  ξ̂x::AbstractArray{T,N}, # conserved metric, could be ξ̂x, η̂z, ζ̂y, etc..
  y_ηz::AbstractArray{T,N}, # 1st term, (∂y/∂η)⋅z, in the 2-term eq to determine ξ̂x
  ζ::Int, # outer dimension to differentiate the 1st term on, e.g. ξ:i±1/2, η:j±1/2, ζ:k±1/2
  yζ_z::AbstractArray{T,N}, # 2nd term, (∂y/∂ζ)⋅z, in the 2-term eq to determine ξ̂x
  η::Int, # outer dimension to differentiate the 2nd term on, e.g. ξ:i±1/2, η:j±1/2, ζ:k±1/2
) where {T,N}
  a = 2282 / 2880
  b = -546 / 2880
  c = 89 / 2880
  d = -3 / 2880
  e = -1 / 2880

  domain = CartesianIndices(ξ̂x)
  inner_domain = expand(domain, -5) # shrink due to the stencil size

  @batch for i in inner_domain
    ₖ₊₁ = up(i, ζ, 1)
    ₖ₊₂ = up(i, ζ, 2)
    ₖ₊₃ = up(i, ζ, 3)
    ₖ₊₄ = up(i, ζ, 4)
    ₖ₊₅ = up(i, ζ, 5)
    ₖ₋₁ = down(i, ζ, 1)
    ₖ₋₂ = down(i, ζ, 2)
    ₖ₋₃ = down(i, ζ, 3)
    ₖ₋₄ = down(i, ζ, 4)
    ₖ₋₅ = down(i, ζ, 5)

    ⱼ₊₁ = up(i, η, 1)
    ⱼ₊₂ = up(i, η, 2)
    ⱼ₊₃ = up(i, η, 3)
    ⱼ₊₄ = up(i, η, 4)
    ⱼ₊₅ = up(i, η, 5)
    ⱼ₋₁ = down(i, η, 1)
    ⱼ₋₂ = down(i, η, 2)
    ⱼ₋₃ = down(i, η, 3)
    ⱼ₋₄ = down(i, η, 4)
    ⱼ₋₅ = down(i, η, 5)

    ξ̂x[i] =
      (
        a * (y_ηz[ₖ₊₁] - y_ηz[ₖ₋₁]) +
        b * (y_ηz[ₖ₊₂] - y_ηz[ₖ₋₂]) +
        c * (y_ηz[ₖ₊₃] - y_ηz[ₖ₋₃]) +
        d * (y_ηz[ₖ₊₄] - y_ηz[ₖ₋₄]) +
        e * (y_ηz[ₖ₊₅] - y_ηz[ₖ₋₅])
      ) - (
        a * (yζ_z[ⱼ₊₁] - yζ_z[ⱼ₋₁]) +
        b * (yζ_z[ⱼ₊₂] - yζ_z[ⱼ₋₂]) +
        c * (yζ_z[ⱼ₊₃] - yζ_z[ⱼ₋₃]) +
        d * (yζ_z[ⱼ₊₄] - yζ_z[ⱼ₋₄]) +
        e * (yζ_z[ⱼ₊₅] - yζ_z[ⱼ₋₅])
      )
  end

  return nothing
end

function to_edge!(ϕᵢ₊½::AbstractArray{T,N}, ϕ::AbstractArray{T,N}, dim::Int) where {T,N}
  a = 1821 / 2880
  b = -461 / 2880
  c = 85 / 2880
  d = -4 / 2880
  e = -1 / 2880

  domain = CartesianIndices(ϕ)
  inner_domain = expand(domain, -5) # shrink due to the stencil size

  @batch for i in inner_domain
    ᵢ₊₁ = up(i, dim, 1)
    ᵢ₊₂ = up(i, dim, 2)
    ᵢ₊₃ = up(i, dim, 3)
    ᵢ₊₄ = up(i, dim, 4)
    ᵢ₊₅ = up(i, dim, 5)
    ᵢ₋₁ = down(i, dim, 1)
    ᵢ₋₂ = down(i, dim, 2)
    ᵢ₋₃ = down(i, dim, 3)
    ᵢ₋₄ = down(i, dim, 4)

    ϕᵢ₊½[i] = (
      a * (ϕ[ᵢ₊₁] + ϕ[i]) +
      b * (ϕ[ᵢ₊₂] + ϕ[ᵢ₋₁]) +
      c * (ϕ[ᵢ₊₃] + ϕ[ᵢ₋₂]) +
      d * (ϕ[ᵢ₊₄] + ϕ[ᵢ₋₃]) +
      e * (ϕ[ᵢ₊₅] + ϕ[ᵢ₋₄])
    )
  end

  return nothing
end

function conserved_metric(
  ξ̂_x::AbstractArray{T,N}, (y, z), (η, ζ)::NTuple{2,Int}, cache
) where {T,N}
  yη_z = @views cache.mixed_deriv1
  yζ_z = @views cache.mixed_deriv2

  inner∂prod!(yη_z, y, z, η) # -> yη_z
  inner∂prod!(yζ_z, y, z, ζ) # -> yζ_z
  conserved!(ξ̂_x, yη_z, ζ, yζ_z, η) # -> ξ̂x

  return ξ̂_x
end
# to_edge!(ξ̂xᵢ₊½, ξ̂x, ξ) # -> ξ̂xᵢ₊½

# --------------------------------------------------------
# --------------------------------------------------------

# # Outer partial of product, i.e. ∂((∂x/∂ξ)z)/∂ζ
# function outermixed∂!(
#   yηz_ζ::AbstractArray{T,N},  # ∂((∂x/∂ξ)z)/∂ζ
#   yηz::AbstractArray{T,N},  # (∂x/∂ξ)z
#   ζ::Int, # outer partial dimension
# ) where {T,N}
#   f1 = 2282 / 2880
#   f1 = -546 / 2880
#   f1 = 89 / 2880
#   f1 = -3 / 2880
#   f1 = -1 / 2880

#   domain = CartesianIndices(yηz)
#   inner_domain = expand(domain, dim, -5)

#   @batch for i in inner_domain
#     ᵢ₊₁ = up(i, ζ, 1)
#     ᵢ₊₂ = up(i, ζ, 2)
#     ᵢ₊₃ = up(i, ζ, 3)
#     ᵢ₊₄ = up(i, ζ, 4)
#     ᵢ₊₅ = up(i, ζ, 5)
#     ᵢ₋₁ = down(i, ζ, 1)
#     ᵢ₋₂ = down(i, ζ, 2)
#     ᵢ₋₃ = down(i, ζ, 3)
#     ᵢ₋₄ = down(i, ζ, 4)
#     ᵢ₋₅ = down(i, ζ, 5)

#     yηz_ζ[i] = (
#       f1 * (yηz[ᵢ₊₁] - yηz[ᵢ₋₁]) + #
#       f2 * (yηz[ᵢ₊₂] - yηz[ᵢ₋₂]) + #
#       f3 * (yηz[ᵢ₊₃] - yηz[ᵢ₋₃]) + #
#       f4 * (yηz[ᵢ₊₄] - yηz[ᵢ₋₄]) + #
#       f5 * (yηz[ᵢ₊₅] - yηz[ᵢ₋₅])  #
#     )
#   end

#   return nothing
# end

# function ϕ_to_edge(ϕᵢ₊½, ϕ, dim, cache)
#   ∂²ϕ = @views cache.∂²ϕ
#   ∂ϕ = @views cache.∂ϕ

#   ∂!(∂ϕ, ϕ, dim)
#   ∂²!(∂²ϕ, ∂ϕ, ϕ, dim)

#   domain = CartesianIndices(ϕ)
#   inner_domain = expand(domain, dim, -1)

#   @batch for i in inner_domain
#     ᵢ₊₁ = up(i, dim, 1)
#     xᴸᵢ₊½ = ϕ[i] + 0.5∂ϕ[i] + ∂²ϕ[i] / 12
#     xᴿᵢ₊½ = ϕ[ᵢ₊₁] - 0.5∂ϕ[ᵢ₊₁] + ∂²ϕ[ᵢ₊₁] / 12
#     ϕᵢ₊½[i] = (xᴿᵢ₊½ + xᴸᵢ₊½) / 2
#   end

#   return nothing
# end

# function mixed_∂!(yηz_ζ, y, η, z, ζ, cache_arrays)
#   ∂²ϕ = @views cache_arrays.∂²ϕ
#   ∂ϕ = @views cache_arrays.∂ϕ
#   yη = @views cache_arrays.inv_metric

#   #y_η
#   ∂!(∂ϕ, y, η)
#   ∂²!(∂²ϕ, ∂ϕ, y, η)
#   ∂x∂ξ!(yη, ∂²ϕ, ∂ϕ, y, η)
#   yηz = mappedarray(*, yη, z)
#   # (y_η⋅z)_ζ
#   ∂!(∂ϕ, yηz, ζ)
#   ∂²!(∂²ϕ, ∂ϕ, yηz, ζ)
#   ∂x∂ξ!(yηz_ζ, ∂²ϕ, ∂ϕ, yηz, ζ)
#   return nothing
# end

# function conserved_metric(ξ̂x, (y, z), (η, ζ), cache)
#   yηz_ζ = @views cache.mixed_deriv1
#   yζz_η = @views cache.mixed_deriv2
#   mixed_∂!(yηz_ζ, y, η, z, ζ, cache)
#   mixed_∂!(yζz_η, y, ζ, z, η, cache)

#   @. ξ̂x = yηz_ζ - yζz_η

#   return nothing
# end

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

# --------------------------------------------------------
ni = nj = nk = 20
nhalo = 5
x, y, z = wavy_grid(ni, nj, nk)
mesh = CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo);

xc = zeros(size(mesh.iterators.cell.full));
yc = zeros(size(mesh.iterators.cell.full));
zc = zeros(size(mesh.iterators.cell.full));

for k in axes(xc, 3)
  for j in axes(xc, 2)
    for i in axes(xc, 1)
      xc[i, j, k] = mesh.x_func(i - nhalo + 0.5, j - nhalo + 0.5, k - nhalo + 0.5)
      yc[i, j, k] = mesh.y_func(i - nhalo + 0.5, j - nhalo + 0.5, k - nhalo + 0.5)
      zc[i, j, k] = mesh.z_func(i - nhalo + 0.5, j - nhalo + 0.5, k - nhalo + 0.5)
    end
  end
end

# --------------------------------------------------------

ξ, η, ζ = (1, 2, 3)
cache = (
  ∂²ϕ=zeros(size(xc)),
  ∂ϕ=zeros(size(xc)),
  mixed_deriv1=zeros(size(xc)),
  mixed_deriv2=zeros(size(xc)),
);

domain = mesh.iterators.cell.domain

begin
  ξ̂_x = zeros(size(xc))
  ξ̂xᵢ₊½ = similar(yc)
  conserved_metric(ξ̂_x, (yc, zc), (η, ζ), cache)
  to_edge!(ξ̂xᵢ₊½, ξ̂_x, ξ)

  # η̂_x = zeros(size(xc))
  # ζ̂_x = zeros(size(xc))
  # ξ̂_y = zeros(size(xc))
  # η̂_y = zeros(size(xc))
  # ζ̂_y = zeros(size(xc))
  # ξ̂_z = zeros(size(xc))
  # η̂_z = zeros(size(xc))
  # ζ̂_z = zeros(size(xc))

  # conserved_metric(η̂_x, (yc, zc), (ζ, ξ), cache)
  # conserved_metric(ζ̂_x, (yc, zc), (ξ, η), cache)

  # conserved_metric(ξ̂_y, (zc, xc), (η, ζ), cache)
  # conserved_metric(η̂_y, (zc, xc), (ζ, ξ), cache)
  # conserved_metric(ζ̂_y, (zc, xc), (ξ, η), cache)

  # conserved_metric(ξ̂_z, (xc, yc), (η, ζ), cache)
  # conserved_metric(η̂_z, (xc, yc), (ζ, ξ), cache)
  # conserved_metric(ζ̂_z, (xc, yc), (ξ, η), cache)

  # η̂xⱼ₊½ = similar(yc)
  # ζ̂xₖ₊½ = similar(yc)

  # ξ̂yᵢ₊½ = similar(yc)
  # η̂yⱼ₊½ = similar(yc)
  # ζ̂yₖ₊½ = similar(yc)

  # ξ̂zᵢ₊½ = similar(yc)
  # η̂zⱼ₊½ = similar(yc)
  # ζ̂zₖ₊½ = similar(yc)

  # to_edge!(η̂xⱼ₊½, η̂_x, η)
  # to_edge!(ζ̂xₖ₊½, ζ̂_x, ζ)

  # to_edge!(ξ̂yᵢ₊½, ξ̂_y, ξ)
  # to_edge!(η̂yⱼ₊½, η̂_y, η)
  # to_edge!(ζ̂yₖ₊½, ζ̂_y, ζ)

  # to_edge!(ξ̂zᵢ₊½, ξ̂_z, ξ)
  # to_edge!(η̂zⱼ₊½, η̂_z, η)
  # to_edge!(ζ̂zₖ₊½, ζ̂_z, ζ)
end

begin
  ξ̂_x = zeros(size(xc))
  η̂_x = zeros(size(xc))
  ζ̂_x = zeros(size(xc))
  ξ̂_y = zeros(size(xc))
  η̂_y = zeros(size(xc))
  ζ̂_y = zeros(size(xc))
  ξ̂_z = zeros(size(xc))
  η̂_z = zeros(size(xc))
  ζ̂_z = zeros(size(xc))

  conserved_metric(ξ̂_x, (yc, zc), (η, ζ), cache)
  conserved_metric(η̂_x, (yc, zc), (ζ, ξ), cache)
  conserved_metric(ζ̂_x, (yc, zc), (ξ, η), cache)

  conserved_metric(ξ̂_y, (zc, xc), (η, ζ), cache)
  conserved_metric(η̂_y, (zc, xc), (ζ, ξ), cache)
  conserved_metric(ζ̂_y, (zc, xc), (ξ, η), cache)

  conserved_metric(ξ̂_z, (xc, yc), (η, ζ), cache)
  conserved_metric(η̂_z, (xc, yc), (ζ, ξ), cache)
  conserved_metric(ζ̂_z, (xc, yc), (ξ, η), cache)

  ξ̂xᵢ₊½ = similar(yc)
  η̂xⱼ₊½ = similar(yc)
  ζ̂xₖ₊½ = similar(yc)

  ξ̂yᵢ₊½ = similar(yc)
  η̂yⱼ₊½ = similar(yc)
  ζ̂yₖ₊½ = similar(yc)

  ξ̂zᵢ₊½ = similar(yc)
  η̂zⱼ₊½ = similar(yc)
  ζ̂zₖ₊½ = similar(yc)

  to_edge!(ξ̂xᵢ₊½, ξ̂_x, ξ)
  to_edge!(η̂xⱼ₊½, η̂_x, η)
  to_edge!(ζ̂xₖ₊½, ζ̂_x, ζ)

  to_edge!(ξ̂yᵢ₊½, ξ̂_y, ξ)
  to_edge!(η̂yⱼ₊½, η̂_y, η)
  to_edge!(ζ̂yₖ₊½, ζ̂_y, ζ)

  to_edge!(ξ̂zᵢ₊½, ξ̂_z, ξ)
  to_edge!(η̂zⱼ₊½, η̂_z, η)
  to_edge!(ζ̂zₖ₊½, ζ̂_z, ζ)
end

for idx in domain
  i, j, k = idx.I
  ∂ξ̂x∂ξ = ξ̂xᵢ₊½[i, j, k] - ξ̂xᵢ₊½[i - 1, j, k]
  ∂η̂x∂η = η̂xⱼ₊½[i, j, k] - η̂xⱼ₊½[i, j - 1, k]
  ∂ζ̂x∂ζ = ζ̂xₖ₊½[i, j, k] - ζ̂xₖ₊½[i, j, k - 1]

  ∂ξ̂y∂ξ = ξ̂yᵢ₊½[i, j, k] - ξ̂yᵢ₊½[i - 1, j, k]
  ∂η̂y∂η = η̂yⱼ₊½[i, j, k] - η̂yⱼ₊½[i, j - 1, k]
  ∂ζ̂y∂ζ = ζ̂yₖ₊½[i, j, k] - ζ̂yₖ₊½[i, j, k - 1]

  ∂ξ̂z∂ξ = ξ̂zᵢ₊½[i, j, k] - ξ̂zᵢ₊½[i - 1, j, k]
  ∂η̂z∂η = η̂zⱼ₊½[i, j, k] - η̂zⱼ₊½[i, j - 1, k]
  ∂ζ̂z∂ζ = ζ̂zₖ₊½[i, j, k] - ζ̂zₖ₊½[i, j, k - 1]

  I₁ = ∂ξ̂x∂ξ + ∂η̂x∂η + ∂ζ̂x∂ζ
  I₂ = ∂ξ̂y∂ξ + ∂η̂y∂η + ∂ζ̂y∂ζ
  I₃ = ∂ξ̂z∂ξ + ∂η̂z∂η + ∂ζ̂z∂ζ

  @show I₁, I₂, I₃
end

# using Plots

# ξ̂x_dom = @views ξ̂x[domain]

# heatmap(ξ̂x_dom[:, 10, :])

# # serial ~5 ms
# # 16 threads ~195 μs
# @benchmark $mixed_∂!($yηz_ζ, $y, $η, $z, $ζ, $cache)

# #yηz_ζ
# @benchmark $mixed_∂!($yηz_ζ, $y, $η, $cache)

# @benchmark $up($idx, 1, 1)

# @benchmark $∂!($∂A, $A, $1)
# @benchmark $e∂!($∂A, $A, $1)

# @code_warntype ∂!(∂A, A, 1)

# domain = CartesianIndices(y)
# inner_domain = expand(domain, 1, -3)
# @benchmark $d∂!($∂A, $A, $inner_domain)

# @benchmark $∂²!($∂A, $A, $A, $1)

function ∂KA!(∂ϕ, ϕ, dim)
  backend = get_backend(ϕ)

  domain = CartesianIndices(ϕ)
  inner_domain = expand(domain, -5)
  nd_range = size(inner_domain)
  # @show nd_range
  # kernel = ∂_kernel!(backend)
  kernel = ∂_kernel2!(backend)

  offset = 4
  kernel(∂ϕ, ϕ, dim, offset; ndrange=nd_range)

  KernelAbstractions.synchronize(backend)
  return nothing
end

@kernel function ∂_kernel!(
  ∂ϕ::AbstractArray{T,N},
  @Const(ϕ::AbstractArray{T,N}),
  @Const(dim::Int),
  @Const(offset::Int)
) where {T,N}
  idx = @index(Global, Cartesian)

  i = CartesianIndex(idx.I .+ offset)

  # convert to Float32 if need be
  f1 = 3 / 4
  f2 = 3 / 20
  f3 = 1 / 60

  ᵢ₊₁ = up(i, dim, 1)
  ᵢ₊₂ = up(i, dim, 2)
  ᵢ₊₃ = up(i, dim, 3)
  ᵢ₋₁ = down(i, dim, 1)
  ᵢ₋₂ = down(i, dim, 2)
  ᵢ₋₃ = down(i, dim, 3)

  ∂ϕ[i] = (
    f1 * (ϕ[ᵢ₊₁] - ϕ[ᵢ₋₁]) -  # 
    f2 * (ϕ[ᵢ₊₂] - ϕ[ᵢ₋₂]) +  # 
    f3 * (ϕ[ᵢ₊₃] - ϕ[ᵢ₋₃]) # 
  )
  return nothing
end

@kernel function ∂_kernel2!(
  ∂ϕ::AbstractArray{T,N},
  @Const(ϕ::AbstractArray{T,N}),
  @Const(dim::Int),
  @Const(offset::Int)
) where {T,N}
  idx = @index(Global, Cartesian)
  i, j, k = idx.I .+ offset
  # i = CartesianIndex(idx.I .+ offset)

  # convert to Float32 if need be
  f1 = 3 / 4
  f2 = 3 / 20
  f3 = 1 / 60

  ∂ϕ[i, j, k] = (
    f1 * (ϕ[i, j + 1, k] - ϕ[i, j - 1, k]) -  # 
    f2 * (ϕ[i, j + 2, k] - ϕ[i, j - 2, k]) +  # 
    f3 * (ϕ[i, j + 3, k] - ϕ[i, j - 3, k]) # 
  )
  return nothing
end

using Adapt
using StaticArrays

if Base.find_package("CUDA") !== nothing
  using CUDA
  using CUDA.CUDAKernels
  const backend = CUDABackend()
  CUDA.allowscalar(false)
else
  const backend = CPU()
end

CI = CartesianIndices((400, 200, 200))

T = Float64
∂A = KernelAbstractions.ones(backend, T, size(CI)...);
A = KernelAbstractions.ones(backend, T, size(CI)...);

# Av = Array{SVector{6,T}}(undef, size(CI)...);
# ∂Av = Array{SVector{6,T}}(undef, size(CI)...);
Av = CuArray{SVector{6,T}}(undef, size(CI)...);
∂Av = CuArray{SVector{6,T}}(undef, size(CI)...);

@benchmark begin
  $∂KA!($∂A, $A, $2)
  # GPU
  # 1.575 ms F64
  # 1.169 ms F32

  # CPU
  # 191 μs
end

@benchmark begin
  $∂KA!($∂Av, $Av, $2)
  # 4.418 ms with SVector{4,Float64}
  # 5.532 ms with SVector{6,Float64}

  # 2.239 ms F32 with SVector{4,Float32}
  # 2.592 ms F32 with SVector{6,Float32}

  # CPU
  # 154 μs
end