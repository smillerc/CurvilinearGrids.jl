module MetricDiscretizationSchemes

using Polyester
using ChunkSplitters, KernelAbstractions
using OffsetArrays, StaticArrays
using .Threads, LinearAlgebra

export update_metrics!

include("6th_order_explicit/MonotoneExplicit6thOrderScheme.jl")
using .MonotoneExplicit6thOrderScheme
export MonotoneExplicit6thOrderDiscretization
export update_edge_conserved_metrics!, update_cell_center_metrics!

function update_metrics!(scheme, centroids, cell_center_metrics, edge_metrics, domain)
  # update_cell_center_metrics!(scheme, cell_center_metrics, centroids, domain)
  update_edge_conserved_metrics!(scheme, edge_metrics, centroids, domain)

  return nothing
end

# include("tuple_utils.jl")
# include("indexing_fun.jl")

# export MEG6Scheme, update_metrics!
# export J
# export ξx, ηx, ζx, ξy, ηy, ζy, ξz, ηz, ζz

# export xξ, xη, xζ, yξ, yη, yζ, zξ, zη, zζ
# export ξ̂x, η̂x, ζ̂x, ξ̂y, η̂y, ζ̂y, ξ̂z, η̂z, ζ̂z

# export ξ̂xᵢ₊½, η̂xᵢ₊½, ζ̂xᵢ₊½, ξ̂yᵢ₊½, η̂yᵢ₊½, ζ̂yᵢ₊½, ξ̂zᵢ₊½, η̂zᵢ₊½, ζ̂zᵢ₊½
# export ξ̂xⱼ₊½, η̂xⱼ₊½, ζ̂xⱼ₊½, ξ̂yⱼ₊½, η̂yⱼ₊½, ζ̂yⱼ₊½, ξ̂zⱼ₊½, η̂zⱼ₊½, ζ̂zⱼ₊½
# export ξ̂xₖ₊½, η̂xₖ₊½, ζ̂xₖ₊½, ξ̂yₖ₊½, η̂yₖ₊½, ζ̂yₖ₊½, ξ̂zₖ₊½, η̂zₖ₊½, ζ̂zₖ₊½

# # Convienient getter functions
# @inline @views J(m::AbstractMetricScheme) = m.Jᵢ
# @inline @views ξx(m::AbstractMetricScheme) = m.∂ξ∂xᵢ.∂ξ₁.∂x₁
# @inline @views ηx(m::AbstractMetricScheme) = m.∂ξ∂xᵢ.∂ξ₂.∂x₁
# @inline @views ζx(m::AbstractMetricScheme) = m.∂ξ∂xᵢ.∂ξ₃.∂x₁
# @inline @views ξy(m::AbstractMetricScheme) = m.∂ξ∂xᵢ.∂ξ₁.∂x₂
# @inline @views ηy(m::AbstractMetricScheme) = m.∂ξ∂xᵢ.∂ξ₂.∂x₂
# @inline @views ζy(m::AbstractMetricScheme) = m.∂ξ∂xᵢ.∂ξ₃.∂x₂
# @inline @views ξz(m::AbstractMetricScheme) = m.∂ξ∂xᵢ.∂ξ₁.∂x₃
# @inline @views ηz(m::AbstractMetricScheme) = m.∂ξ∂xᵢ.∂ξ₂.∂x₃
# @inline @views ζz(m::AbstractMetricScheme) = m.∂ξ∂xᵢ.∂ξ₃.∂x₃

# @inline @views xξ(m::AbstractMetricScheme) = m.∂x∂ξᵢ.∂x₁.∂ξ₁
# @inline @views xη(m::AbstractMetricScheme) = m.∂x∂ξᵢ.∂x₁.∂ξ₂
# @inline @views xζ(m::AbstractMetricScheme) = m.∂x∂ξᵢ.∂x₁.∂ξ₃
# @inline @views yξ(m::AbstractMetricScheme) = m.∂x∂ξᵢ.∂x₂.∂ξ₁
# @inline @views yη(m::AbstractMetricScheme) = m.∂x∂ξᵢ.∂x₂.∂ξ₂
# @inline @views yζ(m::AbstractMetricScheme) = m.∂x∂ξᵢ.∂x₂.∂ξ₃
# @inline @views zξ(m::AbstractMetricScheme) = m.∂x∂ξᵢ.∂x₃.∂ξ₁
# @inline @views zη(m::AbstractMetricScheme) = m.∂x∂ξᵢ.∂x₃.∂ξ₂
# @inline @views zζ(m::AbstractMetricScheme) = m.∂x∂ξᵢ.∂x₃.∂ξ₃

# @inline @views ξ̂x(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ.∂ξ̂₁.∂x₁
# @inline @views η̂x(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ.∂ξ̂₂.∂x₁
# @inline @views ζ̂x(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ.∂ξ̂₃.∂x₁
# @inline @views ξ̂y(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ.∂ξ̂₁.∂x₂
# @inline @views η̂y(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ.∂ξ̂₂.∂x₂
# @inline @views ζ̂y(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ.∂ξ̂₃.∂x₂
# @inline @views ξ̂z(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ.∂ξ̂₁.∂x₃
# @inline @views η̂z(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ.∂ξ̂₂.∂x₃
# @inline @views ζ̂z(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ.∂ξ̂₃.∂x₃

# @inline @views ξ̂xᵢ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₁.∂x₁.i₊½
# @inline @views η̂xᵢ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₁.i₊½
# @inline @views ζ̂xᵢ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₃.∂x₁.i₊½
# @inline @views ξ̂yᵢ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₁.∂x₂.i₊½
# @inline @views η̂yᵢ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₂.i₊½
# @inline @views ζ̂yᵢ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₃.∂x₂.i₊½
# @inline @views ξ̂zᵢ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₁.∂x₃.i₊½
# @inline @views η̂zᵢ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₃.i₊½
# @inline @views ζ̂zᵢ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₃.∂x₃.i₊½
# @inline @views ξ̂xⱼ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₁.∂x₁.j₊½
# @inline @views η̂xⱼ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₁.j₊½
# @inline @views ζ̂xⱼ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₃.∂x₁.j₊½
# @inline @views ξ̂yⱼ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₁.∂x₂.j₊½
# @inline @views η̂yⱼ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₂.j₊½
# @inline @views ζ̂yⱼ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₃.∂x₂.j₊½
# @inline @views ξ̂zⱼ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₁.∂x₃.j₊½
# @inline @views η̂zⱼ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₃.j₊½
# @inline @views ζ̂zⱼ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₃.∂x₃.j₊½
# @inline @views ξ̂xₖ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₁.∂x₁.k₊½
# @inline @views η̂xₖ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₁.k₊½
# @inline @views ζ̂xₖ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₃.∂x₁.k₊½
# @inline @views ξ̂yₖ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₁.∂x₂.k₊½
# @inline @views η̂yₖ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₂.k₊½
# @inline @views ζ̂yₖ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₃.∂x₂.k₊½
# @inline @views ξ̂zₖ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₁.∂x₃.k₊½
# @inline @views η̂zₖ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₃.k₊½
# @inline @views ζ̂zₖ₊½(m::AbstractMetricScheme) = m.∂ξ̂∂xᵢ₊½.∂ξ₃.∂x₃.k₊½

# include("monotone_explicit_6th_order_gradient.jl")

end
