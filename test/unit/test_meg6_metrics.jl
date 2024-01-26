using Test

using CurvilinearGrids
using StaticArrays
using Test
using BenchmarkTools

function wavy_grid(nx, ny)
  x0, x1 = (0, 1)
  y0, y1 = (0, 1)
  a0 = 0.1

  function x(i, j)
    x1d = x0 + (x1 - x0) * ((i - 1) / (nx - 1))
    y1d = y0 + (y1 - y0) * ((j - 1) / (ny - 1))
    return x1d + a0 * sin(2 * pi * x1d) * sin(2 * pi * y1d)
  end

  function y(i, j)
    x1d = x0 + (x1 - x0) * ((i - 1) / (nx - 1))
    y1d = y0 + (y1 - y0) * ((j - 1) / (ny - 1))
    return y1d + a0 * sin(2 * pi * x1d) * sin(2 * pi * y1d)
  end

  return (x, y)
end

ni, nj = (41, 41)
nhalo = 4
x, y = wavy_grid(ni, nj)
mesh = CurvilinearGrid2D(x, y, (ni, nj), nhalo)

xy = coords(mesh)
xn = @view xy[1, :, :]
yn = @view xy[2, :, :]

metrics = MEG6Scheme(cellsize_withhalo(mesh))
domain = mesh.iterators.cell.domain
update_metrics!(metrics, mesh.x_coord, mesh.y_coord, domain)

display(heatmap(J(metrics); title="J"))
display(heatmap(xξ(metrics); title="xξ"))
display(heatmap(ξx(metrics); title="ξx"))
# display(heatmap(ξ̂x(metrics); title="ξ̂x"))
display(heatmap(ξ̂xᵢ₊½(metrics); title="ξ̂xᵢ₊½"))
# display(heatmap(ξ̂xⱼ₊½(metrics); title="ξ̂xⱼ₊½"))

# @testset "Indexing functions" begin
#   include("../../src/metric_schemes/indexing_fun.jl")
#   domain = CartesianIndices((1:10, 4:8))

#   @test lower_boundary_indices(domain, 1, +1) == CartesianIndices((2:2, 4:8))
#   @test lower_boundary_indices(domain, 1, 0) == CartesianIndices((1:1, 4:8))
#   @test lower_boundary_indices(domain, 1, -1) == CartesianIndices((0:0, 4:8))

#   @test upper_boundary_indices(domain, 2, +1) == CartesianIndices((1:10, 9:9))
#   @test upper_boundary_indices(domain, 2, 0) == CartesianIndices((1:10, 8:8))
#   @test upper_boundary_indices(domain, 2, -1) == CartesianIndices((1:10, 7:7))

#   @test expand(domain, 2, -1) == CartesianIndices((1:10, 5:7))
#   @test expand(domain, 2, 0) == CartesianIndices((1:10, 4:8))
#   @test expand(domain, 2, 2) == CartesianIndices((1:10, 2:10))

#   @test expand(domain, 2) == CartesianIndices((-1:12, 2:10))

#   @test expand_lower(domain, 2) == CartesianIndices((-1:10, 2:8))
#   @test expand_upper(domain, 2) == CartesianIndices((1:12, 4:10))

#   @test up(CartesianIndex((4, 5, 6)), 3, 2) == CartesianIndex(4, 5, 8)
#   @test down(CartesianIndex((4, 5, 6)), 1, 2) == CartesianIndex(2, 5, 6)
#   @test down(CartesianIndex((4, 5, 6)), 0, 2) == CartesianIndex(4, 5, 6)

#   @test plus_minus(CartesianIndex((4, 5, 6)), 2, 2) == CartesianIndices((4:4, 3:7, 6:6))

#   # non-uniform +/-: -1:+2
#   @test plus_minus(CartesianIndex((4, 5, 6)), 2, (1, 2)) ==
#     CartesianIndices((4:4, 4:7, 6:6))
#   @test plus_minus(CartesianIndex((4, 5, 6)), 2, (1, 0)) ==
#     CartesianIndices((4:4, 4:5, 6:6))

#   @test δ(1, CartesianIndex((4, 5, 6))) == CartesianIndex((1, 0, 0))
# end

# using ForwardDiff

# const n_i = 20
# function f(i)
#   return cospi((i - 1) / (n_i - 1))
# end

# ∂f(i) = ForwardDiff.derivative(f, i)
# ∂²f(i) = ForwardDiff.derivative(∂f, i)

# xn = -2:52
# xc = -1.5:1:51.5

# fdata = f.(xc);
# ∂fdata = similar(fdata)
# ∂²fdata = similar(fdata)
# fill!(∂fdata, NaN)
# fill!(∂²fdata, NaN)
# fullCI = CartesianIndices(fdata)
# nhalo = 2
# domain = fullCI[(begin + nhalo):(end - nhalo)]
# ∂ϕ∂ξ = similar(fdata)
# fill!(∂ϕ∂ξ, NaN)
# # CurvilinearGrids.MetricDiscretizationSchemes._get_boundary_gradient(
# #   fdata, ∂fdata, ∂²fdata, 1, domain
# # )

# CurvilinearGrids.MetricDiscretizationSchemes._cell_center_metric!(
#   ∂ϕ∂ξ, fdata, ∂fdata, ∂²fdata, 1, domain
# )

# # CurvilinearGrids.MetricDiscretizationSchemes._boundary_∂!(∂fdata, fdata, 1, domain)
# # CurvilinearGrids.MetricDiscretizationSchemes._boundary_∂²!(
# #   ∂²fdata, fdata, ∂fdata, 1, domain
# # )

# # plot(xn, f.(xn); marker=:square, ms=2)
# # plot(xc, f.(xc); marker=:diamond, ms=2)
# plot(xc, ∂f.(xc); marker=:circle, ms=2, label="∂f")
# plot!(xc, ∂²f.(xc); marker=:circle, ms=2, label="∂²f")
# plot!(xc, ∂fdata; marker=:star, label="∂fdata")
# plot!(xc, ∂²fdata; marker=:star, label="∂²fdata")
# # ylims!(-1, 1)

# err = abs.(∂²fdata .- ∂²f.(xc))#[(begin + 2):(end - 2)]
# plot!(twinx(), xc, err)