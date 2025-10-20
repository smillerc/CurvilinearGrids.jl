
using ForwardDiff, StructArrays, KernelAbstractions, CartesianDomains, UnPack
using DifferentiationInterface

using StaticArrays, LinearAlgebra
using WriteVTK, Polyester, .Threads

abstract type AbstractContinuousCurvilinearGrid1D <: AbstractCurvilinearGrid2D end
abstract type AbstractContinuousCurvilinearGrid2D <: AbstractCurvilinearGrid2D end
abstract type AbstractContinuousCurvilinearGrid3D <: AbstractCurvilinearGrid3D end

struct ContinuousCurvilinearGrid3D{A,B,C,D,E,BE,DBE} <: AbstractContinuousCurvilinearGrid3D
  node_coordinates::A
  centroid_coordinates::A
  mapping_functions::B
  metric_functions_cache::C
  metrics::D
  iterators::E
  backend::BE
  diff_backend::DBE
end

struct MetricCache{F1,F2}
  forward::F1
  inverse::F2
end

function ContinuousCurvilinearGrid3D(
  x::Function,
  y::Function,
  z::Function,
  celldims::NTuple,
  nhalo::Int,
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  T=Float64;
)
  iterators = get_iterators(celldims, nhalo)

  cell_center_metrics, edge_metrics = get_metric_soa(celldims .+ 2nhalo, backend, T)
  metrics = (; cell=cell_center_metrics, edge=edge_metrics)
  xn, yn, zn = node_coordinates(x, y, z, iterators)
  xc, yc, zc = centroid_coordinates(x, y, z, iterators)
  mesh = ContinuousCurvilinearGrid3D(
    (; x=xn, y=yn, z=zn),
    (; x=xc, y=yc, z=zc),
    (; x, y, z),
    MetricCache(x, y, z, diff_backend),
    metrics,
    iterators,
    backend,
    diff_backend,
  )

  compute_cell_metrics!(mesh)
  compute_edge_metrics!(mesh)

  return mesh
end

function MetricCache(x::Function, y::Function, z::Function, backend)
  xξ(i, j, k) = derivative(ξ -> x(ξ, j, k), backend, i)
  xη(i, j, k) = derivative(η -> x(i, η, k), backend, j)
  xζ(i, j, k) = derivative(ζ -> x(i, j, ζ), backend, k)
  yξ(i, j, k) = derivative(ξ -> y(ξ, j, k), backend, i)
  yη(i, j, k) = derivative(η -> y(i, η, k), backend, j)
  yζ(i, j, k) = derivative(ζ -> y(i, j, ζ), backend, k)
  zξ(i, j, k) = derivative(ξ -> z(ξ, j, k), backend, i)
  zη(i, j, k) = derivative(η -> z(i, η, k), backend, j)
  zζ(i, j, k) = derivative(ζ -> z(i, j, ζ), backend, k)

  function jacobian_matrix(i, j, k)
    DifferentiationInterface.jacobian(
      u -> SVector(x(u...), y(u...), z(u...)), backend, @SVector [i, j, k]
    )
  end

  forward_metrics = (;
    jacobian=jacobian_matrix, xξ=xξ, xη=xη, xζ=xζ, yξ=yξ, yη=yη, yζ=yζ, zξ=zξ, zη=zη, zζ=zζ
  )

  x_ξ_y(i, j, k) = xξ(i, j, k) * y(i, j, k)
  x_η_y(i, j, k) = xη(i, j, k) * y(i, j, k)
  x_ζ_y(i, j, k) = xζ(i, j, k) * y(i, j, k)

  y_ξ_z(i, j, k) = yξ(i, j, k) * z(i, j, k)
  y_η_z(i, j, k) = yη(i, j, k) * z(i, j, k)
  y_ζ_z(i, j, k) = yζ(i, j, k) * z(i, j, k)

  z_ξ_x(i, j, k) = zξ(i, j, k) * x(i, j, k)
  z_η_x(i, j, k) = zη(i, j, k) * x(i, j, k)
  z_ζ_x(i, j, k) = zζ(i, j, k) * x(i, j, k)

  _, x_ξ_y_η, x_ξ_y_ζ = cell_center_derivative(x_ξ_y, backend)
  x_η_y_ξ, _, x_η_y_ζ = cell_center_derivative(x_η_y, backend)
  x_ζ_y_ξ, x_ζ_y_η, _ = cell_center_derivative(x_ζ_y, backend)

  _, y_ξ_z_η, y_ξ_z_ζ = cell_center_derivative(y_ξ_z, backend)
  y_η_z_ξ, _, y_η_z_ζ = cell_center_derivative(y_η_z, backend)
  y_ζ_z_ξ, y_ζ_z_η, _ = cell_center_derivative(y_ζ_z, backend)

  _, z_ξ_x_η, z_ξ_x_ζ = cell_center_derivative(z_ξ_x, backend)
  z_η_x_ξ, _, z_η_x_ζ = cell_center_derivative(z_η_x, backend)
  z_ζ_x_ξ, z_ζ_x_η, _ = cell_center_derivative(z_ζ_x, backend)

  ξ̂x(i, j, k) = y_η_z_ζ(i, j, k) − y_ζ_z_η(i, j, k)
  η̂x(i, j, k) = y_ζ_z_ξ(i, j, k) − y_ξ_z_ζ(i, j, k)
  ζ̂x(i, j, k) = y_ξ_z_η(i, j, k) − y_η_z_ξ(i, j, k)
  ξ̂y(i, j, k) = z_η_x_ζ(i, j, k) − z_ζ_x_η(i, j, k)
  η̂y(i, j, k) = z_ζ_x_ξ(i, j, k) − z_ξ_x_ζ(i, j, k)
  ζ̂y(i, j, k) = z_ξ_x_η(i, j, k) − z_η_x_ξ(i, j, k)
  ξ̂z(i, j, k) = x_η_y_ζ(i, j, k) − x_ζ_y_η(i, j, k)
  η̂z(i, j, k) = x_ζ_y_ξ(i, j, k) − x_ξ_y_ζ(i, j, k)
  ζ̂z(i, j, k) = x_ξ_y_η(i, j, k) − x_η_y_ξ(i, j, k)

  inverse_metrics = (;
    ξ̂x=ξ̂x, ξ̂y=ξ̂y, ξ̂z=ξ̂z, η̂x=η̂x, η̂y=η̂y, η̂z=η̂z, ζ̂x=ζ̂x, ζ̂y=ζ̂y, ζ̂z=ζ̂z
  )

  return MetricCache(forward_metrics, inverse_metrics)
end

function node_coordinates(x, y, z, iterators)
  ni, nj, nk = size(iterators.node.full)
  nhalo = iterators.nhalo
  xnodes = zeros(ni, nj, nk)
  ynodes = zeros(ni, nj, nk)
  znodes = zeros(ni, nj, nk)

  @batch for I in iterators.node.full
    i, j, k = I.I
    xnodes[I] = x(i - nhalo, j - nhalo, k - nhalo)
    ynodes[I] = y(i - nhalo, j - nhalo, k - nhalo)
    znodes[I] = z(i - nhalo, j - nhalo, k - nhalo)
  end

  return xnodes, ynodes, znodes
end

function centroid_coordinates(x, y, z, iterators)
  ni, nj, nk = size(iterators.node.full)
  nhalo = iterators.nhalo
  xnodes = zeros(ni, nj, nk)
  ynodes = zeros(ni, nj, nk)
  znodes = zeros(ni, nj, nk)

  @batch for I in iterators.node.full
    i, j, k = I.I
    xnodes[I] = x(i - nhalo + 0.5, j - nhalo + 0.5, k - nhalo + 0.5)
    ynodes[I] = y(i - nhalo + 0.5, j - nhalo + 0.5, k - nhalo + 0.5)
    znodes[I] = z(i - nhalo + 0.5, j - nhalo + 0.5, k - nhalo + 0.5)
  end

  return xnodes, ynodes, znodes
end

function cell_center_derivative(ϕ, backend)
  ϕᵢ₊½, ϕⱼ₊½, ϕₖ₊½ = edge_functions(ϕ, backend)

  ∂ϕ_∂ξ(i, j, k) = ϕᵢ₊½(i, j, k) - ϕᵢ₊½(i - 1, j, k)
  ∂ϕ_∂η(i, j, k) = ϕⱼ₊½(i, j, k) - ϕⱼ₊½(i, j - 1, k)
  ∂ϕ_∂ζ(i, j, k) = ϕₖ₊½(i, j, k) - ϕₖ₊½(i, j, k - 1)

  return (; ∂ϕ_∂ξ, ∂ϕ_∂η, ∂ϕ_∂ζ)
end

function edge_functions(ϕ, backend)

  # returns val, ∂, ∂²
  ξ_derivs(i, j, k) = value_derivative_and_second_derivative(ξ -> ϕ(ξ, j, k), backend, i)
  η_derivs(i, j, k) = value_derivative_and_second_derivative(η -> ϕ(i, η, k), backend, j)
  ζ_derivs(i, j, k) = value_derivative_and_second_derivative(ζ -> ϕ(i, j, ζ), backend, k)

  function ϕᵢ₊½(i, j, k)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(i, j, k)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(i + 1, j, k)

    ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
    ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

    return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
  end

  function ϕⱼ₊½(i, j, k)
    ϕⱼ, ∂ϕ_∂ξⱼ, ∂²ϕ_∂ξ²ⱼ = η_derivs(i, j, k)
    ϕⱼ₊₁, ∂ϕ_∂ξⱼ₊₁, ∂²ϕ_∂ξ²ⱼ₊₁ = η_derivs(i, j + 1, k)

    ϕᴸⱼ₊½ = ϕⱼ + (1 / 2) * ∂ϕ_∂ξⱼ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ
    ϕᴿⱼ₊½ = ϕⱼ₊₁ - (1 / 2) * ∂ϕ_∂ξⱼ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ₊₁

    return (ϕᴸⱼ₊½ + ϕᴿⱼ₊½) / 2
  end

  function ϕₖ₊½(i, j, k)
    ϕₖ, ∂ϕ_∂ξₖ, ∂²ϕ_∂ξ²ₖ = ζ_derivs(i, j, k)
    ϕₖ₊₁, ∂ϕ_∂ξₖ₊₁, ∂²ϕ_∂ξ²ₖ₊₁ = ζ_derivs(i, j, k + 1)

    ϕᴸₖ₊½ = ϕₖ + (1 / 2) * ∂ϕ_∂ξₖ + (1 / 12) * ∂²ϕ_∂ξ²ₖ
    ϕᴿₖ₊½ = ϕₖ₊₁ - (1 / 2) * ∂ϕ_∂ξₖ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ₖ₊₁

    return (ϕᴸₖ₊½ + ϕᴿₖ₊½) / 2
  end

  return (; ϕᵢ₊½, ϕⱼ₊½, ϕₖ₊½)
end

function get_iterators(celldims::NTuple{N,Int}, nhalo::Int) where {N}
  cellCI = CartesianIndices(celldims .+ 2nhalo)
  nodeCI = CartesianIndices(celldims .+ 1 .+ 2nhalo)

  node = (full=nodeCI, domain=expand(nodeCI, -nhalo))
  cell = (full=cellCI, domain=expand(cellCI, -nhalo))
  return (; node, cell, nhalo)
end

function compute_cell_metrics!(mesh)
  nhalo = mesh.iterators.nhalo

  for I in mesh.iterators.cell.full
    i, j, k = I.I

    # account for halo cells and centroid offset
    ξηζ = I.I .- nhalo .+ 0.5 # centroid

    # @unpack xξ, yξ, zξ, xη, yη, zη, xζ, yζ, zζ, J = forward(metric_functions_cache, ξηζ)

    # Jacobian
    J = det(mesh.metric_functions_cache.forward.jacobian(ξηζ...))
    mesh.metrics.cell.J[i, j, k] = J
    mesh.metrics.cell.x₁.ξ[i, j, k] = mesh.metric_functions_cache.forward.xξ(ξηζ...)
    mesh.metrics.cell.x₁.η[i, j, k] = mesh.metric_functions_cache.forward.xη(ξηζ...)
    mesh.metrics.cell.x₁.ζ[i, j, k] = mesh.metric_functions_cache.forward.xζ(ξηζ...)
    mesh.metrics.cell.x₂.ξ[i, j, k] = mesh.metric_functions_cache.forward.yξ(ξηζ...)
    mesh.metrics.cell.x₂.η[i, j, k] = mesh.metric_functions_cache.forward.yη(ξηζ...)
    mesh.metrics.cell.x₂.ζ[i, j, k] = mesh.metric_functions_cache.forward.yζ(ξηζ...)
    mesh.metrics.cell.x₃.ξ[i, j, k] = mesh.metric_functions_cache.forward.zξ(ξηζ...)
    mesh.metrics.cell.x₃.η[i, j, k] = mesh.metric_functions_cache.forward.zη(ξηζ...)
    mesh.metrics.cell.x₃.ζ[i, j, k] = mesh.metric_functions_cache.forward.zζ(ξηζ...)

    # hatted inverse metrics
    ξ̂_x = mesh.metric_functions_cache.inverse.ξ̂x(ξηζ...)
    ξ̂_y = mesh.metric_functions_cache.inverse.ξ̂y(ξηζ...)
    ξ̂_z = mesh.metric_functions_cache.inverse.ξ̂z(ξηζ...)
    η̂_x = mesh.metric_functions_cache.inverse.η̂x(ξηζ...)
    η̂_y = mesh.metric_functions_cache.inverse.η̂y(ξηζ...)
    η̂_z = mesh.metric_functions_cache.inverse.η̂z(ξηζ...)
    ζ̂_x = mesh.metric_functions_cache.inverse.ζ̂x(ξηζ...)
    ζ̂_y = mesh.metric_functions_cache.inverse.ζ̂y(ξηζ...)
    ζ̂_z = mesh.metric_functions_cache.inverse.ζ̂z(ξηζ...)

    mesh.metrics.cell.ξ̂.x₁[i, j, k] = ξ̂_x
    mesh.metrics.cell.ξ̂.x₂[i, j, k] = ξ̂_y
    mesh.metrics.cell.ξ̂.x₃[i, j, k] = ξ̂_z
    mesh.metrics.cell.η̂.x₁[i, j, k] = η̂_x
    mesh.metrics.cell.η̂.x₂[i, j, k] = η̂_y
    mesh.metrics.cell.η̂.x₃[i, j, k] = η̂_z
    mesh.metrics.cell.ζ̂.x₁[i, j, k] = ζ̂_x
    mesh.metrics.cell.ζ̂.x₂[i, j, k] = ζ̂_y
    mesh.metrics.cell.ζ̂.x₃[i, j, k] = ζ̂_z

    mesh.metrics.cell.ξ.x₁[i, j, k] = ξ̂_x * J
    mesh.metrics.cell.ξ.x₂[i, j, k] = ξ̂_y * J
    mesh.metrics.cell.ξ.x₃[i, j, k] = ξ̂_z * J
    mesh.metrics.cell.η.x₁[i, j, k] = η̂_x * J
    mesh.metrics.cell.η.x₂[i, j, k] = η̂_y * J
    mesh.metrics.cell.η.x₃[i, j, k] = η̂_z * J
    mesh.metrics.cell.ζ.x₁[i, j, k] = ζ̂_x * J
    mesh.metrics.cell.ζ.x₂[i, j, k] = ζ̂_y * J
    mesh.metrics.cell.ζ.x₃[i, j, k] = ζ̂_z * J
  end
end

function compute_edge_metrics!(mesh)
  nhalo = mesh.iterators.nhalo

  ξ̂xᵢ₊½, ξ̂xⱼ₊½, ξ̂xₖ₊½ = edge_functions(
    mesh.metric_functions_cache.inverse.ξ̂x, mesh.diff_backend
  )
  η̂xᵢ₊½, η̂xⱼ₊½, η̂xₖ₊½ = edge_functions(
    mesh.metric_functions_cache.inverse.η̂x, mesh.diff_backend
  )
  ζ̂xᵢ₊½, ζ̂xⱼ₊½, ζ̂xₖ₊½ = edge_functions(
    mesh.metric_functions_cache.inverse.ζ̂x, mesh.diff_backend
  )
  ξ̂yᵢ₊½, ξ̂yⱼ₊½, ξ̂yₖ₊½ = edge_functions(
    mesh.metric_functions_cache.inverse.ξ̂y, mesh.diff_backend
  )
  η̂yᵢ₊½, η̂yⱼ₊½, η̂yₖ₊½ = edge_functions(
    mesh.metric_functions_cache.inverse.η̂y, mesh.diff_backend
  )
  ζ̂yᵢ₊½, ζ̂yⱼ₊½, ζ̂yₖ₊½ = edge_functions(
    mesh.metric_functions_cache.inverse.ζ̂y, mesh.diff_backend
  )
  ξ̂zᵢ₊½, ξ̂zⱼ₊½, ξ̂zₖ₊½ = edge_functions(
    mesh.metric_functions_cache.inverse.ξ̂z, mesh.diff_backend
  )
  η̂zᵢ₊½, η̂zⱼ₊½, η̂zₖ₊½ = edge_functions(
    mesh.metric_functions_cache.inverse.η̂z, mesh.diff_backend
  )
  ζ̂zᵢ₊½, ζ̂zⱼ₊½, ζ̂zₖ₊½ = edge_functions(
    mesh.metric_functions_cache.inverse.ζ̂z, mesh.diff_backend
  )

  # for (iedge, edge) in enumerate(mesh.metrics.edge_metrics)
  i₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 1, +1)
  j₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 2, +1)
  k₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 3, +1)

  @threads for I in i₊½_edge_domain
    i, j, k = I.I
    ξηζ = I.I .- nhalo .+ (1 / 2)

    mesh.metrics.edge.i₊½.ξ̂.x₁[i, j, k] = ξ̂xᵢ₊½(ξηζ...)
    mesh.metrics.edge.i₊½.ξ̂.x₂[i, j, k] = ξ̂yᵢ₊½(ξηζ...)
    mesh.metrics.edge.i₊½.ξ̂.x₃[i, j, k] = ξ̂zᵢ₊½(ξηζ...)
    mesh.metrics.edge.i₊½.η̂.x₁[i, j, k] = η̂xᵢ₊½(ξηζ...)
    mesh.metrics.edge.i₊½.η̂.x₂[i, j, k] = η̂yᵢ₊½(ξηζ...)
    mesh.metrics.edge.i₊½.η̂.x₃[i, j, k] = η̂zᵢ₊½(ξηζ...)
    mesh.metrics.edge.i₊½.ζ̂.x₁[i, j, k] = ζ̂xᵢ₊½(ξηζ...)
    mesh.metrics.edge.i₊½.ζ̂.x₂[i, j, k] = ζ̂yᵢ₊½(ξηζ...)
    mesh.metrics.edge.i₊½.ζ̂.x₃[i, j, k] = ζ̂zᵢ₊½(ξηζ...)
  end

  @threads for I in j₊½_edge_domain
    i, j, k = I.I
    ξηζ = I.I .- nhalo .+ (1 / 2)

    mesh.metrics.edge.j₊½.ξ̂.x₁[i, j, k] = ξ̂xⱼ₊½(ξηζ...)
    mesh.metrics.edge.j₊½.ξ̂.x₂[i, j, k] = ξ̂yⱼ₊½(ξηζ...)
    mesh.metrics.edge.j₊½.ξ̂.x₃[i, j, k] = ξ̂zⱼ₊½(ξηζ...)
    mesh.metrics.edge.j₊½.η̂.x₁[i, j, k] = η̂xⱼ₊½(ξηζ...)
    mesh.metrics.edge.j₊½.η̂.x₂[i, j, k] = η̂yⱼ₊½(ξηζ...)
    mesh.metrics.edge.j₊½.η̂.x₃[i, j, k] = η̂zⱼ₊½(ξηζ...)
    mesh.metrics.edge.j₊½.ζ̂.x₁[i, j, k] = ζ̂xⱼ₊½(ξηζ...)
    mesh.metrics.edge.j₊½.ζ̂.x₂[i, j, k] = ζ̂yⱼ₊½(ξηζ...)
    mesh.metrics.edge.j₊½.ζ̂.x₃[i, j, k] = ζ̂zⱼ₊½(ξηζ...)
  end

  @threads for I in k₊½_edge_domain
    i, j, k = I.I
    ξηζ = I.I .- nhalo .+ (1 / 2)

    mesh.metrics.edge.k₊½.ξ̂.x₁[i, j, k] = ξ̂xₖ₊½(ξηζ...)
    mesh.metrics.edge.k₊½.ξ̂.x₂[i, j, k] = ξ̂yₖ₊½(ξηζ...)
    mesh.metrics.edge.k₊½.ξ̂.x₃[i, j, k] = ξ̂zₖ₊½(ξηζ...)
    mesh.metrics.edge.k₊½.η̂.x₁[i, j, k] = η̂xₖ₊½(ξηζ...)
    mesh.metrics.edge.k₊½.η̂.x₂[i, j, k] = η̂yₖ₊½(ξηζ...)
    mesh.metrics.edge.k₊½.η̂.x₃[i, j, k] = η̂zₖ₊½(ξηζ...)
    mesh.metrics.edge.k₊½.ζ̂.x₁[i, j, k] = ζ̂xₖ₊½(ξηζ...)
    mesh.metrics.edge.k₊½.ζ̂.x₂[i, j, k] = ζ̂yₖ₊½(ξηζ...)
    mesh.metrics.edge.k₊½.ζ̂.x₃[i, j, k] = ζ̂zₖ₊½(ξηζ...)
  end

  return nothing
end

function cell_center_forward_metrics(mesh, (i, j, k))

  # instead of grabbing all the individual forward metrics, just 
  # call the jacobian and take advantage of the vector-mode AD

  ξηζ = (i, j, k) .- mesh.iterators.nhalo .+ (1 / 2) # centroid
  J⃗ = mesh.metric_functions_cache.forward.jacobian(ξηζ...)

  #   J⃗ = [
  #     xξ xη xζ
  #     yξ yη yζ
  #     zξ zη zζ
  #   ]

  return (;
    xξ=J⃗[1], xη=J⃗[4], xζ=J⃗[7], yξ=J⃗[2], yη=J⃗[5], yζ=J⃗[8], zξ=J⃗[3], zη=J⃗[6], zζ=J⃗[9]
  )
end

function cell_center_jacobian_matrix(mesh, (i, j, k))
  ξηζ = (i, j, k) .- mesh.iterators.nhalo .+ (1 / 2) # centroid
  return mesh.metric_functions_cache.forward.jacobian(ξηζ...)
end

function cell_center_conservative_metrics(mesh, (i, j, k))

  # instead of grabbing all the individual forward metrics, just 
  # call the jacobian and take advantage of the vector-mode AD

  ξηζ = (i, j, k) .- mesh.iterators.nhalo .+ (1 / 2) # centroid

  return (;
    ξ̂x=mesh.metric_functions_cache.inverse.ξ̂x(ξηζ...),
    ξ̂y=mesh.metric_functions_cache.inverse.ξ̂y(ξηζ...),
    ξ̂z=mesh.metric_functions_cache.inverse.ξ̂z(ξηζ...),
    η̂x=mesh.metric_functions_cache.inverse.η̂x(ξηζ...),
    η̂y=mesh.metric_functions_cache.inverse.η̂y(ξηζ...),
    η̂z=mesh.metric_functions_cache.inverse.η̂z(ξηζ...),
    ζ̂x=mesh.metric_functions_cache.inverse.ζ̂x(ξηζ...),
    ζ̂y=mesh.metric_functions_cache.inverse.ζ̂y(ξηζ...),
    ζ̂z=mesh.metric_functions_cache.inverse.ζ̂z(ξηζ...),
  )
end