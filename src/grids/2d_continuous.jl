
using ForwardDiff, StructArraysernelAbstractions, CartesianDomains, UnPack
using DifferentiationInterface

using StaticArrays, LinearAlgebra
using WriteVTK, Polyester, .Threads

abstract type AbstractContinuousCurvilinearGrid2D <: AbstractCurvilinearGrid2D end

struct ContinuousCurvilinearGrid2D{A,B,C,EM,CM,E,BE,DBE,DS} <:
       AbstractContinuousCurvilinearGrid2D
  node_coordinates::A
  centroid_coordinates::A
  mapping_functions::B
  metric_functions_cache::C
  edge_metrics::EM
  cell_center_metrics::CM
  iterators::E
  backend::BE
  diff_backend::DBE
  nhalo::Int
  discretization_scheme::DS
  discretization_scheme_name::Symbol
end

struct MetricCache{F1,F2}
  forward::F1
  inverse::F2
end

function ContinuousCurvilinearGrid2D(
  x::Function,
  y::Function,
  z::Function,
  celldims::NTuple,
  discretization_scheme::Symbol,
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  T=Float64;
  compute_metrics=true,
  global_node_indices::Union{Nothing,CartesianIndices{3}}=nothing,
)
  GradientDiscretizationScheme, order, _, nhalo, scheme_name = get_gradient_discretization_scheme(
    discretization_scheme
  )
  iterators = get_iterators(celldims, nhalo)

  cell_center_metrics, edge_metrics = get_metric_soa(celldims .+ 2nhalo, backend, T)
  xyz_n = node_coordinates(x, y, z, iterators, backend, T)
  xyz_c = centroid_coordinates(x, y, z, iterators, backend, T)

  discr_scheme = GradientDiscretizationScheme(order; use_cache=false)

  mesh = ContinuousCurvilinearGrid2D(
    xyz_n,
    xyz_c,
    (; x, y, z),
    MetricCache(x, y, z, diff_backend),
    edge_metrics,
    cell_center_metrics,
    iterators,
    backend,
    diff_backend,
    nhalo,
    discr_scheme,
    scheme_name,
  )

  if compute_metrics
    compute_cell_metrics!(mesh)
    compute_edge_metrics!(mesh)
  end

  return mesh
end

function MetricCache2D(x::Function, y::Function, backend)
  xξ(i, j) = derivative(ξ -> x(ξ, j), backend, i)
  xη(i, j) = derivative(η -> x(i, η), backend, j)
  yξ(i, j) = derivative(ξ -> y(ξ, j), backend, i)
  yη(i, j) = derivative(η -> y(i, η), backend, j)

  function jacobian_matrix(i, j)

    # compute the jacobian matrix w/o any extra logic
    function compute_jacobian_matrix(i, j)
      jac = DifferentiationInterface.jacobian(
        u -> SVector(x(u...), y(u...)), backend, @SVector [i, j]
      )

      jac = @. jac * (abs(jac) >= eps())
      return jac
    end

    jac_matrix = compute_jacobian_matrix(i, j)
    J = det(jac_matrix)

    # sometimes the Jacobian can be zero for various reasons 
    # (derivs are zero)

    if !iszero(J) && isfinite(J)
      return jac_matrix
    else

      # perturb the coordinates if we're at a singularity

      pert = sqrt(eps())
      iter = 1
      itermax = 25

      last_valid_jacobian_matrix = jac_matrix

      # perturb until we get zero again, and then use the last nonzero value
      while true
        if iter > itermax
          error(
            "Maximum iteration count reached ($itermax) when trying to iteratively calculate the Jacobian near a singularity",
          )
        end

        # perturb the coordinates
        ξηζ = (i, j) .+ pert

        # compute the jacobian
        _jac = compute_jacobian_matrix(ξηζ...)
        _J = det(_jac)

        if !isfinite(_J) || iszero(_J)
          # we've made the perturbation too small
          break
        else
          last_valid_jacobian_matrix = _jac
        end

        # make the perturbation smaller
        pert = pert / 10
        iter += 1
      end

      last_valid_jacobian_matrix = @. last_valid_jacobian_matrix *
        (abs(last_valid_jacobian_matrix) >= eps())
      return last_valid_jacobian_matrix
    end
  end

  jacobian(i, j) = det(jacobian_matrix(i, j))
  function jinv(i, j)
    J_inv = inv(jacobian_matrix(i, j))

    J_inv = @. J_inv * (abs(J_inv) >= eps())
    return J_inv
  end

  forward_metrics = (; jacobian=jacobian_matrix, J=jacobian, xξ=xξ, xη=xη, yξ=yξ, yη=yη)

  # x_ξ_y(i, j) = xξ(i, j) * y(i, j)
  # x_η_y(i, j) = xη(i, j) * y(i, j)
  # x_ζ_y(i, j) = xζ(i, j) * y(i, j)

  # y_ξ_z(i, j) = yξ(i, j) * z(i, j)
  # y_η_z(i, j) = yη(i, j) * z(i, j)
  # y_ζ_z(i, j) = yζ(i, j) * z(i, j)

  # z_ξ_x(i, j) = zξ(i, j) * x(i, j)
  # z_η_x(i, j) = zη(i, j) * x(i, j)
  # z_ζ_x(i, j) = zζ(i, j) * x(i, j)

  # _, x_ξ_y_η, x_ξ_y_ζ = cell_center_derivative(x_ξ_y, backend)
  # x_η_y_ξ, _, x_η_y_ζ = cell_center_derivative(x_η_y, backend)
  # x_ζ_y_ξ, x_ζ_y_η, _ = cell_center_derivative(x_ζ_y, backend)

  # _, y_ξ_z_η, y_ξ_z_ζ = cell_center_derivative(y_ξ_z, backend)
  # y_η_z_ξ, _, y_η_z_ζ = cell_center_derivative(y_η_z, backend)
  # y_ζ_z_ξ, y_ζ_z_η, _ = cell_center_derivative(y_ζ_z, backend)

  # _, z_ξ_x_η, z_ξ_x_ζ = cell_center_derivative(z_ξ_x, backend)
  # z_η_x_ξ, _, z_η_x_ζ = cell_center_derivative(z_η_x, backend)
  # z_ζ_x_ξ, z_ζ_x_η, _ = cell_center_derivative(z_ζ_x, backend)

  # # Do NOT put eps() tolerance checks on these! It will create GCL-related errors
  # ξ̂x(i, j) = y_η_z_ζ(i, j) − y_ζ_z_η(i, j)
  # η̂x(i, j) = y_ζ_z_ξ(i, j) − y_ξ_z_ζ(i, j)
  # ζ̂x(i, j) = y_ξ_z_η(i, j) − y_η_z_ξ(i, j)
  # ξ̂y(i, j) = z_η_x_ζ(i, j) − z_ζ_x_η(i, j)
  # η̂y(i, j) = z_ζ_x_ξ(i, j) − z_ξ_x_ζ(i, j)
  # ζ̂y(i, j) = z_ξ_x_η(i, j) − z_η_x_ξ(i, j)
  # ξ̂z(i, j) = x_η_y_ζ(i, j) − x_ζ_y_η(i, j)
  # η̂z(i, j) = x_ζ_y_ξ(i, j) − x_ξ_y_ζ(i, j)
  # ζ̂z(i, j) = x_ξ_y_η(i, j) − x_η_y_ξ(i, j)

  inverse_metrics = (;
    # ξ̂x=ξ̂x,
    # ξ̂y=ξ̂y,
    # η̂x=η̂x,
    # η̂y=η̂y,
    Jinv=jinv,
  )

  return MetricCache(forward_metrics, inverse_metrics)
end

function node_coordinates(x, y, iterators, backend, T)
  ni, nj = size(iterators.node.full)
  nhalo = iterators.nhalo

  coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, (ni, nj)),
    y=KernelAbstractions.zeros(backend, T, (ni, nj)),
  ))

  @batch for I in iterators.node.full
    i, j = I.I
    coords.x[I] = x(i - nhalo, j - nhalo)
    coords.y[I] = y(i - nhalo, j - nhalo)
  end

  return coords
end

function centroid_coordinates(x, y, iterators, backend, T)
  ni, nj = size(iterators.cell.full)
  nhalo = iterators.nhalo
  coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, (ni, nj, nk)),
    y=KernelAbstractions.zeros(backend, T, (ni, nj, nk)),
  ))

  @batch for I in iterators.cell.full
    i, j = I.I
    coords.x[I] = x(i - nhalo + 0.5, j - nhalo + 0.5)
    coords.y[I] = y(i - nhalo + 0.5, j - nhalo + 0.5)
  end

  return coords
end

@inline function tol_diff(a::T, b::T) where {T}
  a_m_b = a - b
  return a_m_b * !isapprox(a, b; rtol=sqrt(eps(T)), atol=eps(T))
end

function cell_center_derivative(ϕ, backend)
  ϕᵢ₊½, ϕⱼ₊½ = edge_functions_2d(ϕ, backend)

  ∂ϕ_∂ξ(i, j) = ϕᵢ₊½(i, j) - ϕᵢ₊½(i - 1, j)
  ∂ϕ_∂η(i, j) = ϕⱼ₊½(i, j) - ϕⱼ₊½(i, j - 1)

  return (; ∂ϕ_∂ξ, ∂ϕ_∂η)
end

function edge_functions_2d(ϕ, backend)

  # returns val, ∂, ∂²
  ξ_derivs(i, j) = value_derivative_and_second_derivative(ξ -> ϕ(ξ, j), backend, i)
  η_derivs(i, j) = value_derivative_and_second_derivative(η -> ϕ(i, η), backend, j)

  function ϕᵢ₊½(i, j)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(i, j)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(i + 1, j)

    ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
    ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

    return (ϕᴸᵢ₊½ .+ ϕᴿᵢ₊½) ./ 2
  end

  function ϕⱼ₊½(i, j)
    ϕⱼ, ∂ϕ_∂ξⱼ, ∂²ϕ_∂ξ²ⱼ = η_derivs(i, j)
    ϕⱼ₊₁, ∂ϕ_∂ξⱼ₊₁, ∂²ϕ_∂ξ²ⱼ₊₁ = η_derivs(i, j + 1)

    ϕᴸⱼ₊½ = ϕⱼ + (1 / 2) * ∂ϕ_∂ξⱼ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ
    ϕᴿⱼ₊½ = ϕⱼ₊₁ - (1 / 2) * ∂ϕ_∂ξⱼ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ₊₁

    return (ϕᴸⱼ₊½ + ϕᴿⱼ₊½) / 2
  end

  return (; ϕᵢ₊½, ϕⱼ₊½)
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

  @threads for I in mesh.iterators.cell.full
    i, j = I.I

    # account for halo cells and centroid offset
    ξηζ = I.I .- nhalo .+ 0.5 # centroid

    # @unpack xξ, yξ, zξ, xη, yη, zη, xζ, yζ, zζ, J = forward(metric_functions_cache, ξηζ)

    # Jacobian
    J = det(mesh.metric_functions_cache.forward.jacobian(ξηζ...))
    # if iszero(J)
    #   ξηζ = ξηζ .+ sqrt(eps())
    #   J = det(mesh.metric_functions_cache.forward.jacobian(ξηζ...))
    # end

    mesh.cell_center_metrics.J[i, j] = J
    mesh.cell_center_metrics.x₁.ξ[i, j] = mesh.metric_functions_cache.forward.xξ(ξηζ...)
    mesh.cell_center_metrics.x₁.η[i, j] = mesh.metric_functions_cache.forward.xη(ξηζ...)
    mesh.cell_center_metrics.x₁.ζ[i, j] = mesh.metric_functions_cache.forward.xζ(ξηζ...)
    mesh.cell_center_metrics.x₂.ξ[i, j] = mesh.metric_functions_cache.forward.yξ(ξηζ...)
    mesh.cell_center_metrics.x₂.η[i, j] = mesh.metric_functions_cache.forward.yη(ξηζ...)
    mesh.cell_center_metrics.x₂.ζ[i, j] = mesh.metric_functions_cache.forward.yζ(ξηζ...)
    mesh.cell_center_metrics.x₃.ξ[i, j] = mesh.metric_functions_cache.forward.zξ(ξηζ...)
    mesh.cell_center_metrics.x₃.η[i, j] = mesh.metric_functions_cache.forward.zη(ξηζ...)
    mesh.cell_center_metrics.x₃.ζ[i, j] = mesh.metric_functions_cache.forward.zζ(ξηζ...)

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

    mesh.cell_center_metrics.ξ̂.x₁[i, j] = ξ̂_x
    mesh.cell_center_metrics.ξ̂.x₂[i, j] = ξ̂_y
    mesh.cell_center_metrics.ξ̂.x₃[i, j] = ξ̂_z
    mesh.cell_center_metrics.η̂.x₁[i, j] = η̂_x
    mesh.cell_center_metrics.η̂.x₂[i, j] = η̂_y
    mesh.cell_center_metrics.η̂.x₃[i, j] = η̂_z
    mesh.cell_center_metrics.ζ̂.x₁[i, j] = ζ̂_x
    mesh.cell_center_metrics.ζ̂.x₂[i, j] = ζ̂_y
    mesh.cell_center_metrics.ζ̂.x₃[i, j] = ζ̂_z

    mesh.cell_center_metrics.ξ.x₁[i, j] = ξ̂_x / J
    mesh.cell_center_metrics.ξ.x₂[i, j] = ξ̂_y / J
    mesh.cell_center_metrics.ξ.x₃[i, j] = ξ̂_z / J
    mesh.cell_center_metrics.η.x₁[i, j] = η̂_x / J
    mesh.cell_center_metrics.η.x₂[i, j] = η̂_y / J
    mesh.cell_center_metrics.η.x₃[i, j] = η̂_z / J
    mesh.cell_center_metrics.ζ.x₁[i, j] = ζ̂_x / J
    mesh.cell_center_metrics.ζ.x₂[i, j] = ζ̂_y / J
    mesh.cell_center_metrics.ζ.x₃[i, j] = ζ̂_z / J
  end

  if any(mesh.cell_center_metrics.J .<= 0)
    @error "Invalid Jacobians detected (i.e. negative or zero volumes)"
  end
end

function compute_edge_metrics!(mesh)
  nhalo = mesh.iterators.nhalo

  ξ̂xᵢ₊½, ξ̂xⱼ₊½ = edge_functions_2d(
    mesh.metric_functions_cache.inverse.ξ̂x, mesh.diff_backend
  )
  η̂xᵢ₊½, η̂xⱼ₊½ = edge_functions_2d(
    mesh.metric_functions_cache.inverse.η̂x, mesh.diff_backend
  )

  ξ̂yᵢ₊½, ξ̂yⱼ₊½ = edge_functions_2d(
    mesh.metric_functions_cache.inverse.ξ̂y, mesh.diff_backend
  )
  η̂yᵢ₊½, η̂yⱼ₊½ = edge_functions_2d(
    mesh.metric_functions_cache.inverse.η̂y, mesh.diff_backend
  )

  Jinv_ᵢ₊½, Jinv_ⱼ₊½ = edge_functions_2d(
    mesh.metric_functions_cache.inverse.Jinv, mesh.diff_backend
  )

  # for (iedge, edge) in enumerate(mesh.metrics.edge_metrics)
  i₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 1, +1)
  j₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 2, +1)

  @threads for I in i₊½_edge_domain
    i, j = I.I
    ξηζ = I.I .- nhalo .+ (1 / 2) # centroid index

    # J = Jᵢ₊½(ξηζ...)
    ξ̂_xᵢ₊½ = ξ̂xᵢ₊½(ξηζ...)
    ξ̂_yᵢ₊½ = ξ̂yᵢ₊½(ξηζ...)
    η̂_xᵢ₊½ = η̂xᵢ₊½(ξηζ...)
    η̂_yᵢ₊½ = η̂yᵢ₊½(ξηζ...)

    mesh.edge_metrics.i₊½.ξ̂.x₁[i, j] = ξ̂_xᵢ₊½
    mesh.edge_metrics.i₊½.ξ̂.x₂[i, j] = ξ̂_yᵢ₊½
    mesh.edge_metrics.i₊½.η̂.x₁[i, j] = η̂_xᵢ₊½
    mesh.edge_metrics.i₊½.η̂.x₂[i, j] = η̂_yᵢ₊½

    ξ_xᵢ₊½, η_xᵢ₊½, ξ_yᵢ₊½, η_yᵢ₊½ = Jinv_ᵢ₊½(ξηζ...)

    mesh.edge_metrics.i₊½.ξ.x₁[i, j] = ifelse(isfinite(ξ_xᵢ₊½), ξ_xᵢ₊½, zero(ξ_xᵢ₊½))
    mesh.edge_metrics.i₊½.ξ.x₂[i, j] = ifelse(isfinite(ξ_yᵢ₊½), ξ_yᵢ₊½, zero(ξ_yᵢ₊½))
    mesh.edge_metrics.i₊½.η.x₁[i, j] = ifelse(isfinite(η_xᵢ₊½), η_xᵢ₊½, zero(η_xᵢ₊½))
    mesh.edge_metrics.i₊½.η.x₂[i, j] = ifelse(isfinite(η_yᵢ₊½), η_yᵢ₊½, zero(η_yᵢ₊½))
  end

  @threads for I in j₊½_edge_domain
    i, j = I.I
    ξηζ = I.I .- nhalo .+ (1 / 2) # centroid index

    # J = Jⱼ₊½(ξηζ...)
    ξ̂_xⱼ₊½ = ξ̂xⱼ₊½(ξηζ...)
    ξ̂_yⱼ₊½ = ξ̂yⱼ₊½(ξηζ...)
    η̂_xⱼ₊½ = η̂xⱼ₊½(ξηζ...)
    η̂_yⱼ₊½ = η̂yⱼ₊½(ξηζ...)

    mesh.edge_metrics.j₊½.ξ̂.x₁[i, j] = ξ̂_xⱼ₊½
    mesh.edge_metrics.j₊½.ξ̂.x₂[i, j] = ξ̂_yⱼ₊½
    mesh.edge_metrics.j₊½.η̂.x₁[i, j] = η̂_xⱼ₊½
    mesh.edge_metrics.j₊½.η̂.x₂[i, j] = η̂_yⱼ₊½

    ξ_xⱼ₊½, η_xⱼ₊½, ξ_yⱼ₊½, η_yⱼ₊½ = Jinv_ⱼ₊½(ξηζ...)

    mesh.edge_metrics.j₊½.ξ.x₁[i, j] = ifelse(isfinite(ξ_xⱼ₊½), ξ_xⱼ₊½, zero(ξ_xⱼ₊½))
    mesh.edge_metrics.j₊½.ξ.x₂[i, j] = ifelse(isfinite(ξ_yⱼ₊½), ξ_yⱼ₊½, zero(ξ_yⱼ₊½))
    mesh.edge_metrics.j₊½.η.x₁[i, j] = ifelse(isfinite(η_xⱼ₊½), η_xⱼ₊½, zero(η_xⱼ₊½))
    mesh.edge_metrics.j₊½.η.x₂[i, j] = ifelse(isfinite(η_yⱼ₊½), η_yⱼ₊½, zero(η_yⱼ₊½))
  end

  return nothing
end
