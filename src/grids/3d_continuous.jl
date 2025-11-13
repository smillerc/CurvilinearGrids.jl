
abstract type AbstractContinuousCurvilinearGrid3D <: AbstractCurvilinearGrid3D end

struct ContinuousCurvilinearGrid3D{A,B,C,EM,CM,E,BE,DBE,DS} <:
       AbstractContinuousCurvilinearGrid3D
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

function ContinuousCurvilinearGrid3D(
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

  mesh = ContinuousCurvilinearGrid3D(
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

function node_coordinates(x, y, z, iterators, backend, T)
  ni, nj, nk = size(iterators.node.full)
  nhalo = iterators.nhalo

  coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, (ni, nj, nk)),
    y=KernelAbstractions.zeros(backend, T, (ni, nj, nk)),
    z=KernelAbstractions.zeros(backend, T, (ni, nj, nk)),
  ))

  @batch for I in iterators.node.full
    i, j, k = I.I
    coords.x[I] = x(i - nhalo, j - nhalo, k - nhalo)
    coords.y[I] = y(i - nhalo, j - nhalo, k - nhalo)
    coords.z[I] = z(i - nhalo, j - nhalo, k - nhalo)
  end

  return coords
end

function centroid_coordinates(x, y, z, iterators, backend, T)
  ni, nj, nk = size(iterators.cell.full)
  nhalo = iterators.nhalo
  coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, (ni, nj, nk)),
    y=KernelAbstractions.zeros(backend, T, (ni, nj, nk)),
    z=KernelAbstractions.zeros(backend, T, (ni, nj, nk)),
  ))

  @batch for I in iterators.cell.full
    i, j, k = I.I
    coords.x[I] = x(i - nhalo + 0.5, j - nhalo + 0.5, k - nhalo + 0.5)
    coords.y[I] = y(i - nhalo + 0.5, j - nhalo + 0.5, k - nhalo + 0.5)
    coords.z[I] = z(i - nhalo + 0.5, j - nhalo + 0.5, k - nhalo + 0.5)
  end

  return coords
end

@inline function tol_diff(a::T, b::T) where {T}
  a_m_b = a - b
  return a_m_b * !isapprox(a, b; rtol=sqrt(eps(T)), atol=eps(T))
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
    i, j, k = I.I

    # account for halo cells and centroid offset
    ξηζ = I.I .- nhalo .+ 0.5 # centroid

    # @unpack xξ, yξ, zξ, xη, yη, zη, xζ, yζ, zζ, J = forward(metric_functions_cache, ξηζ)

    # Jacobian
    J = det(mesh.metric_functions_cache.forward.jacobian(ξηζ...))
    mesh.cell_center_metrics.J[i, j, k] = J

    xξ = mesh.metric_functions_cache.forward.xξ(ξηζ...)
    xη = mesh.metric_functions_cache.forward.xη(ξηζ...)
    xζ = mesh.metric_functions_cache.forward.xζ(ξηζ...)
    yξ = mesh.metric_functions_cache.forward.yξ(ξηζ...)
    yη = mesh.metric_functions_cache.forward.yη(ξηζ...)
    yζ = mesh.metric_functions_cache.forward.yζ(ξηζ...)
    zξ = mesh.metric_functions_cache.forward.zξ(ξηζ...)
    zη = mesh.metric_functions_cache.forward.zη(ξηζ...)
    zζ = mesh.metric_functions_cache.forward.zζ(ξηζ...)

    mesh.cell_center_metrics.x₁.ξ[i, j, k] = xξ * (abs(xξ) >= eps())
    mesh.cell_center_metrics.x₁.η[i, j, k] = xη * (abs(xη) >= eps())
    mesh.cell_center_metrics.x₁.ζ[i, j, k] = xζ * (abs(xζ) >= eps())
    mesh.cell_center_metrics.x₂.ξ[i, j, k] = yξ * (abs(yξ) >= eps())
    mesh.cell_center_metrics.x₂.η[i, j, k] = yη * (abs(yη) >= eps())
    mesh.cell_center_metrics.x₂.ζ[i, j, k] = yζ * (abs(yζ) >= eps())
    mesh.cell_center_metrics.x₃.ξ[i, j, k] = zξ * (abs(zξ) >= eps())
    mesh.cell_center_metrics.x₃.η[i, j, k] = zη * (abs(zη) >= eps())
    mesh.cell_center_metrics.x₃.ζ[i, j, k] = zζ * (abs(zζ) >= eps())

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

    ξ̂_x = ξ̂_x # * (abs(ξ̂_x) >= eps())
    ξ̂_y = ξ̂_y # * (abs(ξ̂_y) >= eps())
    ξ̂_z = ξ̂_z # * (abs(ξ̂_z) >= eps())
    η̂_x = η̂_x # * (abs(η̂_x) >= eps())
    η̂_y = η̂_y # * (abs(η̂_y) >= eps())
    η̂_z = η̂_z # * (abs(η̂_z) >= eps())
    ζ̂_x = ζ̂_x # * (abs(ζ̂_x) >= eps())
    ζ̂_y = ζ̂_y # * (abs(ζ̂_y) >= eps())
    ζ̂_z = ζ̂_z # * (abs(ζ̂_z) >= eps())

    mesh.cell_center_metrics.ξ̂.x₁[i, j, k] = ξ̂_x
    mesh.cell_center_metrics.ξ̂.x₂[i, j, k] = ξ̂_y
    mesh.cell_center_metrics.ξ̂.x₃[i, j, k] = ξ̂_z
    mesh.cell_center_metrics.η̂.x₁[i, j, k] = η̂_x
    mesh.cell_center_metrics.η̂.x₂[i, j, k] = η̂_y
    mesh.cell_center_metrics.η̂.x₃[i, j, k] = η̂_z
    mesh.cell_center_metrics.ζ̂.x₁[i, j, k] = ζ̂_x
    mesh.cell_center_metrics.ζ̂.x₂[i, j, k] = ζ̂_y
    mesh.cell_center_metrics.ζ̂.x₃[i, j, k] = ζ̂_z

    ξx, ηx, ζx, ξy, ηy, ζy, ξz, ηz, ζz = mesh.metric_functions_cache.inverse.Jinv(ξηζ...)
    mesh.cell_center_metrics.ξ.x₁[i, j, k] = ξx
    mesh.cell_center_metrics.ξ.x₂[i, j, k] = ξy
    mesh.cell_center_metrics.ξ.x₃[i, j, k] = ξz
    mesh.cell_center_metrics.η.x₁[i, j, k] = ηx
    mesh.cell_center_metrics.η.x₂[i, j, k] = ηy
    mesh.cell_center_metrics.η.x₃[i, j, k] = ηz
    mesh.cell_center_metrics.ζ.x₁[i, j, k] = ζx
    mesh.cell_center_metrics.ζ.x₂[i, j, k] = ζy
    mesh.cell_center_metrics.ζ.x₃[i, j, k] = ζz
  end

  if any(mesh.cell_center_metrics.J .<= 0)
    @error "Invalid Jacobians detected (i.e. negative or zero volumes)"
  end
end

# function compute_edge_metrics!(mesh)
#   nhalo = mesh.iterators.nhalo

#   ξ̂xᵢ₊½, ξ̂xⱼ₊½, ξ̂xₖ₊½ = edge_functions_3d(
#     mesh.metric_functions_cache.inverse.ξ̂x, mesh.diff_backend
#   )
#   η̂xᵢ₊½, η̂xⱼ₊½, η̂xₖ₊½ = edge_functions_3d(
#     mesh.metric_functions_cache.inverse.η̂x, mesh.diff_backend
#   )
#   ζ̂xᵢ₊½, ζ̂xⱼ₊½, ζ̂xₖ₊½ = edge_functions_3d(
#     mesh.metric_functions_cache.inverse.ζ̂x, mesh.diff_backend
#   )
#   ξ̂yᵢ₊½, ξ̂yⱼ₊½, ξ̂yₖ₊½ = edge_functions_3d(
#     mesh.metric_functions_cache.inverse.ξ̂y, mesh.diff_backend
#   )
#   η̂yᵢ₊½, η̂yⱼ₊½, η̂yₖ₊½ = edge_functions_3d(
#     mesh.metric_functions_cache.inverse.η̂y, mesh.diff_backend
#   )
#   ζ̂yᵢ₊½, ζ̂yⱼ₊½, ζ̂yₖ₊½ = edge_functions_3d(
#     mesh.metric_functions_cache.inverse.ζ̂y, mesh.diff_backend
#   )
#   ξ̂zᵢ₊½, ξ̂zⱼ₊½, ξ̂zₖ₊½ = edge_functions_3d(
#     mesh.metric_functions_cache.inverse.ξ̂z, mesh.diff_backend
#   )
#   η̂zᵢ₊½, η̂zⱼ₊½, η̂zₖ₊½ = edge_functions_3d(
#     mesh.metric_functions_cache.inverse.η̂z, mesh.diff_backend
#   )
#   ζ̂zᵢ₊½, ζ̂zⱼ₊½, ζ̂zₖ₊½ = edge_functions_3d(
#     mesh.metric_functions_cache.inverse.ζ̂z, mesh.diff_backend
#   )

#   Jinv_ᵢ₊½, Jinv_ⱼ₊½, Jinv_ₖ₊½ = edge_functions_3d(
#     mesh.metric_functions_cache.inverse.Jinv, mesh.diff_backend
#   )

#   normJinv_ᵢ₊½, normJinv_ⱼ₊½, normJinv_ₖ₊½ = edge_functions_3d(
#     mesh.metric_functions_cache.inverse.Jinv_norm, mesh.diff_backend
#   )

#   Jᵢ₊½, Jⱼ₊½, Jₖ₊½ = edge_functions_3d(
#     mesh.metric_functions_cache.forward.J, mesh.diff_backend
#   )

#   # for (iedge, edge) in enumerate(mesh.metrics.edge_metrics)
#   i₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 1, +1)
#   j₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 2, +1)
#   k₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 3, +1)

#   @threads for I in mesh.iterators.cell.full
#     i, j, k = I.I
#     ξηζ = I.I .- nhalo .+ (1 / 2) # centroid index

#     J = Jᵢ₊½(ξηζ...)

#     ξ̂_xᵢ₊½, η̂_xᵢ₊½, ζ̂_xᵢ₊½, ξ̂_yᵢ₊½, η̂_yᵢ₊½, ζ̂_yᵢ₊½, ξ̂_zᵢ₊½, η̂_zᵢ₊½, ζ̂_zᵢ₊½ = normJinv_ᵢ₊½(
#       ξηζ...
#     )

#     mesh.edge_metrics.i₊½.J[i, j, k] = J
#     mesh.edge_metrics.i₊½.ξ̂.x₁[i, j, k] = ξ̂_xᵢ₊½
#     mesh.edge_metrics.i₊½.ξ̂.x₂[i, j, k] = ξ̂_yᵢ₊½
#     mesh.edge_metrics.i₊½.ξ̂.x₃[i, j, k] = ξ̂_zᵢ₊½
#     mesh.edge_metrics.i₊½.η̂.x₁[i, j, k] = η̂_xᵢ₊½
#     mesh.edge_metrics.i₊½.η̂.x₂[i, j, k] = η̂_yᵢ₊½
#     mesh.edge_metrics.i₊½.η̂.x₃[i, j, k] = η̂_zᵢ₊½
#     mesh.edge_metrics.i₊½.ζ̂.x₁[i, j, k] = ζ̂_xᵢ₊½
#     mesh.edge_metrics.i₊½.ζ̂.x₂[i, j, k] = ζ̂_yᵢ₊½
#     mesh.edge_metrics.i₊½.ζ̂.x₃[i, j, k] = ζ̂_zᵢ₊½

#     ξ_xᵢ₊½, η_xᵢ₊½, ζ_xᵢ₊½, ξ_yᵢ₊½, η_yᵢ₊½, ζ_yᵢ₊½, ξ_zᵢ₊½, η_zᵢ₊½, ζ_zᵢ₊½ = Jinv_ᵢ₊½(
#       ξηζ...
#     )

#     mesh.edge_metrics.i₊½.ξ.x₁[i, j, k] = ξ_xᵢ₊½
#     mesh.edge_metrics.i₊½.ξ.x₂[i, j, k] = ξ_yᵢ₊½
#     mesh.edge_metrics.i₊½.ξ.x₃[i, j, k] = ξ_zᵢ₊½
#     mesh.edge_metrics.i₊½.η.x₁[i, j, k] = η_xᵢ₊½
#     mesh.edge_metrics.i₊½.η.x₂[i, j, k] = η_yᵢ₊½
#     mesh.edge_metrics.i₊½.η.x₃[i, j, k] = η_zᵢ₊½
#     mesh.edge_metrics.i₊½.ζ.x₁[i, j, k] = ζ_xᵢ₊½
#     mesh.edge_metrics.i₊½.ζ.x₂[i, j, k] = ζ_yᵢ₊½
#     mesh.edge_metrics.i₊½.ζ.x₃[i, j, k] = ζ_zᵢ₊½
#   end

#   @threads for I in mesh.iterators.cell.full
#     i, j, k = I.I
#     ξηζ = I.I .- nhalo .+ (1 / 2) # centroid index

#     J = Jⱼ₊½(ξηζ...)

#     ξ̂_xⱼ₊½, η̂_xⱼ₊½, ζ̂_xⱼ₊½, ξ̂_yⱼ₊½, η̂_yⱼ₊½, ζ̂_yⱼ₊½, ξ̂_zⱼ₊½, η̂_zⱼ₊½, ζ̂_zⱼ₊½ = normJinv_ⱼ₊½(
#       ξηζ...
#     )

#     mesh.edge_metrics.j₊½.J[i, j, k] = J
#     mesh.edge_metrics.j₊½.ξ̂.x₁[i, j, k] = ξ̂_xⱼ₊½
#     mesh.edge_metrics.j₊½.ξ̂.x₂[i, j, k] = ξ̂_yⱼ₊½
#     mesh.edge_metrics.j₊½.ξ̂.x₃[i, j, k] = ξ̂_zⱼ₊½
#     mesh.edge_metrics.j₊½.η̂.x₁[i, j, k] = η̂_xⱼ₊½
#     mesh.edge_metrics.j₊½.η̂.x₂[i, j, k] = η̂_yⱼ₊½
#     mesh.edge_metrics.j₊½.η̂.x₃[i, j, k] = η̂_zⱼ₊½
#     mesh.edge_metrics.j₊½.ζ̂.x₁[i, j, k] = ζ̂_xⱼ₊½
#     mesh.edge_metrics.j₊½.ζ̂.x₂[i, j, k] = ζ̂_yⱼ₊½
#     mesh.edge_metrics.j₊½.ζ̂.x₃[i, j, k] = ζ̂_zⱼ₊½

#     ξ_xⱼ₊½, η_xⱼ₊½, ζ_xⱼ₊½, ξ_yⱼ₊½, η_yⱼ₊½, ζ_yⱼ₊½, ξ_zⱼ₊½, η_zⱼ₊½, ζ_zⱼ₊½ = Jinv_ⱼ₊½(
#       ξηζ...
#     )

#     mesh.edge_metrics.j₊½.ξ.x₁[i, j, k] = ξ_xⱼ₊½
#     mesh.edge_metrics.j₊½.ξ.x₂[i, j, k] = ξ_yⱼ₊½
#     mesh.edge_metrics.j₊½.ξ.x₃[i, j, k] = ξ_zⱼ₊½
#     mesh.edge_metrics.j₊½.η.x₁[i, j, k] = η_xⱼ₊½
#     mesh.edge_metrics.j₊½.η.x₂[i, j, k] = η_yⱼ₊½
#     mesh.edge_metrics.j₊½.η.x₃[i, j, k] = η_zⱼ₊½
#     mesh.edge_metrics.j₊½.ζ.x₁[i, j, k] = ζ_xⱼ₊½
#     mesh.edge_metrics.j₊½.ζ.x₂[i, j, k] = ζ_yⱼ₊½
#     mesh.edge_metrics.j₊½.ζ.x₃[i, j, k] = ζ_zⱼ₊½
#   end

#   @threads for I in mesh.iterators.cell.full
#     i, j, k = I.I
#     ξηζ = I.I .- nhalo .+ (1 / 2) # centroid index

#     J = Jₖ₊½(ξηζ...)

#     ξ̂_xₖ₊½, η̂_xₖ₊½, ζ̂_xₖ₊½, ξ̂_yₖ₊½, η̂_yₖ₊½, ζ̂_yₖ₊½, ξ̂_zₖ₊½, η̂_zₖ₊½, ζ̂_zₖ₊½ = normJinv_ₖ₊½(
#       ξηζ...
#     )

#     mesh.edge_metrics.k₊½.J[i, j, k] = J
#     mesh.edge_metrics.k₊½.ξ̂.x₁[i, j, k] = ξ̂_xₖ₊½
#     mesh.edge_metrics.k₊½.ξ̂.x₂[i, j, k] = ξ̂_yₖ₊½
#     mesh.edge_metrics.k₊½.ξ̂.x₃[i, j, k] = ξ̂_zₖ₊½
#     mesh.edge_metrics.k₊½.η̂.x₁[i, j, k] = η̂_xₖ₊½
#     mesh.edge_metrics.k₊½.η̂.x₂[i, j, k] = η̂_yₖ₊½
#     mesh.edge_metrics.k₊½.η̂.x₃[i, j, k] = η̂_zₖ₊½
#     mesh.edge_metrics.k₊½.ζ̂.x₁[i, j, k] = ζ̂_xₖ₊½
#     mesh.edge_metrics.k₊½.ζ̂.x₂[i, j, k] = ζ̂_yₖ₊½
#     mesh.edge_metrics.k₊½.ζ̂.x₃[i, j, k] = ζ̂_zₖ₊½

#     ξ_xₖ₊½, η_xₖ₊½, ζ_xₖ₊½, ξ_yₖ₊½, η_yₖ₊½, ζ_yₖ₊½, ξ_zₖ₊½, η_zₖ₊½, ζ_zₖ₊½ = Jinv_ₖ₊½(
#       ξηζ...
#     )

#     mesh.edge_metrics.k₊½.ξ.x₁[i, j, k] = ξ_xₖ₊½
#     mesh.edge_metrics.k₊½.ξ.x₂[i, j, k] = ξ_yₖ₊½
#     mesh.edge_metrics.k₊½.ξ.x₃[i, j, k] = ξ_zₖ₊½
#     mesh.edge_metrics.k₊½.η.x₁[i, j, k] = η_xₖ₊½
#     mesh.edge_metrics.k₊½.η.x₂[i, j, k] = η_yₖ₊½
#     mesh.edge_metrics.k₊½.η.x₃[i, j, k] = η_zₖ₊½
#     mesh.edge_metrics.k₊½.ζ.x₁[i, j, k] = ζ_xₖ₊½
#     mesh.edge_metrics.k₊½.ζ.x₂[i, j, k] = ζ_yₖ₊½
#     mesh.edge_metrics.k₊½.ζ.x₃[i, j, k] = ζ_zₖ₊½
#   end

#   return nothing
# end

function compute_edge_metrics!(mesh)
  nhalo = mesh.iterators.nhalo

  ξ̂xᵢ₊½, ξ̂xⱼ₊½, ξ̂xₖ₊½ = edge_functions_3d(
    mesh.metric_functions_cache.inverse.ξ̂x, mesh.diff_backend
  )
  η̂xᵢ₊½, η̂xⱼ₊½, η̂xₖ₊½ = edge_functions_3d(
    mesh.metric_functions_cache.inverse.η̂x, mesh.diff_backend
  )
  ζ̂xᵢ₊½, ζ̂xⱼ₊½, ζ̂xₖ₊½ = edge_functions_3d(
    mesh.metric_functions_cache.inverse.ζ̂x, mesh.diff_backend
  )
  ξ̂yᵢ₊½, ξ̂yⱼ₊½, ξ̂yₖ₊½ = edge_functions_3d(
    mesh.metric_functions_cache.inverse.ξ̂y, mesh.diff_backend
  )
  η̂yᵢ₊½, η̂yⱼ₊½, η̂yₖ₊½ = edge_functions_3d(
    mesh.metric_functions_cache.inverse.η̂y, mesh.diff_backend
  )
  ζ̂yᵢ₊½, ζ̂yⱼ₊½, ζ̂yₖ₊½ = edge_functions_3d(
    mesh.metric_functions_cache.inverse.ζ̂y, mesh.diff_backend
  )
  ξ̂zᵢ₊½, ξ̂zⱼ₊½, ξ̂zₖ₊½ = edge_functions_3d(
    mesh.metric_functions_cache.inverse.ξ̂z, mesh.diff_backend
  )
  η̂zᵢ₊½, η̂zⱼ₊½, η̂zₖ₊½ = edge_functions_3d(
    mesh.metric_functions_cache.inverse.η̂z, mesh.diff_backend
  )
  ζ̂zᵢ₊½, ζ̂zⱼ₊½, ζ̂zₖ₊½ = edge_functions_3d(
    mesh.metric_functions_cache.inverse.ζ̂z, mesh.diff_backend
  )

  Jinv_ᵢ₊½, Jinv_ⱼ₊½, Jinv_ₖ₊½ = edge_functions_3d(
    mesh.metric_functions_cache.inverse.Jinv, mesh.diff_backend
  )

  Jᵢ₊½, Jⱼ₊½, Jₖ₊½ = edge_functions_3d(
    mesh.metric_functions_cache.forward.J, mesh.diff_backend
  )

  # for (iedge, edge) in enumerate(mesh.metrics.edge_metrics)
  i₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 1, +1)
  j₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 2, +1)
  k₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 3, +1)

  @threads for I in i₊½_edge_domain
    i, j, k = I.I
    ξηζ = I.I .- nhalo .+ (1 / 2) # centroid index

    J = Jᵢ₊½(ξηζ...)
    ξ̂_xᵢ₊½ = ξ̂xᵢ₊½(ξηζ...)
    ξ̂_yᵢ₊½ = ξ̂yᵢ₊½(ξηζ...)
    ξ̂_zᵢ₊½ = ξ̂zᵢ₊½(ξηζ...)
    η̂_xᵢ₊½ = η̂xᵢ₊½(ξηζ...)
    η̂_yᵢ₊½ = η̂yᵢ₊½(ξηζ...)
    η̂_zᵢ₊½ = η̂zᵢ₊½(ξηζ...)
    ζ̂_xᵢ₊½ = ζ̂xᵢ₊½(ξηζ...)
    ζ̂_yᵢ₊½ = ζ̂yᵢ₊½(ξηζ...)
    ζ̂_zᵢ₊½ = ζ̂zᵢ₊½(ξηζ...)

    mesh.edge_metrics.i₊½.J[i, j, k] = J
    mesh.edge_metrics.i₊½.ξ̂.x₁[i, j, k] = ξ̂_xᵢ₊½ * (abs(ξ̂_xᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.ξ̂.x₂[i, j, k] = ξ̂_yᵢ₊½ * (abs(ξ̂_yᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.ξ̂.x₃[i, j, k] = ξ̂_zᵢ₊½ * (abs(ξ̂_zᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.η̂.x₁[i, j, k] = η̂_xᵢ₊½ * (abs(η̂_xᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.η̂.x₂[i, j, k] = η̂_yᵢ₊½ * (abs(η̂_yᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.η̂.x₃[i, j, k] = η̂_zᵢ₊½ * (abs(η̂_zᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.ζ̂.x₁[i, j, k] = ζ̂_xᵢ₊½ * (abs(ζ̂_xᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.ζ̂.x₂[i, j, k] = ζ̂_yᵢ₊½ * (abs(ζ̂_yᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.ζ̂.x₃[i, j, k] = ζ̂_zᵢ₊½ * (abs(ζ̂_zᵢ₊½) >= eps())

    ξ_xᵢ₊½, η_xᵢ₊½, ζ_xᵢ₊½, ξ_yᵢ₊½, η_yᵢ₊½, ζ_yᵢ₊½, ξ_zᵢ₊½, η_zᵢ₊½, ζ_zᵢ₊½ = Jinv_ᵢ₊½(
      ξηζ...
    )

    mesh.edge_metrics.i₊½.ξ.x₁[i, j, k] = ifelse(
      isfinite(ξ_xᵢ₊½), ξ_xᵢ₊½ * (abs(ξ_xᵢ₊½) >= eps()), zero(ξ_xᵢ₊½)
    )
    mesh.edge_metrics.i₊½.ξ.x₂[i, j, k] = ifelse(
      isfinite(ξ_yᵢ₊½), ξ_yᵢ₊½ * (abs(ξ_yᵢ₊½) >= eps()), zero(ξ_yᵢ₊½)
    )
    mesh.edge_metrics.i₊½.ξ.x₃[i, j, k] = ifelse(
      isfinite(ξ_zᵢ₊½), ξ_zᵢ₊½ * (abs(ξ_zᵢ₊½) >= eps()), zero(ξ_zᵢ₊½)
    )
    mesh.edge_metrics.i₊½.η.x₁[i, j, k] = ifelse(
      isfinite(η_xᵢ₊½), η_xᵢ₊½ * (abs(η_xᵢ₊½) >= eps()), zero(η_xᵢ₊½)
    )
    mesh.edge_metrics.i₊½.η.x₂[i, j, k] = ifelse(
      isfinite(η_yᵢ₊½), η_yᵢ₊½ * (abs(η_yᵢ₊½) >= eps()), zero(η_yᵢ₊½)
    )
    mesh.edge_metrics.i₊½.η.x₃[i, j, k] = ifelse(
      isfinite(η_zᵢ₊½), η_zᵢ₊½ * (abs(η_zᵢ₊½) >= eps()), zero(η_zᵢ₊½)
    )
    mesh.edge_metrics.i₊½.ζ.x₁[i, j, k] = ifelse(
      isfinite(ζ_xᵢ₊½), ζ_xᵢ₊½ * (abs(ζ_xᵢ₊½) >= eps()), zero(ζ_xᵢ₊½)
    )
    mesh.edge_metrics.i₊½.ζ.x₂[i, j, k] = ifelse(
      isfinite(ζ_yᵢ₊½), ζ_yᵢ₊½ * (abs(ζ_yᵢ₊½) >= eps()), zero(ζ_yᵢ₊½)
    )
    mesh.edge_metrics.i₊½.ζ.x₃[i, j, k] = ifelse(
      isfinite(ζ_zᵢ₊½), ζ_zᵢ₊½ * (abs(ζ_zᵢ₊½) >= eps()), zero(ζ_zᵢ₊½)
    )
  end

  @threads for I in j₊½_edge_domain
    i, j, k = I.I
    ξηζ = I.I .- nhalo .+ (1 / 2) # centroid index

    J = Jⱼ₊½(ξηζ...)
    ξ̂_xⱼ₊½ = ξ̂xⱼ₊½(ξηζ...)
    ξ̂_yⱼ₊½ = ξ̂yⱼ₊½(ξηζ...)
    ξ̂_zⱼ₊½ = ξ̂zⱼ₊½(ξηζ...)
    η̂_xⱼ₊½ = η̂xⱼ₊½(ξηζ...)
    η̂_yⱼ₊½ = η̂yⱼ₊½(ξηζ...)
    η̂_zⱼ₊½ = η̂zⱼ₊½(ξηζ...)
    ζ̂_xⱼ₊½ = ζ̂xⱼ₊½(ξηζ...)
    ζ̂_yⱼ₊½ = ζ̂yⱼ₊½(ξηζ...)
    ζ̂_zⱼ₊½ = ζ̂zⱼ₊½(ξηζ...)

    mesh.edge_metrics.j₊½.J[i, j, k] = J
    mesh.edge_metrics.j₊½.ξ̂.x₁[i, j, k] = ξ̂_xⱼ₊½ * (abs(ξ̂_xⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.ξ̂.x₂[i, j, k] = ξ̂_yⱼ₊½ * (abs(ξ̂_yⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.ξ̂.x₃[i, j, k] = ξ̂_zⱼ₊½ * (abs(ξ̂_zⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.η̂.x₁[i, j, k] = η̂_xⱼ₊½ * (abs(η̂_xⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.η̂.x₂[i, j, k] = η̂_yⱼ₊½ * (abs(η̂_yⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.η̂.x₃[i, j, k] = η̂_zⱼ₊½ * (abs(η̂_zⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.ζ̂.x₁[i, j, k] = ζ̂_xⱼ₊½ * (abs(ζ̂_xⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.ζ̂.x₂[i, j, k] = ζ̂_yⱼ₊½ * (abs(ζ̂_yⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.ζ̂.x₃[i, j, k] = ζ̂_zⱼ₊½ * (abs(ζ̂_zⱼ₊½) >= eps())

    ξ_xⱼ₊½, η_xⱼ₊½, ζ_xⱼ₊½, ξ_yⱼ₊½, η_yⱼ₊½, ζ_yⱼ₊½, ξ_zⱼ₊½, η_zⱼ₊½, ζ_zⱼ₊½ = Jinv_ⱼ₊½(
      ξηζ...
    )

    mesh.edge_metrics.j₊½.ξ.x₁[i, j, k] = ifelse(
      isfinite(ξ_xⱼ₊½), ξ_xⱼ₊½ * (abs(ξ_xⱼ₊½) >= eps()), zero(ξ_xⱼ₊½)
    )
    mesh.edge_metrics.j₊½.ξ.x₂[i, j, k] = ifelse(
      isfinite(ξ_yⱼ₊½), ξ_yⱼ₊½ * (abs(ξ_yⱼ₊½) >= eps()), zero(ξ_yⱼ₊½)
    )
    mesh.edge_metrics.j₊½.ξ.x₃[i, j, k] = ifelse(
      isfinite(ξ_zⱼ₊½), ξ_zⱼ₊½ * (abs(ξ_zⱼ₊½) >= eps()), zero(ξ_zⱼ₊½)
    )
    mesh.edge_metrics.j₊½.η.x₁[i, j, k] = ifelse(
      isfinite(η_xⱼ₊½), η_xⱼ₊½ * (abs(η_xⱼ₊½) >= eps()), zero(η_xⱼ₊½)
    )
    mesh.edge_metrics.j₊½.η.x₂[i, j, k] = ifelse(
      isfinite(η_yⱼ₊½), η_yⱼ₊½ * (abs(η_yⱼ₊½) >= eps()), zero(η_yⱼ₊½)
    )
    mesh.edge_metrics.j₊½.η.x₃[i, j, k] = ifelse(
      isfinite(η_zⱼ₊½), η_zⱼ₊½ * (abs(η_zⱼ₊½) >= eps()), zero(η_zⱼ₊½)
    )
    mesh.edge_metrics.j₊½.ζ.x₁[i, j, k] = ifelse(
      isfinite(ζ_xⱼ₊½), ζ_xⱼ₊½ * (abs(ζ_xⱼ₊½) >= eps()), zero(ζ_xⱼ₊½)
    )
    mesh.edge_metrics.j₊½.ζ.x₂[i, j, k] = ifelse(
      isfinite(ζ_yⱼ₊½), ζ_yⱼ₊½ * (abs(ζ_yⱼ₊½) >= eps()), zero(ζ_yⱼ₊½)
    )
    mesh.edge_metrics.j₊½.ζ.x₃[i, j, k] = ifelse(
      isfinite(ζ_zⱼ₊½), ζ_zⱼ₊½ * (abs(ζ_zⱼ₊½) >= eps()), zero(ζ_zⱼ₊½)
    )
  end

  @threads for I in k₊½_edge_domain
    i, j, k = I.I
    ξηζ = I.I .- nhalo .+ (1 / 2) # centroid index

    J = Jₖ₊½(ξηζ...)
    ξ̂_xₖ₊½ = ξ̂xₖ₊½(ξηζ...)
    ξ̂_yₖ₊½ = ξ̂yₖ₊½(ξηζ...)
    ξ̂_zₖ₊½ = ξ̂zₖ₊½(ξηζ...)
    η̂_xₖ₊½ = η̂xₖ₊½(ξηζ...)
    η̂_yₖ₊½ = η̂yₖ₊½(ξηζ...)
    η̂_zₖ₊½ = η̂zₖ₊½(ξηζ...)
    ζ̂_xₖ₊½ = ζ̂xₖ₊½(ξηζ...)
    ζ̂_yₖ₊½ = ζ̂yₖ₊½(ξηζ...)
    ζ̂_zₖ₊½ = ζ̂zₖ₊½(ξηζ...)

    mesh.edge_metrics.k₊½.J[i, j, k] = J
    mesh.edge_metrics.k₊½.ξ̂.x₁[i, j, k] = ξ̂_xₖ₊½ * (abs(ξ̂_xₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.ξ̂.x₂[i, j, k] = ξ̂_yₖ₊½ * (abs(ξ̂_yₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.ξ̂.x₃[i, j, k] = ξ̂_zₖ₊½ * (abs(ξ̂_zₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.η̂.x₁[i, j, k] = η̂_xₖ₊½ * (abs(η̂_xₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.η̂.x₂[i, j, k] = η̂_yₖ₊½ * (abs(η̂_yₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.η̂.x₃[i, j, k] = η̂_zₖ₊½ * (abs(η̂_zₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.ζ̂.x₁[i, j, k] = ζ̂_xₖ₊½ * (abs(ζ̂_xₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.ζ̂.x₂[i, j, k] = ζ̂_yₖ₊½ * (abs(ζ̂_yₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.ζ̂.x₃[i, j, k] = ζ̂_zₖ₊½ * (abs(ζ̂_zₖ₊½) >= eps())

    ξ_xₖ₊½, η_xₖ₊½, ζ_xₖ₊½, ξ_yₖ₊½, η_yₖ₊½, ζ_yₖ₊½, ξ_zₖ₊½, η_zₖ₊½, ζ_zₖ₊½ = Jinv_ₖ₊½(
      ξηζ...
    )

    mesh.edge_metrics.k₊½.ξ.x₁[i, j, k] = ifelse(
      isfinite(ξ_xₖ₊½), ξ_xₖ₊½ * (abs(ξ_xₖ₊½) >= eps()), zero(ξ_xₖ₊½)
    )
    mesh.edge_metrics.k₊½.ξ.x₂[i, j, k] = ifelse(
      isfinite(ξ_yₖ₊½), ξ_yₖ₊½ * (abs(ξ_yₖ₊½) >= eps()), zero(ξ_yₖ₊½)
    )
    mesh.edge_metrics.k₊½.ξ.x₃[i, j, k] = ifelse(
      isfinite(ξ_zₖ₊½), ξ_zₖ₊½ * (abs(ξ_zₖ₊½) >= eps()), zero(ξ_zₖ₊½)
    )
    mesh.edge_metrics.k₊½.η.x₁[i, j, k] = ifelse(
      isfinite(η_xₖ₊½), η_xₖ₊½ * (abs(η_xₖ₊½) >= eps()), zero(η_xₖ₊½)
    )
    mesh.edge_metrics.k₊½.η.x₂[i, j, k] = ifelse(
      isfinite(η_yₖ₊½), η_yₖ₊½ * (abs(η_yₖ₊½) >= eps()), zero(η_yₖ₊½)
    )
    mesh.edge_metrics.k₊½.η.x₃[i, j, k] = ifelse(
      isfinite(η_zₖ₊½), η_zₖ₊½ * (abs(η_zₖ₊½) >= eps()), zero(η_zₖ₊½)
    )
    mesh.edge_metrics.k₊½.ζ.x₁[i, j, k] = ifelse(
      isfinite(ζ_xₖ₊½), ζ_xₖ₊½ * (abs(ζ_xₖ₊½) >= eps()), zero(ζ_xₖ₊½)
    )
    mesh.edge_metrics.k₊½.ζ.x₂[i, j, k] = ifelse(
      isfinite(ζ_yₖ₊½), ζ_yₖ₊½ * (abs(ζ_yₖ₊½) >= eps()), zero(ζ_yₖ₊½)
    )
    mesh.edge_metrics.k₊½.ζ.x₃[i, j, k] = ifelse(
      isfinite(ζ_zₖ₊½), ζ_zₖ₊½ * (abs(ζ_zₖ₊½) >= eps()), zero(ζ_zₖ₊½)
    )
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

function metric_derivative_functions(mesh)

  #! format: off
  ∂Jinv_∂ξ(i,j,k) = derivative(ξ -> mesh.metric_functions_cache.inverse.Jinv(ξ, j, k), mesh.diff_backend, i)
  ∂Jinv_∂η(i,j,k) = derivative(η -> mesh.metric_functions_cache.inverse.Jinv(i, η, k), mesh.diff_backend, j)
  ∂Jinv_∂ζ(i,j,k) = derivative(ζ -> mesh.metric_functions_cache.inverse.Jinv(i, j, ζ), mesh.diff_backend, k)
  #! format: on

  #! format: off
  # ∂ξx_∂ξ(i,j,k) = derivative(ξ -> mesh.metric_functions_cache.forward.ξx(ξ, j, k), mesh.diff_backend, i)
  # ∂ξy_∂ξ(i,j,k) = derivative(ξ -> mesh.metric_functions_cache.forward.ξy(ξ, j, k), mesh.diff_backend, i)
  # ∂ξz_∂ξ(i,j,k) = derivative(ξ -> mesh.metric_functions_cache.forward.ξz(ξ, j, k), mesh.diff_backend, i)
  # ∂ξx_∂η(i,j,k) = derivative(η -> mesh.metric_functions_cache.forward.ξx(i, η, k), mesh.diff_backend, j)
  # ∂ξy_∂η(i,j,k) = derivative(η -> mesh.metric_functions_cache.forward.ξy(i, η, k), mesh.diff_backend, j)
  # ∂ξz_∂η(i,j,k) = derivative(η -> mesh.metric_functions_cache.forward.ξz(i, η, k), mesh.diff_backend, j)
  # ∂ξx_∂ζ(i,j,k) = derivative(ζ -> mesh.metric_functions_cache.forward.ξx(i, j, ζ), mesh.diff_backend, i)
  # ∂ξy_∂ζ(i,j,k) = derivative(ζ -> mesh.metric_functions_cache.forward.ξy(i, j, ζ), mesh.diff_backend, k)
  # ∂ξz_∂ζ(i,j,k) = derivative(ζ -> mesh.metric_functions_cache.forward.ξz(i, j, ζ), mesh.diff_backend, k)
  # ∂ηx_∂ξ(i,j,k) = derivative(ξ -> mesh.metric_functions_cache.forward.ηx(ξ, j, k), mesh.diff_backend, k)
  # ∂ηy_∂ξ(i,j,k) = derivative(ξ -> mesh.metric_functions_cache.forward.ηy(ξ, j, k), mesh.diff_backend, i)
  # ∂ηz_∂ξ(i,j,k) = derivative(ξ -> mesh.metric_functions_cache.forward.ηz(ξ, j, k), mesh.diff_backend, i)
  # ∂ηx_∂η(i,j,k) = derivative(η -> mesh.metric_functions_cache.forward.ηx(i, η, k), mesh.diff_backend, j)
  # ∂ηy_∂η(i,j,k) = derivative(η -> mesh.metric_functions_cache.forward.ηy(i, η, k), mesh.diff_backend, j)
  # ∂ηz_∂η(i,j,k) = derivative(η -> mesh.metric_functions_cache.forward.ηz(i, η, k), mesh.diff_backend, j)
  # ∂ηx_∂ζ(i,j,k) = derivative(ζ -> mesh.metric_functions_cache.forward.ηx(i, j, ζ), mesh.diff_backend, k)
  # ∂ηy_∂ζ(i,j,k) = derivative(ζ -> mesh.metric_functions_cache.forward.ηy(i, j, ζ), mesh.diff_backend, k)
  # ∂ηz_∂ζ(i,j,k) = derivative(ζ -> mesh.metric_functions_cache.forward.ηz(i, j, ζ), mesh.diff_backend, k)
  # ∂ζx_∂ξ(i,j,k) = derivative(ξ -> mesh.metric_functions_cache.forward.ζx(ξ, j, k), mesh.diff_backend, i)
  # ∂ζy_∂ξ(i,j,k) = derivative(ξ -> mesh.metric_functions_cache.forward.ζy(ξ, j, k), mesh.diff_backend, i)
  # ∂ζz_∂ξ(i,j,k) = derivative(ξ -> mesh.metric_functions_cache.forward.ζz(ξ, j, k), mesh.diff_backend, i)
  # ∂ζx_∂η(i,j,k) = derivative(η -> mesh.metric_functions_cache.forward.ζx(i, η, k), mesh.diff_backend, j)
  # ∂ζy_∂η(i,j,k) = derivative(η -> mesh.metric_functions_cache.forward.ζy(i, η, k), mesh.diff_backend, j)
  # ∂ζz_∂η(i,j,k) = derivative(η -> mesh.metric_functions_cache.forward.ζz(i, η, k), mesh.diff_backend, j)
  # ∂ζx_∂ζ(i,j,k) = derivative(ζ -> mesh.metric_functions_cache.forward.ζx(i, j, ζ), mesh.diff_backend, k)
  # ∂ζy_∂ζ(i,j,k) = derivative(ζ -> mesh.metric_functions_cache.forward.ζy(i, j, ζ), mesh.diff_backend, k)
  # ∂ζz_∂ζ(i,j,k) = derivative(ζ -> mesh.metric_functions_cache.forward.ζz(i, j, ζ), mesh.diff_backend, k)
  #! format: on

  return (; ∂Jinv_∂ξ, ∂Jinv_∂η, ∂Jinv_∂ζ)
end