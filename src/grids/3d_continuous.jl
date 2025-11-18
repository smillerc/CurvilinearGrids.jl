
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
  mapping_function_parameters,
  celldims::NTuple,
  discretization_scheme::Symbol,
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  t=zero(Float64),
  T=Float64;
  compute_metrics=true,
  global_node_indices::Union{Nothing,CartesianIndices{3}}=nothing,
)
  GradientDiscretizationScheme, order, _, nhalo, scheme_name = get_gradient_discretization_scheme(
    discretization_scheme
  )
  iterators = get_iterators(celldims, nhalo)

  cell_center_metrics, edge_metrics = get_metric_soa(celldims .+ 2nhalo, backend, T)
  ni_nodes, nj_nodes, nk_nodes = size(iterators.node.full)
  node_coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes, nk_nodes)),
    y=KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes, nk_nodes)),
    z=KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes, nk_nodes)),
  ))

  ni_cells, nj_cells, nk_cells = size(iterators.cell.full)
  centroid_coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells, nk_cells)),
    y=KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells, nk_cells)),
    z=KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells, nk_cells)),
  ))
  discr_scheme = GradientDiscretizationScheme(order; use_cache=false)
  metric_cache = MetricCache(x, y, z, diff_backend)
  mesh = ContinuousCurvilinearGrid3D(
    node_coords,
    centroid_coords,
    (; x, y, z),
    metric_cache,
    edge_metrics,
    cell_center_metrics,
    iterators,
    backend,
    diff_backend,
    nhalo,
    discr_scheme,
    scheme_name,
  )

  compute_node_coordinates!(mesh, t, mapping_function_parameters)
  compute_centroid_coordinates!(mesh, t, mapping_function_parameters)

  if compute_metrics
    compute_cell_metrics!(mesh, t, mapping_function_parameters)
    compute_edge_metrics!(mesh, t, mapping_function_parameters)
  end

  return mesh
end

function update_mapping_functions!(
  mesh::ContinuousCurvilinearGrid3D, t, new_params, compute_metrics=true
)
  mesh.mapping_function_params = new_params
  compute_node_coordinates!(mesh, t, new_params)
  compute_centroid_coordinates!(mesh, t, new_params)

  if compute_metrics
    compute_cell_metrics!(mesh, t, new_params)
    compute_edge_metrics!(mesh, t, new_params)
  end

  return nothing
end

function compute_node_coordinates!(mesh::ContinuousCurvilinearGrid3D, t, params)
  nhalo = mesh.nhalo

  x = mesh.mapping_functions.x
  y = mesh.mapping_functions.y
  z = mesh.mapping_functions.z
  for I in mesh.iterators.node.full
    ξη = I.I .- nhalo
    mesh.node_coordinates.x[I] = x(t, ξη..., params)
    mesh.node_coordinates.y[I] = y(t, ξη..., params)
    mesh.node_coordinates.z[I] = z(t, ξη..., params)
  end

  return nothing
end

function compute_centroid_coordinates!(mesh::ContinuousCurvilinearGrid3D, t, params)
  nhalo = mesh.nhalo

  x = mesh.mapping_functions.x
  y = mesh.mapping_functions.y
  z = mesh.mapping_functions.z

  for I in mesh.iterators.cell.full
    # account for halo cells and centroid offset
    ξη = I.I .- nhalo .+ 0.5 # centroid
    mesh.centroid_coordinates.x[I] = x(t, ξη..., params)
    mesh.centroid_coordinates.y[I] = y(t, ξη..., params)
    mesh.centroid_coordinates.y[I] = z(t, ξη..., params)
  end

  return nothing
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

function compute_cell_metrics!(mesh, t, params)
  nhalo = mesh.iterators.nhalo

  @threads for I in mesh.iterators.cell.full

    # account for halo cells and centroid offset
    ξηζ = I.I .- nhalo .+ 0.5 # centroid

    # @unpack xξ, yξ, zξ, xη, yη, zη, xζ, yζ, zζ, J = forward(metric_functions_cache, ξηζ)

    # Jacobian
    J = det(mesh.metric_functions_cache.forward.jacobian(t, ξηζ..., params))
    mesh.cell_center_metrics.J[I] = J

    xξ = mesh.metric_functions_cache.forward.xξ(t, ξηζ..., params)
    xη = mesh.metric_functions_cache.forward.xη(t, ξηζ..., params)
    xζ = mesh.metric_functions_cache.forward.xζ(t, ξηζ..., params)
    yξ = mesh.metric_functions_cache.forward.yξ(t, ξηζ..., params)
    yη = mesh.metric_functions_cache.forward.yη(t, ξηζ..., params)
    yζ = mesh.metric_functions_cache.forward.yζ(t, ξηζ..., params)
    zξ = mesh.metric_functions_cache.forward.zξ(t, ξηζ..., params)
    zη = mesh.metric_functions_cache.forward.zη(t, ξηζ..., params)
    zζ = mesh.metric_functions_cache.forward.zζ(t, ξηζ..., params)

    mesh.cell_center_metrics.x₁.ξ[I] = xξ * (abs(xξ) >= eps())
    mesh.cell_center_metrics.x₁.η[I] = xη * (abs(xη) >= eps())
    mesh.cell_center_metrics.x₁.ζ[I] = xζ * (abs(xζ) >= eps())
    mesh.cell_center_metrics.x₂.ξ[I] = yξ * (abs(yξ) >= eps())
    mesh.cell_center_metrics.x₂.η[I] = yη * (abs(yη) >= eps())
    mesh.cell_center_metrics.x₂.ζ[I] = yζ * (abs(yζ) >= eps())
    mesh.cell_center_metrics.x₃.ξ[I] = zξ * (abs(zξ) >= eps())
    mesh.cell_center_metrics.x₃.η[I] = zη * (abs(zη) >= eps())
    mesh.cell_center_metrics.x₃.ζ[I] = zζ * (abs(zζ) >= eps())

    # hatted inverse metrics
    ξ̂_x = mesh.metric_functions_cache.inverse.ξ̂x(t, ξηζ..., params)
    ξ̂_y = mesh.metric_functions_cache.inverse.ξ̂y(t, ξηζ..., params)
    ξ̂_z = mesh.metric_functions_cache.inverse.ξ̂z(t, ξηζ..., params)
    η̂_x = mesh.metric_functions_cache.inverse.η̂x(t, ξηζ..., params)
    η̂_y = mesh.metric_functions_cache.inverse.η̂y(t, ξηζ..., params)
    η̂_z = mesh.metric_functions_cache.inverse.η̂z(t, ξηζ..., params)
    ζ̂_x = mesh.metric_functions_cache.inverse.ζ̂x(t, ξηζ..., params)
    ζ̂_y = mesh.metric_functions_cache.inverse.ζ̂y(t, ξηζ..., params)
    ζ̂_z = mesh.metric_functions_cache.inverse.ζ̂z(t, ξηζ..., params)

    ξ̂_x = ξ̂_x # * (abs(ξ̂_x) >= eps())
    ξ̂_y = ξ̂_y # * (abs(ξ̂_y) >= eps())
    ξ̂_z = ξ̂_z # * (abs(ξ̂_z) >= eps())
    η̂_x = η̂_x # * (abs(η̂_x) >= eps())
    η̂_y = η̂_y # * (abs(η̂_y) >= eps())
    η̂_z = η̂_z # * (abs(η̂_z) >= eps())
    ζ̂_x = ζ̂_x # * (abs(ζ̂_x) >= eps())
    ζ̂_y = ζ̂_y # * (abs(ζ̂_y) >= eps())
    ζ̂_z = ζ̂_z # * (abs(ζ̂_z) >= eps())

    mesh.cell_center_metrics.ξ̂.x₁[I] = ξ̂_x
    mesh.cell_center_metrics.ξ̂.x₂[I] = ξ̂_y
    mesh.cell_center_metrics.ξ̂.x₃[I] = ξ̂_z
    mesh.cell_center_metrics.η̂.x₁[I] = η̂_x
    mesh.cell_center_metrics.η̂.x₂[I] = η̂_y
    mesh.cell_center_metrics.η̂.x₃[I] = η̂_z
    mesh.cell_center_metrics.ζ̂.x₁[I] = ζ̂_x
    mesh.cell_center_metrics.ζ̂.x₂[I] = ζ̂_y
    mesh.cell_center_metrics.ζ̂.x₃[I] = ζ̂_z

    ξx, ηx, ζx, ξy, ηy, ζy, ξz, ηz, ζz = mesh.metric_functions_cache.inverse.Jinv(
      t, ξηζ..., params
    )
    mesh.cell_center_metrics.ξ.x₁[I] = ξx
    mesh.cell_center_metrics.ξ.x₂[I] = ξy
    mesh.cell_center_metrics.ξ.x₃[I] = ξz
    mesh.cell_center_metrics.η.x₁[I] = ηx
    mesh.cell_center_metrics.η.x₂[I] = ηy
    mesh.cell_center_metrics.η.x₃[I] = ηz
    mesh.cell_center_metrics.ζ.x₁[I] = ζx
    mesh.cell_center_metrics.ζ.x₂[I] = ζy
    mesh.cell_center_metrics.ζ.x₃[I] = ζz
  end

  if any(mesh.cell_center_metrics.J .<= 0)
    @error "Invalid Jacobians detected (i.e. negative or zero volumes)"
  end
end

function compute_edge_metrics!(mesh, t, params)
  nhalo = mesh.iterators.nhalo

  # for (iedge, edge) in enumerate(mesh.metrics.edge_metrics)
  i₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 1, +1)
  j₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 2, +1)
  k₊½_edge_domain = expand_lower(mesh.iterators.cell.domain, 3, +1)

  @threads for I in i₊½_edge_domain
    ξηζ = I.I .- nhalo .+ (1 / 2) # centroid index

    J = mesh.metric_functions_cache.edge.Jᵢ₊½(t, ξηζ..., params)
    ξ̂_xᵢ₊½ = mesh.metric_functions_cache.edge.ξ̂xᵢ₊½(t, ξηζ..., params)
    ξ̂_yᵢ₊½ = mesh.metric_functions_cache.edge.ξ̂yᵢ₊½(t, ξηζ..., params)
    ξ̂_zᵢ₊½ = mesh.metric_functions_cache.edge.ξ̂zᵢ₊½(t, ξηζ..., params)
    η̂_xᵢ₊½ = mesh.metric_functions_cache.edge.η̂xᵢ₊½(t, ξηζ..., params)
    η̂_yᵢ₊½ = mesh.metric_functions_cache.edge.η̂yᵢ₊½(t, ξηζ..., params)
    η̂_zᵢ₊½ = mesh.metric_functions_cache.edge.η̂zᵢ₊½(t, ξηζ..., params)
    ζ̂_xᵢ₊½ = mesh.metric_functions_cache.edge.ζ̂xᵢ₊½(t, ξηζ..., params)
    ζ̂_yᵢ₊½ = mesh.metric_functions_cache.edge.ζ̂yᵢ₊½(t, ξηζ..., params)
    ζ̂_zᵢ₊½ = mesh.metric_functions_cache.edge.ζ̂zᵢ₊½(t, ξηζ..., params)

    mesh.edge_metrics.i₊½.J[I] = J
    mesh.edge_metrics.i₊½.ξ̂.x₁[I] = ξ̂_xᵢ₊½ * (abs(ξ̂_xᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.ξ̂.x₂[I] = ξ̂_yᵢ₊½ * (abs(ξ̂_yᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.ξ̂.x₃[I] = ξ̂_zᵢ₊½ * (abs(ξ̂_zᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.η̂.x₁[I] = η̂_xᵢ₊½ * (abs(η̂_xᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.η̂.x₂[I] = η̂_yᵢ₊½ * (abs(η̂_yᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.η̂.x₃[I] = η̂_zᵢ₊½ * (abs(η̂_zᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.ζ̂.x₁[I] = ζ̂_xᵢ₊½ * (abs(ζ̂_xᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.ζ̂.x₂[I] = ζ̂_yᵢ₊½ * (abs(ζ̂_yᵢ₊½) >= eps())
    mesh.edge_metrics.i₊½.ζ̂.x₃[I] = ζ̂_zᵢ₊½ * (abs(ζ̂_zᵢ₊½) >= eps())

    ξ_xᵢ₊½, η_xᵢ₊½, ζ_xᵢ₊½, ξ_yᵢ₊½, η_yᵢ₊½, ζ_yᵢ₊½, ξ_zᵢ₊½, η_zᵢ₊½, ζ_zᵢ₊½ = mesh.metric_functions_cache.edge.Jinv_ᵢ₊½(
      t, ξηζ..., params
    )

    mesh.edge_metrics.i₊½.ξ.x₁[I] = ifelse(
      isfinite(ξ_xᵢ₊½), ξ_xᵢ₊½ * (abs(ξ_xᵢ₊½) >= eps()), zero(ξ_xᵢ₊½)
    )
    mesh.edge_metrics.i₊½.ξ.x₂[I] = ifelse(
      isfinite(ξ_yᵢ₊½), ξ_yᵢ₊½ * (abs(ξ_yᵢ₊½) >= eps()), zero(ξ_yᵢ₊½)
    )
    mesh.edge_metrics.i₊½.ξ.x₃[I] = ifelse(
      isfinite(ξ_zᵢ₊½), ξ_zᵢ₊½ * (abs(ξ_zᵢ₊½) >= eps()), zero(ξ_zᵢ₊½)
    )
    mesh.edge_metrics.i₊½.η.x₁[I] = ifelse(
      isfinite(η_xᵢ₊½), η_xᵢ₊½ * (abs(η_xᵢ₊½) >= eps()), zero(η_xᵢ₊½)
    )
    mesh.edge_metrics.i₊½.η.x₂[I] = ifelse(
      isfinite(η_yᵢ₊½), η_yᵢ₊½ * (abs(η_yᵢ₊½) >= eps()), zero(η_yᵢ₊½)
    )
    mesh.edge_metrics.i₊½.η.x₃[I] = ifelse(
      isfinite(η_zᵢ₊½), η_zᵢ₊½ * (abs(η_zᵢ₊½) >= eps()), zero(η_zᵢ₊½)
    )
    mesh.edge_metrics.i₊½.ζ.x₁[I] = ifelse(
      isfinite(ζ_xᵢ₊½), ζ_xᵢ₊½ * (abs(ζ_xᵢ₊½) >= eps()), zero(ζ_xᵢ₊½)
    )
    mesh.edge_metrics.i₊½.ζ.x₂[I] = ifelse(
      isfinite(ζ_yᵢ₊½), ζ_yᵢ₊½ * (abs(ζ_yᵢ₊½) >= eps()), zero(ζ_yᵢ₊½)
    )
    mesh.edge_metrics.i₊½.ζ.x₃[I] = ifelse(
      isfinite(ζ_zᵢ₊½), ζ_zᵢ₊½ * (abs(ζ_zᵢ₊½) >= eps()), zero(ζ_zᵢ₊½)
    )
  end

  @threads for I in j₊½_edge_domain
    ξηζ = I.I .- nhalo .+ (1 / 2) # centroid index

    J = mesh.metric_functions_cache.edge.Jⱼ₊½(t, ξηζ..., params)
    ξ̂_xⱼ₊½ = mesh.metric_functions_cache.edge.ξ̂xⱼ₊½(t, ξηζ..., params)
    ξ̂_yⱼ₊½ = mesh.metric_functions_cache.edge.ξ̂yⱼ₊½(t, ξηζ..., params)
    ξ̂_zⱼ₊½ = mesh.metric_functions_cache.edge.ξ̂zⱼ₊½(t, ξηζ..., params)
    η̂_xⱼ₊½ = mesh.metric_functions_cache.edge.η̂xⱼ₊½(t, ξηζ..., params)
    η̂_yⱼ₊½ = mesh.metric_functions_cache.edge.η̂yⱼ₊½(t, ξηζ..., params)
    η̂_zⱼ₊½ = mesh.metric_functions_cache.edge.η̂zⱼ₊½(t, ξηζ..., params)
    ζ̂_xⱼ₊½ = mesh.metric_functions_cache.edge.ζ̂xⱼ₊½(t, ξηζ..., params)
    ζ̂_yⱼ₊½ = mesh.metric_functions_cache.edge.ζ̂yⱼ₊½(t, ξηζ..., params)
    ζ̂_zⱼ₊½ = mesh.metric_functions_cache.edge.ζ̂zⱼ₊½(t, ξηζ..., params)

    mesh.edge_metrics.j₊½.J[I] = J
    mesh.edge_metrics.j₊½.ξ̂.x₁[I] = ξ̂_xⱼ₊½ * (abs(ξ̂_xⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.ξ̂.x₂[I] = ξ̂_yⱼ₊½ * (abs(ξ̂_yⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.ξ̂.x₃[I] = ξ̂_zⱼ₊½ * (abs(ξ̂_zⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.η̂.x₁[I] = η̂_xⱼ₊½ * (abs(η̂_xⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.η̂.x₂[I] = η̂_yⱼ₊½ * (abs(η̂_yⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.η̂.x₃[I] = η̂_zⱼ₊½ * (abs(η̂_zⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.ζ̂.x₁[I] = ζ̂_xⱼ₊½ * (abs(ζ̂_xⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.ζ̂.x₂[I] = ζ̂_yⱼ₊½ * (abs(ζ̂_yⱼ₊½) >= eps())
    mesh.edge_metrics.j₊½.ζ̂.x₃[I] = ζ̂_zⱼ₊½ * (abs(ζ̂_zⱼ₊½) >= eps())

    ξ_xⱼ₊½, η_xⱼ₊½, ζ_xⱼ₊½, ξ_yⱼ₊½, η_yⱼ₊½, ζ_yⱼ₊½, ξ_zⱼ₊½, η_zⱼ₊½, ζ_zⱼ₊½ = mesh.metric_functions_cache.edge.Jinv_ⱼ₊½(
      t, ξηζ..., params
    )

    mesh.edge_metrics.j₊½.ξ.x₁[I] = ifelse(
      isfinite(ξ_xⱼ₊½), ξ_xⱼ₊½ * (abs(ξ_xⱼ₊½) >= eps()), zero(ξ_xⱼ₊½)
    )
    mesh.edge_metrics.j₊½.ξ.x₂[I] = ifelse(
      isfinite(ξ_yⱼ₊½), ξ_yⱼ₊½ * (abs(ξ_yⱼ₊½) >= eps()), zero(ξ_yⱼ₊½)
    )
    mesh.edge_metrics.j₊½.ξ.x₃[I] = ifelse(
      isfinite(ξ_zⱼ₊½), ξ_zⱼ₊½ * (abs(ξ_zⱼ₊½) >= eps()), zero(ξ_zⱼ₊½)
    )
    mesh.edge_metrics.j₊½.η.x₁[I] = ifelse(
      isfinite(η_xⱼ₊½), η_xⱼ₊½ * (abs(η_xⱼ₊½) >= eps()), zero(η_xⱼ₊½)
    )
    mesh.edge_metrics.j₊½.η.x₂[I] = ifelse(
      isfinite(η_yⱼ₊½), η_yⱼ₊½ * (abs(η_yⱼ₊½) >= eps()), zero(η_yⱼ₊½)
    )
    mesh.edge_metrics.j₊½.η.x₃[I] = ifelse(
      isfinite(η_zⱼ₊½), η_zⱼ₊½ * (abs(η_zⱼ₊½) >= eps()), zero(η_zⱼ₊½)
    )
    mesh.edge_metrics.j₊½.ζ.x₁[I] = ifelse(
      isfinite(ζ_xⱼ₊½), ζ_xⱼ₊½ * (abs(ζ_xⱼ₊½) >= eps()), zero(ζ_xⱼ₊½)
    )
    mesh.edge_metrics.j₊½.ζ.x₂[I] = ifelse(
      isfinite(ζ_yⱼ₊½), ζ_yⱼ₊½ * (abs(ζ_yⱼ₊½) >= eps()), zero(ζ_yⱼ₊½)
    )
    mesh.edge_metrics.j₊½.ζ.x₃[I] = ifelse(
      isfinite(ζ_zⱼ₊½), ζ_zⱼ₊½ * (abs(ζ_zⱼ₊½) >= eps()), zero(ζ_zⱼ₊½)
    )
  end

  @threads for I in k₊½_edge_domain
    ξηζ = I.I .- nhalo .+ (1 / 2) # centroid index

    J = mesh.metric_functions_cache.edge.Jₖ₊½(t, ξηζ..., params)
    ξ̂_xₖ₊½ = mesh.metric_functions_cache.edge.ξ̂xₖ₊½(t, ξηζ..., params)
    ξ̂_yₖ₊½ = mesh.metric_functions_cache.edge.ξ̂yₖ₊½(t, ξηζ..., params)
    ξ̂_zₖ₊½ = mesh.metric_functions_cache.edge.ξ̂zₖ₊½(t, ξηζ..., params)
    η̂_xₖ₊½ = mesh.metric_functions_cache.edge.η̂xₖ₊½(t, ξηζ..., params)
    η̂_yₖ₊½ = mesh.metric_functions_cache.edge.η̂yₖ₊½(t, ξηζ..., params)
    η̂_zₖ₊½ = mesh.metric_functions_cache.edge.η̂zₖ₊½(t, ξηζ..., params)
    ζ̂_xₖ₊½ = mesh.metric_functions_cache.edge.ζ̂xₖ₊½(t, ξηζ..., params)
    ζ̂_yₖ₊½ = mesh.metric_functions_cache.edge.ζ̂yₖ₊½(t, ξηζ..., params)
    ζ̂_zₖ₊½ = mesh.metric_functions_cache.edge.ζ̂zₖ₊½(t, ξηζ..., params)

    mesh.edge_metrics.k₊½.J[I] = J
    mesh.edge_metrics.k₊½.ξ̂.x₁[I] = ξ̂_xₖ₊½ * (abs(ξ̂_xₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.ξ̂.x₂[I] = ξ̂_yₖ₊½ * (abs(ξ̂_yₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.ξ̂.x₃[I] = ξ̂_zₖ₊½ * (abs(ξ̂_zₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.η̂.x₁[I] = η̂_xₖ₊½ * (abs(η̂_xₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.η̂.x₂[I] = η̂_yₖ₊½ * (abs(η̂_yₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.η̂.x₃[I] = η̂_zₖ₊½ * (abs(η̂_zₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.ζ̂.x₁[I] = ζ̂_xₖ₊½ * (abs(ζ̂_xₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.ζ̂.x₂[I] = ζ̂_yₖ₊½ * (abs(ζ̂_yₖ₊½) >= eps())
    mesh.edge_metrics.k₊½.ζ̂.x₃[I] = ζ̂_zₖ₊½ * (abs(ζ̂_zₖ₊½) >= eps())

    ξ_xₖ₊½, η_xₖ₊½, ζ_xₖ₊½, ξ_yₖ₊½, η_yₖ₊½, ζ_yₖ₊½, ξ_zₖ₊½, η_zₖ₊½, ζ_zₖ₊½ = mesh.metric_functions_cache.edge.Jinv_ₖ₊½(
      t, ξηζ..., params
    )

    mesh.edge_metrics.k₊½.ξ.x₁[I] = ifelse(
      isfinite(ξ_xₖ₊½), ξ_xₖ₊½ * (abs(ξ_xₖ₊½) >= eps()), zero(ξ_xₖ₊½)
    )
    mesh.edge_metrics.k₊½.ξ.x₂[I] = ifelse(
      isfinite(ξ_yₖ₊½), ξ_yₖ₊½ * (abs(ξ_yₖ₊½) >= eps()), zero(ξ_yₖ₊½)
    )
    mesh.edge_metrics.k₊½.ξ.x₃[I] = ifelse(
      isfinite(ξ_zₖ₊½), ξ_zₖ₊½ * (abs(ξ_zₖ₊½) >= eps()), zero(ξ_zₖ₊½)
    )
    mesh.edge_metrics.k₊½.η.x₁[I] = ifelse(
      isfinite(η_xₖ₊½), η_xₖ₊½ * (abs(η_xₖ₊½) >= eps()), zero(η_xₖ₊½)
    )
    mesh.edge_metrics.k₊½.η.x₂[I] = ifelse(
      isfinite(η_yₖ₊½), η_yₖ₊½ * (abs(η_yₖ₊½) >= eps()), zero(η_yₖ₊½)
    )
    mesh.edge_metrics.k₊½.η.x₃[I] = ifelse(
      isfinite(η_zₖ₊½), η_zₖ₊½ * (abs(η_zₖ₊½) >= eps()), zero(η_zₖ₊½)
    )
    mesh.edge_metrics.k₊½.ζ.x₁[I] = ifelse(
      isfinite(ζ_xₖ₊½), ζ_xₖ₊½ * (abs(ζ_xₖ₊½) >= eps()), zero(ζ_xₖ₊½)
    )
    mesh.edge_metrics.k₊½.ζ.x₂[I] = ifelse(
      isfinite(ζ_yₖ₊½), ζ_yₖ₊½ * (abs(ζ_yₖ₊½) >= eps()), zero(ζ_yₖ₊½)
    )
    mesh.edge_metrics.k₊½.ζ.x₃[I] = ifelse(
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