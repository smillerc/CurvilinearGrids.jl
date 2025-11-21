
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
  global_cell_indices::Union{Nothing,CartesianIndices{3}}=nothing,
)
  GradientDiscretizationScheme, order, _, nhalo, scheme_name = get_gradient_discretization_scheme(
    discretization_scheme
  )
  iterators = get_iterators(celldims, nhalo, global_cell_indices)

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
  mapping_funcs = (; x, y, z)
  mesh = ContinuousCurvilinearGrid3D{
    T,
    typeof(node_coords),
    typeof(mapping_funcs),
    typeof(metric_cache),
    typeof(edge_metrics),
    typeof(cell_center_metrics),
    typeof(iterators),
    typeof(backend),
    typeof(diff_backend),
    typeof(discr_scheme),
  }(
    node_coords,
    centroid_coords,
    mapping_funcs,
    metric_cache,
    edge_metrics,
    cell_center_metrics,
    backend,
    diff_backend,
    nhalo,
    discr_scheme,
    scheme_name,
    iterators,
  )

  compute_node_coordinates!(mesh, t, mapping_function_parameters)
  compute_centroid_coordinates!(mesh, t, mapping_function_parameters)

  if compute_metrics
    compute_cell_metrics!(mesh, t, mapping_function_parameters)
    compute_edge_metrics!(mesh, t, mapping_function_parameters)
  end

  return mesh
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

function compute_cell_metrics!(mesh, t, params)
  nhalo = mesh.iterators.nhalo
  ϵ = eps()

  # @info "Computing Cell Metrics"
  Threads.@threads :greedy for I in mesh.iterators.cell.full
    Iglobal = mesh.iterators.global_domain.cell.full[I]
    # account for halo cells and centroid offset
    ξηζ = Iglobal.I .- nhalo .+ 0.5 # centroid

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

    mesh.cell_center_metrics.x₁.ξ[I] = xξ * (abs(xξ) >= ϵ)
    mesh.cell_center_metrics.x₁.η[I] = xη * (abs(xη) >= ϵ)
    mesh.cell_center_metrics.x₁.ζ[I] = xζ * (abs(xζ) >= ϵ)
    mesh.cell_center_metrics.x₂.ξ[I] = yξ * (abs(yξ) >= ϵ)
    mesh.cell_center_metrics.x₂.η[I] = yη * (abs(yη) >= ϵ)
    mesh.cell_center_metrics.x₂.ζ[I] = yζ * (abs(yζ) >= ϵ)
    mesh.cell_center_metrics.x₃.ξ[I] = zξ * (abs(zξ) >= ϵ)
    mesh.cell_center_metrics.x₃.η[I] = zη * (abs(zη) >= ϵ)
    mesh.cell_center_metrics.x₃.ζ[I] = zζ * (abs(zζ) >= ϵ)

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

    # ξ̂_x = ξ̂_x # * (abs(ξ̂_x) >= ϵ)
    # ξ̂_y = ξ̂_y # * (abs(ξ̂_y) >= ϵ)
    # ξ̂_z = ξ̂_z # * (abs(ξ̂_z) >= ϵ)
    # η̂_x = η̂_x # * (abs(η̂_x) >= ϵ)
    # η̂_y = η̂_y # * (abs(η̂_y) >= ϵ)
    # η̂_z = η̂_z # * (abs(η̂_z) >= ϵ)
    # ζ̂_x = ζ̂_x # * (abs(ζ̂_x) >= ϵ)
    # ζ̂_y = ζ̂_y # * (abs(ζ̂_y) >= ϵ)
    # ζ̂_z = ζ̂_z # * (abs(ζ̂_z) >= ϵ)

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
  ϵ = eps()
  # @info "Computing Edge Metrics"
  Threads.@threads :greedy for I in mesh.iterators.cell.full
    Iglobal = mesh.iterators.global_domain.cell.full[I]
    ξηζ = Iglobal.I .- nhalo .+ (1 / 2) # centroid index

    # Jᵢ₊½ = mesh.metric_functions_cache.edge.Jᵢ₊½(t, ξηζ..., params)
    # Jⱼ₊½ = mesh.metric_functions_cache.edge.Jⱼ₊½(t, ξηζ..., params)
    # Jₖ₊½ = mesh.metric_functions_cache.edge.Jₖ₊½(t, ξηζ..., params)

    ξ̂xᵢ₊½ = mesh.metric_functions_cache.edge.ξ̂xᵢ₊½(t, ξηζ..., params)
    ξ̂yᵢ₊½ = mesh.metric_functions_cache.edge.ξ̂yᵢ₊½(t, ξηζ..., params)
    ξ̂zᵢ₊½ = mesh.metric_functions_cache.edge.ξ̂zᵢ₊½(t, ξηζ..., params)
    η̂xᵢ₊½ = mesh.metric_functions_cache.edge.η̂xᵢ₊½(t, ξηζ..., params)
    η̂yᵢ₊½ = mesh.metric_functions_cache.edge.η̂yᵢ₊½(t, ξηζ..., params)
    η̂zᵢ₊½ = mesh.metric_functions_cache.edge.η̂zᵢ₊½(t, ξηζ..., params)
    ζ̂xᵢ₊½ = mesh.metric_functions_cache.edge.ζ̂xᵢ₊½(t, ξηζ..., params)
    ζ̂yᵢ₊½ = mesh.metric_functions_cache.edge.ζ̂yᵢ₊½(t, ξηζ..., params)
    ζ̂zᵢ₊½ = mesh.metric_functions_cache.edge.ζ̂zᵢ₊½(t, ξηζ..., params)

    ξ̂xⱼ₊½ = mesh.metric_functions_cache.edge.ξ̂xⱼ₊½(t, ξηζ..., params)
    ξ̂yⱼ₊½ = mesh.metric_functions_cache.edge.ξ̂yⱼ₊½(t, ξηζ..., params)
    ξ̂zⱼ₊½ = mesh.metric_functions_cache.edge.ξ̂zⱼ₊½(t, ξηζ..., params)
    η̂xⱼ₊½ = mesh.metric_functions_cache.edge.η̂xⱼ₊½(t, ξηζ..., params)
    η̂yⱼ₊½ = mesh.metric_functions_cache.edge.η̂yⱼ₊½(t, ξηζ..., params)
    η̂zⱼ₊½ = mesh.metric_functions_cache.edge.η̂zⱼ₊½(t, ξηζ..., params)
    ζ̂xⱼ₊½ = mesh.metric_functions_cache.edge.ζ̂xⱼ₊½(t, ξηζ..., params)
    ζ̂yⱼ₊½ = mesh.metric_functions_cache.edge.ζ̂yⱼ₊½(t, ξηζ..., params)
    ζ̂zⱼ₊½ = mesh.metric_functions_cache.edge.ζ̂zⱼ₊½(t, ξηζ..., params)

    ξ̂xₖ₊½ = mesh.metric_functions_cache.edge.ξ̂xₖ₊½(t, ξηζ..., params)
    ξ̂yₖ₊½ = mesh.metric_functions_cache.edge.ξ̂yₖ₊½(t, ξηζ..., params)
    ξ̂zₖ₊½ = mesh.metric_functions_cache.edge.ξ̂zₖ₊½(t, ξηζ..., params)
    η̂xₖ₊½ = mesh.metric_functions_cache.edge.η̂xₖ₊½(t, ξηζ..., params)
    η̂yₖ₊½ = mesh.metric_functions_cache.edge.η̂yₖ₊½(t, ξηζ..., params)
    η̂zₖ₊½ = mesh.metric_functions_cache.edge.η̂zₖ₊½(t, ξηζ..., params)
    ζ̂xₖ₊½ = mesh.metric_functions_cache.edge.ζ̂xₖ₊½(t, ξηζ..., params)
    ζ̂yₖ₊½ = mesh.metric_functions_cache.edge.ζ̂yₖ₊½(t, ξηζ..., params)
    ζ̂zₖ₊½ = mesh.metric_functions_cache.edge.ζ̂zₖ₊½(t, ξηζ..., params)

    # mesh.edge_metrics.i₊½.J[I] = J
    mesh.edge_metrics.i₊½.ξ̂.x₁[I] = ξ̂xᵢ₊½ * (abs(ξ̂xᵢ₊½) >= ϵ)
    mesh.edge_metrics.i₊½.ξ̂.x₂[I] = ξ̂yᵢ₊½ * (abs(ξ̂yᵢ₊½) >= ϵ)
    mesh.edge_metrics.i₊½.η̂.x₁[I] = η̂xᵢ₊½ * (abs(η̂xᵢ₊½) >= ϵ)
    mesh.edge_metrics.i₊½.η̂.x₂[I] = η̂yᵢ₊½ * (abs(η̂yᵢ₊½) >= ϵ)
    mesh.edge_metrics.i₊½.ζ̂.x₁[I] = ζ̂xᵢ₊½ * (abs(ζ̂xᵢ₊½) >= ϵ)
    mesh.edge_metrics.i₊½.ζ̂.x₂[I] = ζ̂yᵢ₊½ * (abs(ζ̂yᵢ₊½) >= ϵ)
    mesh.edge_metrics.i₊½.ξ̂.x₃[I] = ξ̂zᵢ₊½ * (abs(ξ̂zᵢ₊½) >= ϵ)
    mesh.edge_metrics.i₊½.η̂.x₃[I] = η̂zᵢ₊½ * (abs(η̂zᵢ₊½) >= ϵ)
    mesh.edge_metrics.i₊½.ζ̂.x₃[I] = ζ̂zᵢ₊½ * (abs(ζ̂zᵢ₊½) >= ϵ)

    # mesh.edge_metrics.j₊½.J[I] = J
    mesh.edge_metrics.j₊½.ξ̂.x₁[I] = ξ̂xⱼ₊½ * (abs(ξ̂xⱼ₊½) >= ϵ)
    mesh.edge_metrics.j₊½.ξ̂.x₂[I] = ξ̂yⱼ₊½ * (abs(ξ̂yⱼ₊½) >= ϵ)
    mesh.edge_metrics.j₊½.η̂.x₁[I] = η̂xⱼ₊½ * (abs(η̂xⱼ₊½) >= ϵ)
    mesh.edge_metrics.j₊½.η̂.x₂[I] = η̂yⱼ₊½ * (abs(η̂yⱼ₊½) >= ϵ)
    mesh.edge_metrics.j₊½.ζ̂.x₁[I] = ζ̂xⱼ₊½ * (abs(ζ̂xⱼ₊½) >= ϵ)
    mesh.edge_metrics.j₊½.ζ̂.x₂[I] = ζ̂yⱼ₊½ * (abs(ζ̂yⱼ₊½) >= ϵ)
    mesh.edge_metrics.j₊½.ξ̂.x₃[I] = ξ̂zⱼ₊½ * (abs(ξ̂zⱼ₊½) >= ϵ)
    mesh.edge_metrics.j₊½.η̂.x₃[I] = η̂zⱼ₊½ * (abs(η̂zⱼ₊½) >= ϵ)
    mesh.edge_metrics.j₊½.ζ̂.x₃[I] = ζ̂zⱼ₊½ * (abs(ζ̂zⱼ₊½) >= ϵ)

    # mesh.edge_metrics.k₊½.J[I] = J
    mesh.edge_metrics.k₊½.ξ̂.x₁[I] = ξ̂xₖ₊½ * (abs(ξ̂xₖ₊½) >= ϵ)
    mesh.edge_metrics.k₊½.ξ̂.x₂[I] = ξ̂yₖ₊½ * (abs(ξ̂yₖ₊½) >= ϵ)
    mesh.edge_metrics.k₊½.η̂.x₁[I] = η̂xₖ₊½ * (abs(η̂xₖ₊½) >= ϵ)
    mesh.edge_metrics.k₊½.η̂.x₂[I] = η̂yₖ₊½ * (abs(η̂yₖ₊½) >= ϵ)
    mesh.edge_metrics.k₊½.ζ̂.x₁[I] = ζ̂xₖ₊½ * (abs(ζ̂xₖ₊½) >= ϵ)
    mesh.edge_metrics.k₊½.ζ̂.x₂[I] = ζ̂yₖ₊½ * (abs(ζ̂yₖ₊½) >= ϵ)
    mesh.edge_metrics.k₊½.ξ̂.x₃[I] = ξ̂zₖ₊½ * (abs(ξ̂zₖ₊½) >= ϵ)
    mesh.edge_metrics.k₊½.η̂.x₃[I] = η̂zₖ₊½ * (abs(η̂zₖ₊½) >= ϵ)
    mesh.edge_metrics.k₊½.ζ̂.x₃[I] = ζ̂zₖ₊½ * (abs(ζ̂zₖ₊½) >= ϵ)

    ξ_xᵢ₊½, η_xᵢ₊½, ζ_xᵢ₊½, ξ_yᵢ₊½, η_yᵢ₊½, ζ_yᵢ₊½, ξ_zᵢ₊½, η_zᵢ₊½, ζ_zᵢ₊½ = mesh.metric_functions_cache.edge.Jinvᵢ₊½(
      t, ξηζ..., params
    )

    ξ_xⱼ₊½, η_xⱼ₊½, ζ_xⱼ₊½, ξ_yⱼ₊½, η_yⱼ₊½, ζ_yⱼ₊½, ξ_zⱼ₊½, η_zⱼ₊½, ζ_zⱼ₊½ = mesh.metric_functions_cache.edge.Jinvⱼ₊½(
      t, ξηζ..., params
    )

    ξ_xₖ₊½, η_xₖ₊½, ζ_xₖ₊½, ξ_yₖ₊½, η_yₖ₊½, ζ_yₖ₊½, ξ_zₖ₊½, η_zₖ₊½, ζ_zₖ₊½ = mesh.metric_functions_cache.edge.Jinvₖ₊½(
      t, ξηζ..., params
    )
    #! format: off
    mesh.edge_metrics.i₊½.ξ.x₁[I] = ifelse(isfinite(ξ_xᵢ₊½), ξ_xᵢ₊½ * (abs(ξ_xᵢ₊½) >= ϵ), zero(ξ_xᵢ₊½))
    mesh.edge_metrics.i₊½.ξ.x₂[I] = ifelse(isfinite(ξ_yᵢ₊½), ξ_yᵢ₊½ * (abs(ξ_yᵢ₊½) >= ϵ), zero(ξ_yᵢ₊½))
    mesh.edge_metrics.i₊½.ξ.x₃[I] = ifelse(isfinite(ξ_zᵢ₊½), ξ_zᵢ₊½ * (abs(ξ_zᵢ₊½) >= ϵ), zero(ξ_zᵢ₊½))
    mesh.edge_metrics.i₊½.η.x₁[I] = ifelse(isfinite(η_xᵢ₊½), η_xᵢ₊½ * (abs(η_xᵢ₊½) >= ϵ), zero(η_xᵢ₊½))
    mesh.edge_metrics.i₊½.η.x₂[I] = ifelse(isfinite(η_yᵢ₊½), η_yᵢ₊½ * (abs(η_yᵢ₊½) >= ϵ), zero(η_yᵢ₊½))
    mesh.edge_metrics.i₊½.η.x₃[I] = ifelse(isfinite(η_zᵢ₊½), η_zᵢ₊½ * (abs(η_zᵢ₊½) >= ϵ), zero(η_zᵢ₊½))
    mesh.edge_metrics.i₊½.ζ.x₁[I] = ifelse(isfinite(ζ_xᵢ₊½), ζ_xᵢ₊½ * (abs(ζ_xᵢ₊½) >= ϵ), zero(ζ_xᵢ₊½))
    mesh.edge_metrics.i₊½.ζ.x₂[I] = ifelse(isfinite(ζ_yᵢ₊½), ζ_yᵢ₊½ * (abs(ζ_yᵢ₊½) >= ϵ), zero(ζ_yᵢ₊½))
    mesh.edge_metrics.i₊½.ζ.x₃[I] = ifelse(isfinite(ζ_zᵢ₊½), ζ_zᵢ₊½ * (abs(ζ_zᵢ₊½) >= ϵ), zero(ζ_zᵢ₊½))


    mesh.edge_metrics.j₊½.ξ.x₁[I] = ifelse(isfinite(ξ_xⱼ₊½), ξ_xⱼ₊½ * (abs(ξ_xⱼ₊½) >= ϵ), zero(ξ_xⱼ₊½))
    mesh.edge_metrics.j₊½.ξ.x₂[I] = ifelse(isfinite(ξ_yⱼ₊½), ξ_yⱼ₊½ * (abs(ξ_yⱼ₊½) >= ϵ), zero(ξ_yⱼ₊½))
    mesh.edge_metrics.j₊½.ξ.x₃[I] = ifelse(isfinite(ξ_zⱼ₊½), ξ_zⱼ₊½ * (abs(ξ_zⱼ₊½) >= ϵ), zero(ξ_zⱼ₊½))
    mesh.edge_metrics.j₊½.η.x₁[I] = ifelse(isfinite(η_xⱼ₊½), η_xⱼ₊½ * (abs(η_xⱼ₊½) >= ϵ), zero(η_xⱼ₊½))
    mesh.edge_metrics.j₊½.η.x₂[I] = ifelse(isfinite(η_yⱼ₊½), η_yⱼ₊½ * (abs(η_yⱼ₊½) >= ϵ), zero(η_yⱼ₊½))
    mesh.edge_metrics.j₊½.η.x₃[I] = ifelse(isfinite(η_zⱼ₊½), η_zⱼ₊½ * (abs(η_zⱼ₊½) >= ϵ), zero(η_zⱼ₊½))
    mesh.edge_metrics.j₊½.ζ.x₁[I] = ifelse(isfinite(ζ_xⱼ₊½), ζ_xⱼ₊½ * (abs(ζ_xⱼ₊½) >= ϵ), zero(ζ_xⱼ₊½))
    mesh.edge_metrics.j₊½.ζ.x₂[I] = ifelse(isfinite(ζ_yⱼ₊½), ζ_yⱼ₊½ * (abs(ζ_yⱼ₊½) >= ϵ), zero(ζ_yⱼ₊½))
    mesh.edge_metrics.j₊½.ζ.x₃[I] = ifelse(isfinite(ζ_zⱼ₊½), ζ_zⱼ₊½ * (abs(ζ_zⱼ₊½) >= ϵ), zero(ζ_zⱼ₊½))


    mesh.edge_metrics.k₊½.ξ.x₁[I] = ifelse(isfinite(ξ_xₖ₊½), ξ_xₖ₊½ * (abs(ξ_xₖ₊½) >= ϵ), zero(ξ_xₖ₊½))
    mesh.edge_metrics.k₊½.ξ.x₂[I] = ifelse(isfinite(ξ_yₖ₊½), ξ_yₖ₊½ * (abs(ξ_yₖ₊½) >= ϵ), zero(ξ_yₖ₊½))
    mesh.edge_metrics.k₊½.ξ.x₃[I] = ifelse(isfinite(ξ_zₖ₊½), ξ_zₖ₊½ * (abs(ξ_zₖ₊½) >= ϵ), zero(ξ_zₖ₊½))
    mesh.edge_metrics.k₊½.η.x₁[I] = ifelse(isfinite(η_xₖ₊½), η_xₖ₊½ * (abs(η_xₖ₊½) >= ϵ), zero(η_xₖ₊½))
    mesh.edge_metrics.k₊½.η.x₂[I] = ifelse(isfinite(η_yₖ₊½), η_yₖ₊½ * (abs(η_yₖ₊½) >= ϵ), zero(η_yₖ₊½))
    mesh.edge_metrics.k₊½.η.x₃[I] = ifelse(isfinite(η_zₖ₊½), η_zₖ₊½ * (abs(η_zₖ₊½) >= ϵ), zero(η_zₖ₊½))
    mesh.edge_metrics.k₊½.ζ.x₁[I] = ifelse(isfinite(ζ_xₖ₊½), ζ_xₖ₊½ * (abs(ζ_xₖ₊½) >= ϵ), zero(ζ_xₖ₊½))
    mesh.edge_metrics.k₊½.ζ.x₂[I] = ifelse(isfinite(ζ_yₖ₊½), ζ_yₖ₊½ * (abs(ζ_yₖ₊½) >= ϵ), zero(ζ_yₖ₊½))
    mesh.edge_metrics.k₊½.ζ.x₃[I] = ifelse(isfinite(ζ_zₖ₊½), ζ_zₖ₊½ * (abs(ζ_zₖ₊½) >= ϵ), zero(ζ_zₖ₊½))
    #! format: on
  end

  return nothing
end

# function cell_center_forward_metrics(mesh, (i, j, k))

#   # instead of grabbing all the individual forward metrics, just 
#   # call the jacobian and take advantage of the vector-mode AD

#   ξηζ = (i, j, k) .- mesh.iterators.nhalo .+ (1 / 2) # centroid
#   J⃗ = mesh.metric_functions_cache.forward.jacobian(ξηζ...)

#   #   J⃗ = [
#   #     xξ xη xζ
#   #     yξ yη yζ
#   #     zξ zη zζ
#   #   ]

#   return (;
#     xξ=J⃗[1], xη=J⃗[4], xζ=J⃗[7], yξ=J⃗[2], yη=J⃗[5], yζ=J⃗[8], zξ=J⃗[3], zη=J⃗[6], zζ=J⃗[9]
#   )
# end

# function cell_center_jacobian_matrix(mesh, (i, j, k))
#   ξηζ = (i, j, k) .- mesh.iterators.nhalo .+ (1 / 2) # centroid
#   return mesh.metric_functions_cache.forward.jacobian(ξηζ...)
# end

# function cell_center_conservative_metrics(mesh, (i, j, k))

#   # instead of grabbing all the individual forward metrics, just 
#   # call the jacobian and take advantage of the vector-mode AD

#   ξηζ = (i, j, k) .- mesh.iterators.nhalo .+ (1 / 2) # centroid

#   return (;
#     ξ̂x=mesh.metric_functions_cache.inverse.ξ̂x(ξηζ...),
#     ξ̂y=mesh.metric_functions_cache.inverse.ξ̂y(ξηζ...),
#     ξ̂z=mesh.metric_functions_cache.inverse.ξ̂z(ξηζ...),
#     η̂x=mesh.metric_functions_cache.inverse.η̂x(ξηζ...),
#     η̂y=mesh.metric_functions_cache.inverse.η̂y(ξηζ...),
#     η̂z=mesh.metric_functions_cache.inverse.η̂z(ξηζ...),
#     ζ̂x=mesh.metric_functions_cache.inverse.ζ̂x(ξηζ...),
#     ζ̂y=mesh.metric_functions_cache.inverse.ζ̂y(ξηζ...),
#     ζ̂z=mesh.metric_functions_cache.inverse.ζ̂z(ξηζ...),
#   )
# end

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