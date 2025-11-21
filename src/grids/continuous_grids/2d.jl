
function ContinuousCurvilinearGrid2D(
  x::Function,
  y::Function,
  mapping_function_parameters::NamedTuple,
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

  ni_nodes, nj_nodes = size(iterators.node.full)
  node_coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes)),
    y=KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes)),
  ))

  ni_cells, nj_cells = size(iterators.cell.full)
  centroid_coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells)),
    y=KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells)),
  ))

  discr_scheme = GradientDiscretizationScheme(order; use_cache=false)
  metric_cache = MetricCache(x, y, diff_backend)
  mapping_funcs = (; x, y)
  mesh = ContinuousCurvilinearGrid2D{
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

function coord(mesh::ContinuousCurvilinearGrid2D, t, i, j, params)
  ξη = (i, j) .- mesh.nhalo

  return @SVector [
    mesh.mapping_functions.x(t, ξη..., params), mesh.mapping_functions.y(t, ξη..., params)
  ]
end

function centroid(mesh::ContinuousCurvilinearGrid2D, t, i, j, params)
  ξη = (i, j) .- mesh.nhalo .+ 0.5

  return @SVector [
    mesh.mapping_functions.x(t, ξη..., params), mesh.mapping_functions.y(t, ξη..., params)
  ]
end

function compute_node_coordinates!(mesh::ContinuousCurvilinearGrid2D, t, params)
  nhalo = mesh.nhalo

  x = mesh.mapping_functions.x
  y = mesh.mapping_functions.y
  @threads for I in mesh.iterators.cell.full
    Iglobal = mesh.iterators.global_domain.cell.full[I]
    # account for halo cells 
    ξη = Iglobal.I .- nhalo
    mesh.node_coordinates.x[I] = x(t, ξη..., params)
    mesh.node_coordinates.y[I] = y(t, ξη..., params)
  end

  return nothing
end

function compute_centroid_coordinates!(mesh::ContinuousCurvilinearGrid2D, t, params)
  nhalo = mesh.nhalo

  x = mesh.mapping_functions.x
  y = mesh.mapping_functions.y

  @threads for I in mesh.iterators.cell.full
    Iglobal = mesh.iterators.global_domain.cell.full[I]
    # account for halo cells and centroid offset
    ξη = Iglobal.I .- nhalo .+ 0.5 # centroid
    mesh.centroid_coordinates.x[I] = x(t, ξη..., params)
    mesh.centroid_coordinates.y[I] = y(t, ξη..., params)
  end

  return nothing
end

function compute_cell_metrics!(mesh::ContinuousCurvilinearGrid2D, t, params)
  nhalo = mesh.iterators.nhalo

  @threads for I in mesh.iterators.cell.full
    Iglobal = mesh.iterators.global_domain.cell.full[I]
    # account for halo cells and centroid offset
    ξη = Iglobal.I .- nhalo .+ 0.5 # centroid

    jac_matrix = mesh.metric_functions_cache.forward.jacobian(t, ξη..., params)
    J = det(jac_matrix)

    xξ, yξ, xη, yη = jac_matrix
    ξx, ηx, ξy, ηy = mesh.metric_functions_cache.inverse.Jinv(t, ξη..., params)

    mesh.cell_center_metrics.J[I] = J
    mesh.cell_center_metrics.x₁.ξ[I] = xξ
    mesh.cell_center_metrics.x₁.η[I] = xη
    mesh.cell_center_metrics.x₂.ξ[I] = yξ
    mesh.cell_center_metrics.x₂.η[I] = yη

    mesh.cell_center_metrics.ξ.x₁[I] = ξx
    mesh.cell_center_metrics.ξ.x₂[I] = ξy
    mesh.cell_center_metrics.η.x₁[I] = ηx
    mesh.cell_center_metrics.η.x₂[I] = ηy

    # hatted inverse metrics
    mesh.cell_center_metrics.ξ̂.x₁[I] = ξx * J
    mesh.cell_center_metrics.ξ̂.x₂[I] = ξy * J
    mesh.cell_center_metrics.η̂.x₁[I] = ηx * J
    mesh.cell_center_metrics.η̂.x₂[I] = ηy * J
  end

  if any(mesh.cell_center_metrics.J .<= 0)
    @error "Invalid Jacobians detected (i.e. negative or zero volumes)"
  end
end

function compute_edge_metrics!(mesh::ContinuousCurvilinearGrid2D, t, params)
  nhalo = mesh.iterators.nhalo

  @threads for I in mesh.iterators.cell.full
    Iglobal = mesh.iterators.global_domain.cell.full[I]
    ξη = Iglobal.I .- nhalo .+ (1 / 2) # centroid index

    # J = Jᵢ₊½(ξη...)
    ξ_xᵢ₊½, η_xᵢ₊½, ξ_yᵢ₊½, η_yᵢ₊½ = mesh.metric_functions_cache.edge.Jinv_ᵢ₊½(
      t, ξη..., params
    )
    ξ̂_xᵢ₊½, η̂_xᵢ₊½, ξ̂_yᵢ₊½, η̂_yᵢ₊½ = mesh.metric_functions_cache.edge.norm_Jinv_ᵢ₊½(
      t, ξη..., params
    )

    mesh.edge_metrics.i₊½.ξ̂.x₁[I] = ξ̂_xᵢ₊½
    mesh.edge_metrics.i₊½.ξ̂.x₂[I] = ξ̂_yᵢ₊½
    mesh.edge_metrics.i₊½.η̂.x₁[I] = η̂_xᵢ₊½
    mesh.edge_metrics.i₊½.η̂.x₂[I] = η̂_yᵢ₊½

    mesh.edge_metrics.i₊½.ξ.x₁[I] = ξ_xᵢ₊½
    mesh.edge_metrics.i₊½.ξ.x₂[I] = ξ_yᵢ₊½

    mesh.edge_metrics.i₊½.η.x₁[I] = η_xᵢ₊½
    mesh.edge_metrics.i₊½.η.x₂[I] = η_yᵢ₊½
    # end

    # @threads for I in mesh.iterators.cell.full
    #   Iglobal = mesh.iterators.global_domain.cell.full[I]
    # ξη = Iglobal.I .- nhalo .+ (1 / 2) # centroid index

    # J = Jⱼ₊½(t, ξη..., params)
    ξ_xⱼ₊½, η_xⱼ₊½, ξ_yⱼ₊½, η_yⱼ₊½ = mesh.metric_functions_cache.edge.Jinv_ⱼ₊½(
      t, ξη..., params
    )
    ξ̂_xⱼ₊½, η̂_xⱼ₊½, ξ̂_yⱼ₊½, η̂_yⱼ₊½ = mesh.metric_functions_cache.edge.norm_Jinv_ⱼ₊½(
      t, ξη..., params
    )

    mesh.edge_metrics.j₊½.ξ̂.x₁[I] = ξ̂_xⱼ₊½
    mesh.edge_metrics.j₊½.ξ̂.x₂[I] = ξ̂_yⱼ₊½
    mesh.edge_metrics.j₊½.η̂.x₁[I] = η̂_xⱼ₊½
    mesh.edge_metrics.j₊½.η̂.x₂[I] = η̂_yⱼ₊½

    mesh.edge_metrics.j₊½.ξ.x₁[I] = ξ_xⱼ₊½
    mesh.edge_metrics.j₊½.ξ.x₂[I] = ξ_yⱼ₊½

    mesh.edge_metrics.j₊½.η.x₁[I] = η_xⱼ₊½
    mesh.edge_metrics.j₊½.η.x₂[I] = η_yⱼ₊½
  end

  return nothing
end

function jacobian_matrix(mesh::ContinuousCurvilinearGrid2D{T}, t, i, j) where {T}
  p = mesh.mapping_function_params
  return mesh.metric_functions_cache.forward.jacobian(t, i, j, p)::SMatrix{2,2,T,4}
end