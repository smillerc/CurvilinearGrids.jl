
function ContinuousCurvilinearGrid1D(
  x::Function,
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

  ni_nodes, = size(iterators.node.full)
  node_coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, (ni_nodes,)),
    y=KernelAbstractions.zeros(backend, T, (ni_nodes,)),
  ))

  ni_cells, = size(iterators.cell.full)
  centroid_coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, (ni_cells,)),
    y=KernelAbstractions.zeros(backend, T, (ni_cells,)),
  ))

  discr_scheme = GradientDiscretizationScheme(order; use_cache=false)
  metric_cache = MetricCache(x, diff_backend)
  mapping_funcs = (; x)

  discr_scheme = GradientDiscretizationScheme(order; use_cache=false)

  mesh = ContinuousCurvilinearGrid1D{
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
    compute_cell_metrics!(mesh)
    compute_edge_metrics!(mesh)
  end

  return mesh
end

function compute_node_coordinates!(mesh::ContinuousCurvilinearGrid1D, t, params)
  nhalo = mesh.nhalo

  x = mesh.mapping_functions.x
  @threads for I in mesh.iterators.cell.full
    Iglobal = mesh.iterators.global_domain.cell.full[I]
    # account for halo cells 
    ξη = Iglobal.I .- nhalo
    mesh.node_coordinates.x[I] = x(t, ξη..., params)
  end

  return nothing
end

function compute_centroid_coordinates!(mesh::ContinuousCurvilinearGrid1D, t, params)
  nhalo = mesh.nhalo

  x = mesh.mapping_functions.x

  @threads for I in mesh.iterators.cell.full
    Iglobal = mesh.iterators.global_domain.cell.full[I]
    # account for halo cells and centroid offset
    ξη = Iglobal.I .- nhalo .+ 0.5 # centroid
    mesh.centroid_coordinates.x[I] = x(t, ξη..., params)
  end

  return nothing
end

function compute_cell_metrics!(mesh::ContinuousCurvilinearGrid1D)
  nhalo = mesh.iterators.nhalo

  @threads for I in mesh.iterators.cell.full
    Iglobal = mesh.iterators.global_domain.cell.full[I]
    # account for halo cells and centroid offset
    ξη = Iglobal.I .- nhalo .+ 0.5 # centroid

    jac_matrix = mesh.metric_functions_cache.forward.jacobian(ξη...)
    J = det(jac_matrix)

    xξ, = jac_matrix
    ξx, = mesh.metric_functions_cache.inverse.Jinv(ξη...)

    mesh.cell_center_metrics.J[i] = J
    mesh.cell_center_metrics.x₁.ξ[i] = xξ
    mesh.cell_center_metrics.ξ.x₁[i] = ξx
    mesh.cell_center_metrics.ξ̂.x₁[i] = ξx * J
  end

  if any(mesh.cell_center_metrics.J .<= 0)
    @error "Invalid Jacobians detected (i.e. negative or zero volumes)"
  end
end

function compute_edge_metrics!(mesh::ContinuousCurvilinearGrid1D)
  nhalo = mesh.iterators.nhalo

  # Jᵢ₊½,  = edge_functions_2d(mesh.metric_functions_cache.forward.J, mesh.diff_backend)
  Jinv_ᵢ₊½, = edge_functions_1d(mesh.metric_functions_cache.inverse.Jinv, mesh.diff_backend)
  norm_Jinv_ᵢ₊½, = edge_functions_1d(
    mesh.metric_functions_cache.inverse.Jinv_norm, mesh.diff_backend
  )

  @threads for I in mesh.iterators.cell.full
    Iglobal = mesh.iterators.global_domain.cell.full[I]
    # account for halo cells and centroid offset
    ξη = Iglobal.I .- nhalo .+ 0.5 # centroid

    # J = Jᵢ₊½(ξη...)
    ξ_xᵢ₊½, = Jinv_ᵢ₊½(ξη...)
    ξ̂_xᵢ₊½, = norm_Jinv_ᵢ₊½(ξη...)

    mesh.edge_metrics.i₊½.ξ̂.x₁[i] = ξ̂_xᵢ₊½
    mesh.edge_metrics.i₊½.ξ.x₁[i] = ξ_xᵢ₊½
  end

  return nothing
end
