
abstract type AbstractContinuousCurvilinearGrid1D <: AbstractCurvilinearGrid1D end

struct ContinuousCurvilinearGrid1D{A,B,C,EM,CM,E,BE,DBE,DS} <:
       AbstractContinuousCurvilinearGrid1D
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

function ContinuousCurvilinearGrid1D(
  x::Function,
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
  xyz_n = node_coordinates(x, iterators, backend, T)
  xyz_c = centroid_coordinates(x, iterators, backend, T)

  discr_scheme = GradientDiscretizationScheme(order; use_cache=false)

  mesh = ContinuousCurvilinearGrid1D(
    xyz_n,
    xyz_c,
    (; x),
    MetricCache(x, diff_backend),
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

function node_coordinates(x, iterators, backend, T)
  ni, = size(iterators.node.full)
  nhalo = iterators.nhalo

  coords = StructArray((x=KernelAbstractions.zeros(backend, T, (ni,)),))

  @batch for I in iterators.node.full
    i, = I.I
    coords.x[I] = x(i - nhalo)
  end

  return coords
end

function centroid_coordinates(x, iterators, backend, T)
  ni, = size(iterators.cell.full)
  nhalo = iterators.nhalo
  coords = StructArray((x=KernelAbstractions.zeros(backend, T, (ni,)),))

  @batch for I in iterators.cell.full
    i, = I.I
    coords.x[I] = x(i - nhalo + 0.5)
  end

  return coords
end

function compute_cell_metrics!(mesh::ContinuousCurvilinearGrid1D)
  nhalo = mesh.iterators.nhalo

  @threads for I in mesh.iterators.cell.full
    i, = I.I

    # account for halo cells and centroid offset
    ξη = I.I .- nhalo .+ 0.5 # centroid

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
    i, = I.I
    ξη = I.I .- nhalo .+ (1 / 2) # centroid index

    # J = Jᵢ₊½(ξη...)
    ξ_xᵢ₊½, = Jinv_ᵢ₊½(ξη...)
    ξ̂_xᵢ₊½, = norm_Jinv_ᵢ₊½(ξη...)

    mesh.edge_metrics.i₊½.ξ̂.x₁[i] = ξ̂_xᵢ₊½
    mesh.edge_metrics.i₊½.ξ.x₁[i] = ξ_xᵢ₊½
  end

  return nothing
end
