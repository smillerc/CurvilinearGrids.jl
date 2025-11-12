
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

function ContinuousCurvilinearGrid2D(
  x::Function,
  y::Function,
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
  xyz_n = node_coordinates(x, y, iterators, backend, T)
  xyz_c = centroid_coordinates(x, y, iterators, backend, T)

  discr_scheme = GradientDiscretizationScheme(order; use_cache=false)

  mesh = ContinuousCurvilinearGrid2D(
    xyz_n,
    xyz_c,
    (; x, y),
    MetricCache(x, y, diff_backend),
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
    x=KernelAbstractions.zeros(backend, T, (ni, nj)),
    y=KernelAbstractions.zeros(backend, T, (ni, nj)),
  ))

  @batch for I in iterators.cell.full
    i, j = I.I
    coords.x[I] = x(i - nhalo + 0.5, j - nhalo + 0.5)
    coords.y[I] = y(i - nhalo + 0.5, j - nhalo + 0.5)
  end

  return coords
end

function compute_cell_metrics!(mesh::ContinuousCurvilinearGrid2D)
  nhalo = mesh.iterators.nhalo

  @threads for I in mesh.iterators.cell.full
    i, j = I.I

    # account for halo cells and centroid offset
    ξη = I.I .- nhalo .+ 0.5 # centroid

    jac_matrix = mesh.metric_functions_cache.forward.jacobian(ξη...)
    J = det(jac_matrix)

    xξ, yξ, xη, yη = jac_matrix
    ξx, ηx, ξy, ηy = mesh.metric_functions_cache.inverse.Jinv(ξη...)

    mesh.cell_center_metrics.J[i, j] = J
    mesh.cell_center_metrics.x₁.ξ[i, j] = xξ
    mesh.cell_center_metrics.x₁.η[i, j] = xη
    mesh.cell_center_metrics.x₂.ξ[i, j] = yξ
    mesh.cell_center_metrics.x₂.η[i, j] = yη

    mesh.cell_center_metrics.ξ.x₁[i, j] = ξx
    mesh.cell_center_metrics.ξ.x₂[i, j] = ξy
    mesh.cell_center_metrics.η.x₁[i, j] = ηx
    mesh.cell_center_metrics.η.x₂[i, j] = ηy

    # hatted inverse metrics
    mesh.cell_center_metrics.ξ̂.x₁[i, j] = ξx * J
    mesh.cell_center_metrics.ξ̂.x₂[i, j] = ξy * J
    mesh.cell_center_metrics.η̂.x₁[i, j] = ηx * J
    mesh.cell_center_metrics.η̂.x₂[i, j] = ηy * J
  end

  if any(mesh.cell_center_metrics.J .<= 0)
    @error "Invalid Jacobians detected (i.e. negative or zero volumes)"
  end
end

function compute_edge_metrics!(mesh::ContinuousCurvilinearGrid2D)
  nhalo = mesh.iterators.nhalo

  # Jᵢ₊½, Jⱼ₊½ = edge_functions_2d(mesh.metric_functions_cache.forward.J, mesh.diff_backend)

  Jinv_ᵢ₊½, Jinv_ⱼ₊½ = edge_functions_2d(
    mesh.metric_functions_cache.inverse.Jinv, mesh.diff_backend
  )
  norm_Jinv_ᵢ₊½, norm_Jinv_ⱼ₊½ = edge_functions_2d(
    mesh.metric_functions_cache.inverse.Jinv_norm, mesh.diff_backend
  )

  @threads for I in mesh.iterators.cell.full
    i, j = I.I
    ξη = I.I .- nhalo .+ (1 / 2) # centroid index

    # J = Jᵢ₊½(ξη...)
    ξ_xᵢ₊½, η_xᵢ₊½, ξ_yᵢ₊½, η_yᵢ₊½ = Jinv_ᵢ₊½(ξη...)
    ξ̂_xᵢ₊½, η̂_xᵢ₊½, ξ̂_yᵢ₊½, η̂_yᵢ₊½ = norm_Jinv_ᵢ₊½(ξη...)

    mesh.edge_metrics.i₊½.ξ̂.x₁[i, j] = ξ̂_xᵢ₊½
    mesh.edge_metrics.i₊½.ξ̂.x₂[i, j] = ξ̂_yᵢ₊½
    mesh.edge_metrics.i₊½.η̂.x₁[i, j] = η̂_xᵢ₊½
    mesh.edge_metrics.i₊½.η̂.x₂[i, j] = η̂_yᵢ₊½

    mesh.edge_metrics.i₊½.ξ.x₁[i, j] = ξ_xᵢ₊½
    mesh.edge_metrics.i₊½.ξ.x₂[i, j] = ξ_yᵢ₊½

    mesh.edge_metrics.i₊½.η.x₁[i, j] = η_xᵢ₊½
    mesh.edge_metrics.i₊½.η.x₂[i, j] = η_yᵢ₊½
  end

  @threads for I in mesh.iterators.cell.full
    i, j = I.I
    ξη = I.I .- nhalo .+ (1 / 2) # centroid index

    # J = Jⱼ₊½(ξη...)
    ξ_xⱼ₊½, η_xⱼ₊½, ξ_yⱼ₊½, η_yⱼ₊½ = Jinv_ⱼ₊½(ξη...)
    ξ̂_xⱼ₊½, η̂_xⱼ₊½, ξ̂_yⱼ₊½, η̂_yⱼ₊½ = norm_Jinv_ⱼ₊½(ξη...)

    mesh.edge_metrics.j₊½.ξ̂.x₁[i, j] = ξ̂_xⱼ₊½
    mesh.edge_metrics.j₊½.ξ̂.x₂[i, j] = ξ̂_yⱼ₊½
    mesh.edge_metrics.j₊½.η̂.x₁[i, j] = η̂_xⱼ₊½
    mesh.edge_metrics.j₊½.η̂.x₂[i, j] = η̂_yⱼ₊½

    mesh.edge_metrics.j₊½.ξ.x₁[i, j] = ξ_xⱼ₊½
    mesh.edge_metrics.j₊½.ξ.x₂[i, j] = ξ_yⱼ₊½

    mesh.edge_metrics.j₊½.η.x₁[i, j] = η_xⱼ₊½
    mesh.edge_metrics.j₊½.η.x₂[i, j] = η_yⱼ₊½
  end

  return nothing
end
