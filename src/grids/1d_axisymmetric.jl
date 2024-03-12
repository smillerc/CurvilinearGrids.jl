
struct CylindricalGrid1D{CO,CE,NV,EM,CM,DL,CI,CF,DX} <: AbstractCurvilinearGrid1D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::Int
  domain_limits::DL
  iterators::CI
  snap_to_axis::Bool
  _coordinate_funcs::CF
  _∂x∂ξ::DX
end

struct SphericalGrid1D{CO,CE,NV,EM,CM,DL,CI,CF,DX} <: AbstractCurvilinearGrid1D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::Int
  domain_limits::DL
  iterators::CI
  snap_to_axis::Bool
  _coordinate_funcs::CF
  _∂x∂ξ::DX
end

AxisymmetricGrid1D = Union{CylindricalGrid1D,SphericalGrid1D}

function CylindricalGrid1D(
  r::Function, ni::Int, nhalo, snap_to_axis; T=Float64, backend=CPU()
)
  dim = 1
  check_nargs(r, dim, :r)
  test_coord_func(r, dim, :r)

  function ∂x∂ξ(ξ, t)
    return ForwardDiff.derivative(r, ξ)
  end

  nnodes = ni
  ncells = nnodes - 1
  ilo = nhalo + 1
  limits = (node=(ilo=ilo, ihi=ni + nhalo), cell=(ilo=ilo, ihi=ncells + nhalo))

  nodeCI = CartesianIndices((nnodes + 2nhalo,))
  cellCI = CartesianIndices((ncells + 2nhalo,))

  domain_iterators = get_node_cell_iterators(nodeCI, cellCI, nhalo)

  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  coordinate_funcs = (; r)
  centroids = StructArray((r=KernelAbstractions.zeros(backend, T, celldims),))
  _axisym_centroid_coordinates!(
    centroids, coordinate_funcs, domain_iterators.cell.full, nhalo
  )

  coords = StructArray((r=KernelAbstractions.zeros(backend, T, nodedims),))
  _axisym_node_coordinates!(coords, coordinate_funcs, domain_iterators.node.full, nhalo)

  node_velocities = StructArray((r=KernelAbstractions.zeros(backend, T, nodedims),))

  m = CylindricalGrid1D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes,
    limits,
    domain_iterators,
    snap_to_axis,
    coordinate_funcs,
    ∂x∂ξ,
  )

  update_metrics!(m)
  # check_for_invalid_metrics(m)
  return m
end

function SphericalGrid1D(
  r::Function, ni::Int, nhalo, snap_to_axis; T=Float64, backend=CPU()
)
  dim = 1
  check_nargs(r, dim, :r)
  test_coord_func(r, dim, :r)

  function ∂x∂ξ(ξ, t)
    return ForwardDiff.derivative(r, ξ)
  end

  nnodes = ni
  ncells = nnodes - 1
  ilo = nhalo + 1
  limits = (node=(ilo=ilo, ihi=ni + nhalo), cell=(ilo=ilo, ihi=ncells + nhalo))

  nodeCI = CartesianIndices((nnodes + 2nhalo,))
  cellCI = CartesianIndices((ncells + 2nhalo,))

  domain_iterators = get_node_cell_iterators(nodeCI, cellCI, nhalo)

  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  coordinate_funcs = (; r)
  centroids = StructArray((r=KernelAbstractions.zeros(backend, T, celldims),))
  _axisym_centroid_coordinates!(
    centroids, coordinate_funcs, domain_iterators.cell.full, nhalo
  )

  coords = StructArray((r=KernelAbstractions.zeros(backend, T, nodedims),))
  _axisym_node_coordinates!(coords, coordinate_funcs, domain_iterators.node.full, nhalo)

  node_velocities = StructArray((r=KernelAbstractions.zeros(backend, T, nodedims),))

  m = SphericalGrid1D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes,
    limits,
    domain_iterators,
    snap_to_axis,
    coordinate_funcs,
    ∂x∂ξ,
  )

  update_metrics!(m)
  # check_for_invalid_metrics(m)
  return m
end

# ------------------------------------------------------------------
# Coordinate Functions
# ------------------------------------------------------------------

function _axisym_node_coordinates!(coordinates, coordinate_functions, domain, nhalo)

  # Populate the node coordinates
  for idx in domain
    i = idx - nhalo
    coordinates.r[idx] = coordinate_functions.r(i)
  end

  return nothing
end

function _axisym_centroid_coordinates!(centroids, coordinate_functions, domain, nhalo)

  # Populate the centroid coordinates
  for idx in domain
    i = idx - nhalo + 0.5
    centroids.r[idx] = coordinate_functions.r(i)
  end

  return nothing
end