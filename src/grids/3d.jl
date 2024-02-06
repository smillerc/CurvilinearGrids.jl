
"""
CurvilinearGrid3D

# Fields
 - `x`: Node function; e.g., x(i,j,k)
 - `y`: Node function; e.g., y(i,j,k)
 - `z`: Node function; e.g., z(i,j,k)
 - `jacobian_matrix_func`: Function to compute the jacobian matrix, e.g., J(i,j,k)
 - `conserv_metric_func`: Function to compute the conservative metrics
 - `nhalo`: Number of halo cells for all dims
 - `nnodes`: Number of nodes/vertices
 - `limits`: Cell loop limits based on halo cells
"""
struct CurvilinearGrid3D{A,B,EM,CCM,F4,L,CI,DS} <: AbstractCurvilinearGrid
  coord_funcs::A
  centroids::B
  edge_metrics::EM
  cell_center_metrics::CCM
  jacobian_matrix_func::F4 # jacobian_matrix(ξ,η,ζ)
  nhalo::Int
  nnodes::NTuple{3,Int}
  domain_limits::L
  iterators::CI
  discretization_scheme::DS
end

function CurvilinearGrid3D(
  x::Function,
  y::Function,
  z::Function,
  (n_ξ, n_η, n_ζ),
  nhalo,
  cache=true,
  use_autodiff=true,
  T=Float64,
)
  dim = 3
  check_nargs(x, dim, :x)
  check_nargs(y, dim, :y)
  check_nargs(z, dim, :z)

  test_coord_func(x, dim, :x)
  test_coord_func(y, dim, :y)
  test_coord_func(z, dim, :z)

  coord(i, j, k) = @SVector [x(i, j, k), y(i, j, k), z(i, j, k)]
  function jacobian_matrix_func(i, j, k)
    return ForwardDiff.jacobian(x -> coord(x[1], x[2], x[3]), @SVector [i, j, k])
  end

  # jacobian_matrix_func = _setup_jacobian_func(x, y, z)
  # cons_metric_func = _setup_conservative_metrics_func(x, y, z)
  nnodes = (n_ξ, n_η, n_ζ)
  ncells = nnodes .- 1
  ni_cells = n_ξ - 1
  nj_cells = n_η - 1
  nk_cells = n_ζ - 1
  lo = nhalo + 1

  limits = (
    ilo=lo, ihi=ni_cells + nhalo, jlo=lo, jhi=nj_cells + nhalo, klo=lo, khi=nk_cells + nhalo
  )

  limits = (
    node=(ilo=lo, ihi=n_ξ + nhalo, jlo=lo, jhi=n_η + nhalo, klo=lo, khi=n_ζ + nhalo),
    cell=(
      ilo=lo,
      ihi=ni_cells + nhalo,
      jlo=lo,
      jhi=nj_cells + nhalo,
      klo=lo,
      khi=nk_cells + nhalo,
    ),
  )

  nodeCI = CartesianIndices(nnodes .+ 2nhalo)
  cellCI = CartesianIndices(ncells .+ 2nhalo)

  _iterators = (
    node=(
      full=nodeCI,
      domain=nodeCI[
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
      ],
      ilo_halo=nodeCI[
        begin:(begin + nhalo), (begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo)
      ],
      jlo_halo=nodeCI[
        (begin + nhalo):(end - nhalo), begin:(begin + nhalo), (begin + nhalo):(end - nhalo)
      ],
      klo_halo=nodeCI[
        (begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo), begin:(begin + nhalo)
      ],
      ihi_halo=nodeCI[(end - nhalo):end, (begin + nhalo):(end - nhalo), (end - nhalo):end],
      jhi_halo=nodeCI[
        (begin + nhalo):(end - nhalo), (end - nhalo):end, (begin + nhalo):(end - nhalo)
      ],
      khi_halo=nodeCI[
        (begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo), (end - nhalo):end
      ],
    ),
    cell=(
      full=cellCI,
      domain=cellCI[
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
      ],
      ilo_halo=cellCI[
        begin:(begin + nhalo - 1),
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
      ],
      jlo_halo=cellCI[
        (begin + nhalo):(end - nhalo),
        begin:(begin + nhalo - 1),
        (begin + nhalo):(end - nhalo),
      ],
      klo_halo=cellCI[
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
        begin:(begin + nhalo - 1),
      ],
      ihi_halo=cellCI[
        (end - nhalo + 1):end, (begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo)
      ],
      jhi_halo=cellCI[
        (begin + nhalo):(end - nhalo), (end - nhalo + 1):end, (begin + nhalo):(end - nhalo)
      ],
      khi_halo=cellCI[
        (begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo), (end - nhalo + 1):end
      ],
    ),
  )

  cell_center_metrics = (
    J=zeros(T, _iterators.cell.full.indices),
    ξx=zeros(T, _iterators.cell.full.indices),
    ξy=zeros(T, _iterators.cell.full.indices),
    ξz=zeros(T, _iterators.cell.full.indices),
    ξt=zeros(T, _iterators.cell.full.indices),
    ηx=zeros(T, _iterators.cell.full.indices),
    ηy=zeros(T, _iterators.cell.full.indices),
    ηz=zeros(T, _iterators.cell.full.indices),
    ηt=zeros(T, _iterators.cell.full.indices),
    ζx=zeros(T, _iterators.cell.full.indices),
    ζy=zeros(T, _iterators.cell.full.indices),
    ζz=zeros(T, _iterators.cell.full.indices),
    ζt=zeros(T, _iterators.cell.full.indices),
  )

  edge_metrics = (
    i₊½=(
      J=zeros(T, _iterators.cell.full.indices),
      ξ̂x=zeros(T, _iterators.cell.full.indices),
      ξ̂y=zeros(T, _iterators.cell.full.indices),
      ξ̂z=zeros(T, _iterators.cell.full.indices),
      ξ̂t=zeros(T, _iterators.cell.full.indices),
      η̂x=zeros(T, _iterators.cell.full.indices),
      η̂y=zeros(T, _iterators.cell.full.indices),
      η̂z=zeros(T, _iterators.cell.full.indices),
      η̂t=zeros(T, _iterators.cell.full.indices),
      ζ̂x=zeros(T, _iterators.cell.full.indices),
      ζ̂y=zeros(T, _iterators.cell.full.indices),
      ζ̂z=zeros(T, _iterators.cell.full.indices),
      ζ̂t=zeros(T, _iterators.cell.full.indices),
    ),
    j₊½=(
      J=zeros(T, _iterators.cell.full.indices),
      ξ̂x=zeros(T, _iterators.cell.full.indices),
      ξ̂y=zeros(T, _iterators.cell.full.indices),
      ξ̂z=zeros(T, _iterators.cell.full.indices),
      ξ̂t=zeros(T, _iterators.cell.full.indices),
      η̂x=zeros(T, _iterators.cell.full.indices),
      η̂y=zeros(T, _iterators.cell.full.indices),
      η̂z=zeros(T, _iterators.cell.full.indices),
      η̂t=zeros(T, _iterators.cell.full.indices),
      ζ̂x=zeros(T, _iterators.cell.full.indices),
      ζ̂y=zeros(T, _iterators.cell.full.indices),
      ζ̂z=zeros(T, _iterators.cell.full.indices),
      ζ̂t=zeros(T, _iterators.cell.full.indices),
    ),
    k₊½=(
      J=zeros(T, _iterators.cell.full.indices),
      ξ̂x=zeros(T, _iterators.cell.full.indices),
      ξ̂y=zeros(T, _iterators.cell.full.indices),
      ξ̂z=zeros(T, _iterators.cell.full.indices),
      ξ̂t=zeros(T, _iterators.cell.full.indices),
      η̂x=zeros(T, _iterators.cell.full.indices),
      η̂y=zeros(T, _iterators.cell.full.indices),
      η̂z=zeros(T, _iterators.cell.full.indices),
      η̂t=zeros(T, _iterators.cell.full.indices),
      ζ̂x=zeros(T, _iterators.cell.full.indices),
      ζ̂y=zeros(T, _iterators.cell.full.indices),
      ζ̂z=zeros(T, _iterators.cell.full.indices),
      ζ̂t=zeros(T, _iterators.cell.full.indices),
    ),
  )

  _centroids = (
    x=zeros(T, _iterators.cell.full.indices),
    y=zeros(T, _iterators.cell.full.indices),
    z=zeros(T, _iterators.cell.full.indices),
  )

  @inbounds for idx in _iterators.cell.full
    cell_idx = @. idx.I - nhalo + 0.5

    _centroids.x[idx] = x(cell_idx...)
    _centroids.y[idx] = y(cell_idx...)
    _centroids.z[idx] = z(cell_idx...)
  end

  discretization_scheme = MetricDiscretizationSchemes.MonotoneExplicit6thOrderDiscretization(
    _iterators.cell.full
  )
  @assert nhalo == 5

  m = CurvilinearGrid3D(
    (; x, y, z),
    _centroids,
    edge_metrics,
    cell_center_metrics,
    jacobian_matrix_func,
    nhalo,
    nnodes,
    limits,
    _iterators,
    discretization_scheme,
  )

  update_metrics!(m)

  return m
end

function update_metrics!(m::CurvilinearGrid3D)
  # Update the metrics within the non-halo region, e.g., the domain
  domain = m.iterators.cell.full

  MetricDiscretizationSchemes.update_metrics!(
    m.discretization_scheme, m.centroids, m.cell_center_metrics, m.edge_metrics, domain
  )

  return nothing
end

# function update_metrics!(m::CurvilinearGrid3D, t=0)

#   # cell metrics
#   @inbounds for idx in m.iterators.cell.full
#     cell_idx = idx.I .+ 0.5
#     @unpack J, ξ, η, ζ = metrics(m, cell_idx, t)

#     m.cell_center_metrics.ξ[idx] = ξ
#     m.cell_center_metrics.η[idx] = η
#     m.cell_center_metrics.ζ[idx] = ζ
#     m.cell_center_metrics.J[idx] = J
#   end

#   # # i₊½ conserved metrics
#   # @inbounds for idx in m.iterators.cell.full
#   #   i, j, k = idx.I .+ 0.5 # centroid index

#   #   # get the conserved metrics at (i₊½, j)
#   #   @unpack ξ̂, η̂, ζ̂ = conservative_metrics(m, (i + 1 / 2, j, k), t)

#   #   m.edge_metrics.i₊½.ξ̂[idx] = ξ̂
#   #   m.edge_metrics.i₊½.η̂[idx] = η̂
#   #   m.edge_metrics.i₊½.ζ̂[idx] = ζ̂
#   # end

#   # # j₊½ conserved metrics
#   # @inbounds for idx in m.iterators.cell.full
#   #   i, j, k = idx.I .+ 0.5 # centroid index

#   #   # get the conserved metrics at (i₊½, j)
#   #   @unpack ξ̂, η̂, ζ̂ = conservative_metrics(m, (i, j + 1 / 2, k), t)

#   #   m.edge_metrics.j₊½.ξ̂[idx] = ξ̂
#   #   m.edge_metrics.j₊½.η̂[idx] = η̂
#   #   m.edge_metrics.j₊½.ζ̂[idx] = ζ̂
#   # end

#   # # k₊½ conserved metrics
#   # @inbounds for idx in m.iterators.cell.full
#   #   i, j, k = idx.I .+ 0.5 # centroid index

#   #   # get the conserved metrics at (i₊½, j)
#   #   @unpack ξ̂, η̂, ζ̂ = conservative_metrics(m, (i, j, k + 1 / 2), t)

#   #   m.edge_metrics.k₊½.ξ̂[idx] = ξ̂
#   #   m.edge_metrics.k₊½.η̂[idx] = η̂
#   #   m.edge_metrics.k₊½.ζ̂[idx] = ζ̂
#   # end

#   return nothing
# end

# function update_metrics_meg!(m::CurvilinearGrid3D)
#   cswh = cellsize_withhalo(m)
#   meg6_metrics = MEG6Scheme(cswh)

#   # centroids
#   xc = zeros(cswh)
#   yc = zeros(cswh)
#   zc = zeros(cswh)
#   for idx in CartesianIndices(xc)
#     i, j, k = idx.I
#     xc[idx] = m.x_func(i - m.nhalo + 0.5, j - m.nhalo + 0.5, k - m.nhalo + 0.5) # x
#     yc[idx] = m.y_func(i - m.nhalo + 0.5, j - m.nhalo + 0.5, k - m.nhalo + 0.5) # y
#     zc[idx] = m.z_func(i - m.nhalo + 0.5, j - m.nhalo + 0.5, k - m.nhalo + 0.5) # y
#   end

#   domain = m.iterators.cell.domain

#   MetricDiscretizationSchemes.update_metrics!(meg6_metrics, xc, yc, zc, domain)

#   _J = J(meg6_metrics)
#   _ξx = ξx(meg6_metrics)
#   _ξy = ξy(meg6_metrics)
#   _ξz = ξz(meg6_metrics)
#   _ηx = ηx(meg6_metrics)
#   _ηy = ηy(meg6_metrics)
#   _ηz = ηz(meg6_metrics)
#   _ζx = ζx(meg6_metrics)
#   _ζy = ζy(meg6_metrics)
#   _ζz = ζz(meg6_metrics)

#   # cell centroid metrics
#   for idx in domain
#     _metric = (
#       J=_J[idx],
#       ξx=_ξx[idx],
#       ξy=_ξy[idx],
#       ξz=_ξz[idx],
#       ηx=_ηx[idx],
#       ηy=_ηy[idx],
#       ηz=_ηz[idx],
#       ζx=_ζx[idx],
#       ζy=_ζy[idx],
#       ζz=_ζz[idx],
#       ξt=zero(Float64),
#       ηt=zero(Float64),
#       ζt=zero(Float64),
#     )

#     m.cell_center_metrics[idx] = _metric
#   end

#   _ξ̂xᵢ₊½ = ξ̂xᵢ₊½(meg6_metrics)
#   _ξ̂yᵢ₊½ = ξ̂yᵢ₊½(meg6_metrics)
#   _ξ̂zᵢ₊½ = ξ̂zᵢ₊½(meg6_metrics)
#   _η̂xᵢ₊½ = η̂xᵢ₊½(meg6_metrics)
#   _η̂yᵢ₊½ = η̂yᵢ₊½(meg6_metrics)
#   _η̂zᵢ₊½ = η̂zᵢ₊½(meg6_metrics)
#   _ζ̂xᵢ₊½ = ζ̂xᵢ₊½(meg6_metrics)
#   _ζ̂yᵢ₊½ = ζ̂yᵢ₊½(meg6_metrics)
#   _ζ̂zᵢ₊½ = ζ̂zᵢ₊½(meg6_metrics)

#   _ξ̂xⱼ₊½ = ξ̂xⱼ₊½(meg6_metrics)
#   _ξ̂yⱼ₊½ = ξ̂yⱼ₊½(meg6_metrics)
#   _ξ̂zⱼ₊½ = ξ̂zⱼ₊½(meg6_metrics)
#   _η̂xⱼ₊½ = η̂xⱼ₊½(meg6_metrics)
#   _η̂yⱼ₊½ = η̂yⱼ₊½(meg6_metrics)
#   _η̂zⱼ₊½ = η̂zⱼ₊½(meg6_metrics)
#   _ζ̂xⱼ₊½ = ζ̂xⱼ₊½(meg6_metrics)
#   _ζ̂yⱼ₊½ = ζ̂yⱼ₊½(meg6_metrics)
#   _ζ̂zⱼ₊½ = ζ̂zⱼ₊½(meg6_metrics)

#   _ξ̂xₖ₊½ = ξ̂xₖ₊½(meg6_metrics)
#   _ξ̂yₖ₊½ = ξ̂yₖ₊½(meg6_metrics)
#   _ξ̂zₖ₊½ = ξ̂zₖ₊½(meg6_metrics)
#   _η̂xₖ₊½ = η̂xₖ₊½(meg6_metrics)
#   _η̂yₖ₊½ = η̂yₖ₊½(meg6_metrics)
#   _η̂zₖ₊½ = η̂zₖ₊½(meg6_metrics)
#   _ζ̂xₖ₊½ = ζ̂xₖ₊½(meg6_metrics)
#   _ζ̂yₖ₊½ = ζ̂yₖ₊½(meg6_metrics)
#   _ζ̂zₖ₊½ = ζ̂zₖ₊½(meg6_metrics)

#   i₊½CI = MetricDiscretizationSchemes.expand_lower(domain, 1, 1)
#   j₊½CI = MetricDiscretizationSchemes.expand_lower(domain, 2, 1)
#   k₊½CI = MetricDiscretizationSchemes.expand_lower(domain, 3, 1)

#   @show domain
#   @show i₊½CI
#   @show j₊½CI
#   @show k₊½CI
#   for idx in i₊½CI
#     i₊½_metric = (
#       ξ̂x=_ξ̂xᵢ₊½[idx],
#       ξ̂y=_ξ̂yᵢ₊½[idx],
#       ξ̂z=_ξ̂zᵢ₊½[idx],
#       η̂x=_η̂xᵢ₊½[idx],
#       η̂y=_η̂yᵢ₊½[idx],
#       η̂z=_η̂zᵢ₊½[idx],
#       ζ̂x=_ζ̂xᵢ₊½[idx],
#       ζ̂y=_ζ̂yᵢ₊½[idx],
#       ζ̂z=_ζ̂zᵢ₊½[idx],
#       ξt=0.0,
#       ηt=0.0,
#       ζt=0.0,
#     )

#     m.edge_metrics.i₊½[idx] = i₊½_metric
#   end

#   for idx in j₊½CI
#     j₊½_metric = (
#       ξ̂x=_ξ̂xⱼ₊½[idx],
#       ξ̂y=_ξ̂yⱼ₊½[idx],
#       ξ̂z=_ξ̂zⱼ₊½[idx],
#       η̂x=_η̂xⱼ₊½[idx],
#       η̂y=_η̂yⱼ₊½[idx],
#       η̂z=_η̂zⱼ₊½[idx],
#       ζ̂x=_ζ̂xⱼ₊½[idx],
#       ζ̂y=_ζ̂yⱼ₊½[idx],
#       ζ̂z=_ζ̂zⱼ₊½[idx],
#       ξt=0.0,
#       ηt=0.0,
#       ζt=0.0,
#     )

#     m.edge_metrics.j₊½[idx] = j₊½_metric
#   end

#   for idx in k₊½CI
#     k₊½_metric = (
#       ξ̂x=_ξ̂xₖ₊½[idx],
#       ξ̂y=_ξ̂yₖ₊½[idx],
#       ξ̂z=_ξ̂zₖ₊½[idx],
#       η̂x=_η̂xₖ₊½[idx],
#       η̂y=_η̂yₖ₊½[idx],
#       η̂z=_η̂zₖ₊½[idx],
#       ζ̂x=_ζ̂xₖ₊½[idx],
#       ζ̂y=_ζ̂yₖ₊½[idx],
#       ζ̂z=_ζ̂zₖ₊½[idx],
#       ξt=0.0,
#       ηt=0.0,
#       ζt=0.0,
#     )

#     m.edge_metrics.k₊½[idx] = k₊½_metric
#   end

#   return nothing
# end

function sanity_checks(m)
  domain = m.iterators.cell.domain

  i₊½CI = MetricDiscretizationSchemes.expand_lower(domain, 1, 1)
  j₊½CI = MetricDiscretizationSchemes.expand_lower(domain, 2, 1)
  k₊½CI = MetricDiscretizationSchemes.expand_lower(domain, 3, 1)

  # cell centroid metrics
  for idx in domain
    for v in m.cell_center_metrics[idx]
      if !(isfinite(v))
        @show idx
        @show m.cell_center_metrics[idx]
        error("Invalid grid metrics")
      end
    end
  end

  for (em, dom) in zip(m.edge_metrics, (i₊½CI, j₊½CI, k₊½CI))
    @show dom
    for idx in dom
      @show em[idx]
      for v in em[idx]
        if !(isfinite(v))
          @show idx
          error("Invalid grid metrics")
        end
      end
    end
  end

  return nothing
end

function update_metrics_old!(m::CurvilinearGrid3D)
  function _get_metrics(m, (i, j, k))
    _jacobian_matrix = m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo)

    inv_jacobian_matrix = inv(_jacobian_matrix)

    return (
      J=det(_jacobian_matrix),
      ξx=inv_jacobian_matrix[1, 1],
      ξy=inv_jacobian_matrix[1, 2],
      ξz=inv_jacobian_matrix[1, 3],
      ηx=inv_jacobian_matrix[2, 1],
      ηy=inv_jacobian_matrix[2, 2],
      ηz=inv_jacobian_matrix[2, 3],
      ζx=inv_jacobian_matrix[3, 1],
      ζy=inv_jacobian_matrix[3, 2],
      ζz=inv_jacobian_matrix[3, 3],
      ξt=zero(eltype(_jacobian_matrix)),
      ηt=zero(eltype(_jacobian_matrix)),
      ζt=zero(eltype(_jacobian_matrix)),
    )
  end

  # cell centroid metrics
  for idx in m.iterators.cell.domain
    node_idx = idx.I .- m.nhalo
    cell_idx = node_idx .+ 0.5
    @unpack J, ξx, ξy, ξz, ηx, ηy, ηz, ζx, ζy, ζz = _get_metrics(m, cell_idx)

    _metric = (
      J=J,
      ξx=ξx,
      ξ̂x=ξx / J,
      ξy=ξy,
      ξ̂y=ξy / J,
      ξz=ξz,
      ξ̂z=ξz / J,
      ηx=ηx,
      η̂x=ηx / J,
      ηy=ηy,
      η̂y=ηy / J,
      ηz=ηz,
      η̂z=ηz / J,
      ζx=ζx,
      ζ̂x=ζx / J,
      ζy=ζy,
      ζ̂y=ζy / J,
      ζz=ζz,
      ζ̂z=ζz / J,
      ξt=zero(J),
      ηt=zero(J),
      ζt=zero(J),
    )

    m.cell_center_metrics[idx] = _metric
  end

  ilo_n, ihi_n, jlo_n, jhi_n, klo_n, khi_n = m.domain_limits.node

  i₊½CI = CartesianIndices((ilo_n:ihi_n, jlo_n:(jhi_n - 1), klo_n:(khi_n - 1)))
  j₊½CI = CartesianIndices((ilo_n:(ihi_n - 1), jlo_n:jhi_n, klo_n:(khi_n - 1)))
  k₊½CI = CartesianIndices((ilo_n:(ihi_n - 1), jlo_n:(jhi_n - 1), klo_n:khi_n))

  # The grid x,y functions define the NODE position
  # of the entire mesh, and do not know about halo cells
  # The `- m.nhalo` accounts for this and the `+ 1/2`
  # is so we get the middle of the edge. 

  #             o-------------o
  #             |             |
  #             |             |
  #     i₊½     X    (i,j)    |   
  #   for the   |             |   
  #   cell to   |             |   
  #   the left  o------X------o   
  #                   j₊½ for the cell below
  # the lower left corner node is (i,j) as well,
  # so this is why the ₊½ node indices seem backwards...
  i₊½_offset = (0.0, 0.5, 0.5)
  i₊½_cell_offset = (-1, 0, 0)
  for idx in i₊½CI
    node_idx = idx.I
    i₊½_node_idx = node_idx .+ i₊½_offset .- m.nhalo
    i₊½_cell_idx = CartesianIndex(node_idx .+ i₊½_cell_offset)

    J_i₊½, ξx_i₊½, ξy_i₊½, ξz_i₊½, ηx_i₊½, ηy_i₊½, ηz_i₊½, ζx_i₊½, ζy_i₊½, ζz_i₊½ = _get_metrics(
      m, i₊½_node_idx
    )

    i₊½_metric = (
      J=J_i₊½,
      ξx=ξx_i₊½,
      ξ̂x=ξx_i₊½ / J_i₊½,
      ξy=ξy_i₊½,
      ξ̂y=ξy_i₊½ / J_i₊½,
      ξz=ξz_i₊½,
      ξ̂z=ξz_i₊½ / J_i₊½,
      ηx=ηx_i₊½,
      η̂x=ηx_i₊½ / J_i₊½,
      ηy=ηy_i₊½,
      η̂y=ηy_i₊½ / J_i₊½,
      ηz=ηz_i₊½,
      η̂z=ηz_i₊½ / J_i₊½,
      ζx=ζx_i₊½,
      ζ̂x=ζx_i₊½ / J_i₊½,
      ζy=ζy_i₊½,
      ζ̂y=ζy_i₊½ / J_i₊½,
      ζz=ζz_i₊½,
      ζ̂z=ζz_i₊½ / J_i₊½,
      ξt=zero(J_i₊½),
      ηt=zero(J_i₊½),
      ζt=zero(J_i₊½),
    )

    m.edge_metrics.i₊½[i₊½_cell_idx] = i₊½_metric
  end

  j₊½_offset = (0.5, 0.0, 0.5)
  j₊½_cell_offset = (0, -1, 0)
  for idx in j₊½CI
    node_idx = idx.I
    j₊½_node_idx = node_idx .+ j₊½_offset .- m.nhalo
    j₊½_cell_idx = CartesianIndex(node_idx .+ j₊½_cell_offset)
    J_j₊½, ξx_j₊½, ξy_j₊½, ξz_j₊½, ηx_j₊½, ηy_j₊½, ηz_j₊½, ζx_j₊½, ζy_j₊½, ζz_j₊½ = _get_metrics(
      m, j₊½_node_idx
    )

    j₊½_metric = (
      J=J_j₊½,
      ξx=ξx_j₊½,
      ξ̂x=ξx_j₊½ / J_j₊½,
      ξy=ξy_j₊½,
      ξ̂y=ξy_j₊½ / J_j₊½,
      ξz=ξz_j₊½,
      ξ̂z=ξz_j₊½ / J_j₊½,
      ηx=ηx_j₊½,
      η̂x=ηx_j₊½ / J_j₊½,
      ηy=ηy_j₊½,
      η̂y=ηy_j₊½ / J_j₊½,
      ηz=ηz_j₊½,
      η̂z=ηz_j₊½ / J_j₊½,
      ζx=ζx_j₊½,
      ζ̂x=ζx_j₊½ / J_j₊½,
      ζy=ζy_j₊½,
      ζ̂y=ζy_j₊½ / J_j₊½,
      ζz=ζz_j₊½,
      ζ̂z=ζz_j₊½ / J_j₊½,
      ξt=zero(J_j₊½),
      ηt=zero(J_j₊½),
      ζt=zero(J_j₊½),
    )

    m.edge_metrics.j₊½[j₊½_cell_idx] = j₊½_metric
  end

  k₊½_offset = (0.5, 0.5, 0.0)
  k₊½_cell_offset = (0, 0, -1)
  for idx in k₊½CI
    node_idx = idx.I
    k₊½_node_idx = node_idx .+ k₊½_offset .- m.nhalo
    k₊½_cell_idx = CartesianIndex(node_idx .+ k₊½_cell_offset)
    J_k₊½, ξx_k₊½, ξy_k₊½, ξz_k₊½, ηx_k₊½, ηy_k₊½, ηz_k₊½, ζx_k₊½, ζy_k₊½, ζz_k₊½ = _get_metrics(
      m, k₊½_node_idx
    )

    k₊½_metric = (
      J=J_k₊½,
      ξx=ξx_k₊½,
      ξ̂x=ξx_k₊½ / J_k₊½,
      ξy=ξy_k₊½,
      ξ̂y=ξy_k₊½ / J_k₊½,
      ξz=ξz_k₊½,
      ξ̂z=ξz_k₊½ / J_k₊½,
      ηx=ηx_k₊½,
      η̂x=ηx_k₊½ / J_k₊½,
      ηy=ηy_k₊½,
      η̂y=ηy_k₊½ / J_k₊½,
      ηz=ηz_k₊½,
      η̂z=ηz_k₊½ / J_k₊½,
      ζx=ζx_k₊½,
      ζ̂x=ζx_k₊½ / J_k₊½,
      ζy=ζy_k₊½,
      ζ̂y=ζy_k₊½ / J_k₊½,
      ζz=ζz_k₊½,
      ζ̂z=ζz_k₊½ / J_k₊½,
      ξt=zero(J_k₊½),
      ηt=zero(J_k₊½),
      ζt=zero(J_k₊½),
    )

    m.edge_metrics.k₊½[k₊½_cell_idx] = k₊½_metric
  end

  return nothing
end

@inline conservative_metrics(m::CurvilinearGrid3D, idx) = conservative_metrics(m, idx, 0)

@inline function conservative_metrics(
  m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real
)
  # hidx = (i - m.nhalo, j - m.nhalo, k - m.nhalo)

  _jacobian_matrix = checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo))
  T = eltype(_jacobian_matrix)
  J = det(_jacobian_matrix)

  inv_jacobian_matrix = inv(_jacobian_matrix)

  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ξz = inv_jacobian_matrix[1, 3]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]
  ηz = inv_jacobian_matrix[2, 3]
  ζx = inv_jacobian_matrix[3, 1]
  ζy = inv_jacobian_matrix[3, 2]
  ζz = inv_jacobian_matrix[3, 3]

  xξ = _jacobian_matrix[1, 1]
  yξ = _jacobian_matrix[2, 1]
  zξ = _jacobian_matrix[3, 1]
  xη = _jacobian_matrix[1, 2]
  yη = _jacobian_matrix[2, 2]
  zη = _jacobian_matrix[3, 2]
  xζ = _jacobian_matrix[1, 3]
  yζ = _jacobian_matrix[2, 3]
  zζ = _jacobian_matrix[3, 3]

  ξt = zero(T)
  ηt = zero(T)
  ζt = zero(T)

  # ξ̂ = Metric3D(ξx * J, ξy * J, ξz * J, ξt * J)
  # η̂ = Metric3D(ηx * J, ηy * J, ηz * J, ηt * J)
  # ζ̂ = Metric3D(ζx * J, ζy * J, ζz * J, ζt * J)

  M1, M2, M3 = m.conserv_metric_func(i - m.nhalo, j - m.nhalo, k - m.nhalo) # get the matrices

  ∂ξ = 1
  ∂η = 2
  ∂ζ = 3
  yξzη = M1[∂ξ, ∂η] # ∂(z ∂y/∂ξ)/∂η
  yξzζ = M1[∂ξ, ∂ζ] # ∂(z ∂y/∂ξ)/∂ζ

  yηzξ = M1[∂η, ∂ξ] # ∂(z ∂y/∂η)/∂ξ
  yηzζ = M1[∂η, ∂ζ] # ∂(z ∂y/∂η)/∂ζ

  yζzξ = M1[∂ζ, ∂ξ] # ∂(z ∂y/∂ζ)/∂ξ
  yζzη = M1[∂ζ, ∂η] # ∂(z ∂y/∂ζ)/∂η

  zξxη = M2[∂ξ, ∂η] # ∂(x ∂z/∂ξ)/∂η
  zξxζ = M2[∂ξ, ∂ζ] # ∂(x ∂z/∂ξ)/∂ζ
  zηxξ = M2[∂η, ∂ξ] # ∂(x ∂z/∂η)/∂ξ
  zηxζ = M2[∂η, ∂ζ] # ∂(x ∂z/∂η)/∂ζ
  zζxξ = M2[∂ζ, ∂ξ] # ∂(x ∂z/∂ζ)/∂ξ
  zζxη = M2[∂ζ, ∂η] # ∂(x ∂z/∂ζ)/∂η

  xξyη = M3[∂ξ, ∂η] # ∂(y ∂x/∂ξ)/∂η
  xξyζ = M3[∂ξ, ∂ζ] # ∂(y ∂x/∂ξ)/∂ζ
  xηyξ = M3[∂η, ∂ξ] # ∂(y ∂x/∂η)/∂ξ
  xηyζ = M3[∂η, ∂ζ] # ∂(y ∂x/∂η)/∂ζ
  xζyξ = M3[∂ζ, ∂ξ] # ∂(y ∂x/∂ζ)/∂ξ
  xζyη = M3[∂ζ, ∂η] # ∂(y ∂x/∂ζ)/∂η

  ξ̂x = yηzζ - yζzη

  ξ̂x1 = yη * zζ - yζ * zη
  # @show abs(ξ̂x1 - ξ̂x)

  ξ̂y = zηxζ - zζxη
  ξ̂z = xηyζ - xζyη

  η̂x = yζzξ - yξzζ
  η̂y = zζxξ - zξxζ
  η̂z = xζyξ - xξyζ

  ζ̂x = yξzη - yηzξ
  ζ̂y = zξxη - zηxξ
  ζ̂z = xξyη - xηyξ
  # @show abs(ξ̂x - ξx * J), abs(ξ̂y - ξy * J), abs(ξ̂z - ξz * J)
  # @show abs(η̂x - ηx * J), abs(η̂y - ηy * J), abs(η̂z - ηz * J)

  ξ̂ = Metric3D(ξ̂x, ξ̂y, ξ̂z, ξt)
  η̂ = Metric3D(η̂x, η̂y, η̂z, ηt)
  ζ̂ = Metric3D(ζ̂x, ζ̂y, ζ̂z, ζt)
  # @show hidx, ξ̂

  return (; ξ̂, η̂, ζ̂)
end

# @inline function conservative_metrics(
#   m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, (vx, vy, vz)
# )
#   static = conservative_metrics(m, (i, j, k))
#   @unpack ξ̂x, ξ̂y, ξ̂z, η̂x, η̂y, η̂z, ζ̂x, ζ̂y, ζ̂z = static

#   return merge(
#     static,
#     (
#       ξt=-(vx * ξ̂x + vy * ξ̂y + vz * ξ̂z), # dynamic / moving mesh terms
#       ηt=-(vx * η̂x + vy * η̂y + vz * η̂z), # dynamic / moving mesh terms
#       ζt=-(vx * ζ̂x + vy * ζ̂y + vz * ζ̂z), # dynamic / moving mesh terms
#     ),
#   )
# end

# Get the grid metrics for a static grid
@inline function metrics(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real)
  _jacobian_matrix = checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo))
  J = det(_jacobian_matrix)

  inv_jacobian_matrix = inv(_jacobian_matrix)

  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ξz = inv_jacobian_matrix[1, 3]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]
  ηz = inv_jacobian_matrix[2, 3]
  ζx = inv_jacobian_matrix[3, 1]
  ζy = inv_jacobian_matrix[3, 2]
  ζz = inv_jacobian_matrix[3, 3]

  ξt = zero(eltype(_jacobian_matrix))
  ηt = zero(eltype(_jacobian_matrix))
  ζt = zero(eltype(_jacobian_matrix))

  ξ = Metric3D(ξx, ξy, ξz, ξt)
  η = Metric3D(ηx, ηy, ηz, ηt)
  ζ = Metric3D(ζx, ζy, ζz, ζt)

  return (; ξ, η, ζ, J)
end

# @inline function metrics(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, (vx, vy, vz))
#   static = metrics(m, (i, j, k))
#   @unpack ξx, ξy, ξz, ηx, ηy, ηz, ζx, ζy, ζz = static

#   return merge(
#     static,
#     (
#       ξt=-(vx * ξx + vy * ξy + vz * ξz), # dynamic / moving mesh terms
#       ηt=-(vx * ηx + vy * ηy + vz * ηz), # dynamic / moving mesh terms
#       ζt=-(vx * ζx + vy * ζy + vz * ζz), # dynamic / moving mesh terms
#     ),
#   )
# end

function jacobian_matrix(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real})
  # return checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo))
  return m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo)
end

function jacobian(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real})
  return det(jacobian_matrix(m, (i, j, k)))
end

function _setup_jacobian_func(x, y, z)
  _x((i, j, k)) = x(i, j, k)
  _y((i, j, k)) = y(i, j, k)
  _z((i, j, k)) = z(i, j, k)

  ∇x(ξ, η, ζ) = ForwardDiff.gradient(_x, @SVector [ξ, η, ζ]) # ∂x∂ξ, ∂x∂η, ∂x∂ζ
  ∇y(ξ, η, ζ) = ForwardDiff.gradient(_y, @SVector [ξ, η, ζ]) # ∂z∂ξ, ∂z∂η, ∂z∂ζ
  ∇z(ξ, η, ζ) = ForwardDiff.gradient(_z, @SVector [ξ, η, ζ]) # ∂y∂ξ, ∂y∂η, ∂y∂ζ

  function jacobian_matrix(ξ, η, ζ)
    return SMatrix{3,3,Float64}(∇x(ξ, η, ζ)..., ∇y(ξ, η, ζ)..., ∇z(ξ, η, ζ)...)
  end

  return jacobian_matrix
end

# Get the grid metrics for a dynamic grid -- grid velocities (vx,vy,vz) must be provided
function _get_metrics(_jacobian_matrix::SMatrix{3,3,T}, (vx, vy, vz)) where {T}
  inv_jacobian_matrix = inv(_jacobian_matrix)
  J⁻¹ = det(inv_jacobian_matrix)

  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ξz = inv_jacobian_matrix[1, 3]

  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]
  ηz = inv_jacobian_matrix[2, 3]

  ζx = inv_jacobian_matrix[3, 1]
  ζy = inv_jacobian_matrix[3, 2]
  ζz = inv_jacobian_matrix[3, 3]

  # temporal metrics
  ξt = -(vx * ξx + vy * ξy + vz * ξz)
  ηt = -(vx * ηx + vy * ηy + vz * ηz)
  ζt = -(vx * ζx + vy * ζy + vz * ζz)

  ξ̂x = ξx * J⁻¹
  ξ̂y = ξy * J⁻¹
  ξ̂z = ξz * J⁻¹
  η̂x = ηx * J⁻¹
  η̂y = ηy * J⁻¹
  η̂z = ηz * J⁻¹
  ζ̂x = ζx * J⁻¹
  ζ̂y = ζy * J⁻¹
  ζ̂z = ζz * J⁻¹

  return (
    ξ̂x=ξ̂x,
    ξ̂y=ξ̂y,
    ξ̂z=ξ̂z,
    η̂x=η̂x,
    η̂y=η̂y,
    η̂z=η̂z,
    ζ̂x=ζ̂x,
    ζ̂y=ζ̂y,
    ζ̂z=ζ̂z,
    ξt=ξt,
    ηt=ηt,
    ζt=ζt,
  )
end

function _setup_conservative_metrics_func(x, y, z)
  _x((i, j, k)) = x(i, j, k)
  _y((i, j, k)) = y(i, j, k)
  _z((i, j, k)) = z(i, j, k)

  gradx_y((i, j, k)) = ForwardDiff.gradient(_x, @SVector [i, j, k]) * y(i, j, k)
  grady_z((i, j, k)) = ForwardDiff.gradient(_y, @SVector [i, j, k]) * z(i, j, k)
  gradz_x((i, j, k)) = ForwardDiff.gradient(_z, @SVector [i, j, k]) * x(i, j, k)

  grady_z_jacobian(i, j, k) = ForwardDiff.jacobian(grady_z, @SVector [i, j, k])
  gradz_x_jacobian(i, j, k) = ForwardDiff.jacobian(gradz_x, @SVector [i, j, k])
  gradx_y_jacobian(i, j, k) = ForwardDiff.jacobian(gradx_y, @SVector [i, j, k])

  # Get all the matrices at once
  function conserv_metric_matricies(ξ, η, ζ)
    # checkeps just zeros out terms that are less than ϵ
    return (
      checkeps(grady_z_jacobian(ξ, η, ζ)),
      checkeps(gradz_x_jacobian(ξ, η, ζ)),
      checkeps(gradx_y_jacobian(ξ, η, ζ)),
    )
  end

  # M1 = ∇y_z_jacobian(ξ, η, ζ)
  # M1 = [(z y_ξ)_ξ (z y_ξ)_η (z y_ξ)_ζ
  #       (z y_η)_ξ (z y_η)_η (z y_η)_ζ
  #       (z y_ζ)_ξ (z y_ζ)_η (z y_ζ)_ζ ]

  # M2 = ∇z_x_jacobian(ξ, η, ζ)

  # M2 = [(x z_ξ)_ξ (x z_ξ)_η (x z_ξ)_ζ
  #       (x z_η)_ξ (x z_η)_η (x z_η)_ζ
  #       (x z_ζ)_ξ (x z_ζ)_η (x z_ζ)_ζ ]

  # M3 = ∇x_y_jacobian(ξ, η, ζ)
  # M3 = [(y x_ξ)_ξ (y x_ξ)_η (y x_ξ)_ζ
  #       (y x_η)_ξ (y x_η)_η (y x_η)_ζ
  #       (y x_ζ)_ξ (y x_ζ)_η (y x_ζ)_ζ ]

  return conserv_metric_matricies
end

"""
    coords(mesh::CurvilinearGrid3D, T=Float64) -> Array{Real}

Return the array of coordinate points, indexed as `[xyz,i,j,k]`.
This does _not_ include halo regions since the geometry can be undefined.
"""
function coords(m::CurvilinearGrid3D, T=Float64)
  dims = m.nnodes
  xyz = zeros(T, 3, dims...)
  @inbounds for I in CartesianIndices(dims)
    i, j, k = Tuple(I)
    xyz[1, i, j, k] = m.coord_funcs.x(i, j, k)
    xyz[2, i, j, k] = m.coord_funcs.y(i, j, k)
    xyz[3, i, j, k] = m.coord_funcs.z(i, j, k)
  end

  return xyz
end

"""
    centroids(m::CurvilinearGrid3D, T=Float64) -> Array{Real}

Return the array of coordinate points, indexed as `[xyz,i,j,k]`.
This does _not_ include halo regions since the geometry can be undefined.
"""
function centroids(m::CurvilinearGrid3D, T=Float64)
  dims = (m.nnodes .- 1)
  xyz = zeros(T, 3, dims...)
  @inbounds for I in CartesianIndices(dims)
    i, j, k = Tuple(I)
    xyz[1, i, j, k] = m.x(i + 0.5, j + 0.5, k + 0.5)
    xyz[2, i, j, k] = m.y(i + 0.5, j + 0.5, k + 0.5)
    xyz[3, i, j, k] = m.z(i + 0.5, j + 0.5, k + 0.5)
  end

  return xyz
end
