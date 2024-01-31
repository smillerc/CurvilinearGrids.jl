
"""
CurvilinearGrid2D

# Fields
 - `x`: Node function; e.g., x(i,j)
 - `y`: Node function; e.g., y(i,j)
 - `jacobian_matrix_func`: jacobian matrix, e.g., J(i,j)
 - `nhalo`: Number of halo cells for all dims
 - `nnodes`: Number of nodes/vertices
 - `limits`: Cell loop limits based on halo cells
"""
struct CurvilinearGrid2D{T1,T2,T3,EM,CCM,L,CI,AA} <: AbstractCurvilinearGrid
  x_func::T1 # functional form x(i,j)
  y_func::T2 # functional form j(i,j)
  x_coord::AA # stored x-coordinate location
  y_coord::AA # stored y-coordinate location
  edge_metrics::EM
  cell_center_metrics::CCM
  jacobian_matrix_func::T3
  nhalo::Int
  nnodes::NTuple{2,Int}
  domain_limits::L
  iterators::CI
  use_autodiff::Bool
end

function CurvilinearGrid2D(
  x::Function, y::Function, (n_ξ, n_η), nhalo, cache=true, use_autodiff=true, T=Float64
)
  dim = 2
  check_nargs(x, dim, :x)
  check_nargs(y, dim, :y)
  test_coord_func(x, dim, :x)
  test_coord_func(y, dim, :y)

  xy(i, j) = @SVector [x(i, j), y(i, j)]
  jacobian_matrix_func(i, j) = ForwardDiff.jacobian(x -> xy(x[1], x[2]), @SVector [i, j])

  # jacobian_matrix_func = _setup_jacobian_func(x, y)
  nnodes = (n_ξ, n_η)
  ncells = nnodes .- 1
  ni_cells, nj_cells = ncells
  lo = nhalo + 1
  limits = (
    node=(ilo=lo, ihi=n_ξ + nhalo, jlo=lo, jhi=n_η + nhalo),
    cell=(ilo=lo, ihi=ni_cells + nhalo, jlo=lo, jhi=nj_cells + nhalo),
  )

  nodeCI = CartesianIndices(nnodes .+ 2nhalo)
  cellCI = CartesianIndices(ncells .+ 2nhalo)
  _iterators = (
    node=(
      full=nodeCI,
      domain=nodeCI[(begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo)],
      ilo_halo=nodeCI[begin:(begin + nhalo), (begin + nhalo):(end - nhalo)],
      ihi_halo=nodeCI[(end - nhalo):end, (begin + nhalo):(end - nhalo)],
      jlo_halo=nodeCI[(begin + nhalo):(end - nhalo), begin:(begin + nhalo)],
      jhi_halo=nodeCI[(begin + nhalo):(end - nhalo), (end - nhalo):end],
    ),
    cell=(
      full=cellCI,
      domain=cellCI[(begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo)],
      ilo_halo=cellCI[begin:(begin + nhalo - 1), (begin + nhalo):(end - nhalo)],
      jlo_halo=cellCI[(begin + nhalo):(end - nhalo), begin:(begin + nhalo - 1)],
      ihi_halo=cellCI[(end - nhalo + 1):end, (begin + nhalo):(end - nhalo)],
      jhi_halo=cellCI[(begin + nhalo):(end - nhalo), (end - nhalo + 1):end],
    ),
  )

  if cache
    _edge_metric = (
      # J=zero(T),
      # ξx=zero(T),
      ξ̂x=zero(T),
      # ξy=zero(T),
      ξ̂y=zero(T),
      # ηx=zero(T),
      η̂x=zero(T),
      # ηy=zero(T),
      η̂y=zero(T),
      ξt=zero(T),
      ηt=zero(T),
    )

    _cell_metric = (
      J=zero(T),
      ξx=zero(T),
      # ξ̂x=zero(T),
      ξy=zero(T),
      # ξ̂y=zero(T),
      ηx=zero(T),
      # η̂x=zero(T),
      ηy=zero(T),
      # η̂y=zero(T),
      ξt=zero(T),
      ηt=zero(T),
    )

    EM = typeof(_edge_metric)
    CM = typeof(_cell_metric)
    cell_center_metrics = Matrix{CM}(undef, ncells .+ 2nhalo)

    edge_metrics = (
      i₊½=Matrix{EM}(undef, ncells .+ 2nhalo), j₊½=Matrix{EM}(undef, ncells .+ 2nhalo)
    )

    xcoords = zeros(T, nnodes .+ 2nhalo)
    ycoords = zeros(T, nnodes .+ 2nhalo)

    @inbounds for j in axes(xcoords, 2)
      for i in axes(xcoords, 1)
        xcoords[i, j] = x(i - nhalo, j - nhalo)
        ycoords[i, j] = y(i - nhalo, j - nhalo)
      end
    end

  else
    edge_metrics = nothing
    cell_center_metrics = nothing
    xcoords = nothing
    ycoords = nothing
  end

  m = CurvilinearGrid2D(
    x,
    y,
    xcoords,
    ycoords,
    edge_metrics,
    cell_center_metrics,
    jacobian_matrix_func,
    nhalo,
    nnodes,
    limits,
    _iterators,
    use_autodiff,
  )

  # update_metrics_meg!(m)
  update_metrics_old!(m)

  return m
end

function update_metrics_meg!(m::CurvilinearGrid2D)
  cswh = cellsize_withhalo(m)
  meg6_metrics = MEG6Scheme(cswh)

  # centroids
  xc = zeros(cswh)
  yc = zeros(cswh)
  for idx in CartesianIndices(xc)
    i, j = idx.I
    xc[idx] = m.x_func(i - m.nhalo + 0.5, j - m.nhalo + 0.5) # x
    yc[idx] = m.y_func(i - m.nhalo + 0.5, j - m.nhalo + 0.5) # y
  end

  domain = m.iterators.cell.domain

  MetricDiscretizationSchemes.update_metrics!(meg6_metrics, xc, yc, domain)

  # cell centroid metrics
  _J = J(meg6_metrics)
  _ξx = ξx(meg6_metrics)
  _ξy = ξy(meg6_metrics)
  _ηx = ηx(meg6_metrics)
  _ηy = ηy(meg6_metrics)
  for idx in m.iterators.cell.domain
    _metric = (
      J=_J[idx], #
      ξx=_ξx[idx], #
      ξy=_ξy[idx], #
      ηx=_ηx[idx], #
      ηy=_ηy[idx], #
      ξt=zero(Float64),
      ηt=zero(Float64),
    )

    m.cell_center_metrics[idx] = _metric
  end

  ilo, ihi, jlo, jhi = m.domain_limits.cell

  j₊½CI = CartesianIndices((ilo:ihi, (jlo - 1):jhi))
  i₊½CI = CartesianIndices(((ilo - 1):ihi, jlo:jhi))

  _ξ̂xᵢ₊½ = ξ̂xᵢ₊½(meg6_metrics)
  _η̂xᵢ₊½ = η̂xᵢ₊½(meg6_metrics)
  _ξ̂yᵢ₊½ = ξ̂yᵢ₊½(meg6_metrics)
  _η̂yᵢ₊½ = η̂yᵢ₊½(meg6_metrics)

  _ξ̂xⱼ₊½ = ξ̂xⱼ₊½(meg6_metrics)
  _η̂xⱼ₊½ = η̂xⱼ₊½(meg6_metrics)
  _ξ̂yⱼ₊½ = ξ̂yⱼ₊½(meg6_metrics)
  _η̂yⱼ₊½ = η̂yⱼ₊½(meg6_metrics)

  for idx in i₊½CI
    i₊½_metric = (
      ξ̂x=_ξ̂xᵢ₊½[idx],
      ξ̂y=_ξ̂yᵢ₊½[idx],
      η̂x=_η̂xᵢ₊½[idx],
      η̂y=_η̂yᵢ₊½[idx],
      ξt=zero(Float64),
      ηt=zero(Float64),
    )

    m.edge_metrics.i₊½[idx] = i₊½_metric
  end

  for idx in j₊½CI
    j₊½_metric = (
      ξ̂x=_ξ̂xⱼ₊½[idx],
      ξ̂y=_ξ̂yⱼ₊½[idx],
      η̂x=_η̂xⱼ₊½[idx],
      η̂y=_η̂yⱼ₊½[idx],
      ξt=zero(Float64),
      ηt=zero(Float64),
    )

    m.edge_metrics.j₊½[idx] = j₊½_metric
  end

  return nothing
end

function update_metrics_old!(m::CurvilinearGrid2D)
  function _get_metrics(m, (i, j))
    _jacobian_matrix = m.jacobian_matrix_func(i, j)
    inv_jacobian_matrix = inv(_jacobian_matrix)
    J = det(_jacobian_matrix)
    ξx = inv_jacobian_matrix[1, 1]
    ξy = inv_jacobian_matrix[1, 2]
    ηx = inv_jacobian_matrix[2, 1]
    ηy = inv_jacobian_matrix[2, 2]
    return (J, ξx, ξy, ηx, ηy)
  end

  # cell centroid metrics
  for idx in m.iterators.cell.domain
    node_idx = idx.I .- m.nhalo
    cell_idx = node_idx .+ 0.5
    J, ξx, ξy, ηx, ηy = _get_metrics(m, cell_idx)

    _metric = (
      J=J,
      ξx=ξx,
      # ξ̂x=ξx / J,
      ξy=ξy,
      # ξ̂y=ξy / J,
      ηx=ηx,
      # η̂x=ηx / J,
      ηy=ηy,
      # η̂y=ηy / J,
      ξt=zero(J),
      ηt=zero(J),
    )

    m.cell_center_metrics[idx] = _metric
  end

  ilo_n, ihi_n, jlo_n, jhi_n = m.domain_limits.node

  j₊½CI = CartesianIndices((ilo_n:(ihi_n - 1), jlo_n:jhi_n))
  i₊½CI = CartesianIndices((ilo_n:ihi_n, jlo_n:(jhi_n - 1)))

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

  for idx in i₊½CI
    node_idx = idx.I
    i₊½_node_idx = (node_idx[1], node_idx[2] + 0.5) .- m.nhalo
    i₊½_cell_idx = CartesianIndex((node_idx[1] - 1, node_idx[2]))

    J_i₊½, ξx_i₊½, ξy_i₊½, ηx_i₊½, ηy_i₊½ = _get_metrics(m, i₊½_node_idx)

    i₊½_metric = (
      # J=J_i₊½,
      # ξx=ξx_i₊½,
      ξ̂x=ξx_i₊½ * J_i₊½,
      # ξy=ξy_i₊½,
      ξ̂y=ξy_i₊½ * J_i₊½,
      # ηx=ηx_i₊½,
      η̂x=ηx_i₊½ * J_i₊½,
      # ηy=ηy_i₊½,
      η̂y=ηy_i₊½ * J_i₊½,
      ξt=zero(J_i₊½),
      ηt=zero(J_i₊½),
    )

    m.edge_metrics.i₊½[i₊½_cell_idx] = i₊½_metric
  end

  for idx in j₊½CI
    node_idx = idx.I
    j₊½_node_idx = (node_idx[1] + 0.5, node_idx[2]) .- m.nhalo
    j₊½_cell_idx = CartesianIndex((node_idx[1], node_idx[2] - 1))
    J_j₊½, ξx_j₊½, ξy_j₊½, ηx_j₊½, ηy_j₊½ = _get_metrics(m, j₊½_node_idx)

    j₊½_metric = (
      # J=J_j₊½,
      # ξx=ξx_j₊½,
      ξ̂x=ξx_j₊½ * J_j₊½,
      # ξy=ξy_j₊½,
      ξ̂y=ξy_j₊½ * J_j₊½,
      # ηx=ηx_j₊½,
      η̂x=ηx_j₊½ * J_j₊½,
      # ηy=ηy_j₊½,
      η̂y=ηy_j₊½ * J_j₊½,
      ξt=zero(J_j₊½),
      ηt=zero(J_j₊½),
    )

    m.edge_metrics.j₊½[j₊½_cell_idx] = j₊½_metric
  end

  return nothing
end

# Get the conservative metrics, e.g. normalized by the Jacobian
@inline function conservative_metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  _jacobian_matrix = checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo))
  inv_jacobian_matrix = inv(_jacobian_matrix)
  J = det(_jacobian_matrix)
  J⁻¹ = det(inv_jacobian_matrix)

  # @show J, J⁻¹
  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]

  return (
    ξ̂x=ξx / J⁻¹,
    ξ̂y=ξy / J⁻¹,
    η̂x=ηx / J⁻¹,
    η̂y=ηy / J⁻¹,
    ξt=zero(eltype(_jacobian_matrix)),
    ηt=zero(eltype(_jacobian_matrix)),
  )
end

@inline function conservative_metrics(
  m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, (vx, vy)
)
  static = conservative_metrics(m, (i, j))
  @unpack ξ̂x, ξ̂y, η̂x, η̂y = static

  return merge(static, (
    ξt=-(vx * ξ̂x + vy * ξ̂y), # dynamic / moving mesh terms
    ηt=-(vx * η̂x + vy * η̂y), # dynamic / moving mesh terms
  ))
end

@inline function metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  _jacobian_matrix = m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo) # -> SMatrix{2,2,T}
  inv_jacobian_matrix = inv(_jacobian_matrix)

  return (
    ξx=inv_jacobian_matrix[1, 1],
    ξy=inv_jacobian_matrix[1, 2],
    ηx=inv_jacobian_matrix[2, 1],
    ηy=inv_jacobian_matrix[2, 2],
    ξt=zero(eltype(_jacobian_matrix)),
    ηt=zero(eltype(_jacobian_matrix)),
    J=det(_jacobian_matrix),
  )
end

# @inline function metrics_with_jacobian(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
#   _jacobian_matrix = checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo)) # -> SMatrix{2,2,T}
#   inv_jacobian_matrix = inv(_jacobian_matrix)

#   return (
#     ξx=inv_jacobian_matrix[1, 1],
#     ξy=inv_jacobian_matrix[1, 2],
#     ηx=inv_jacobian_matrix[2, 1],
#     ηy=inv_jacobian_matrix[2, 2],
#     ξt=zero(eltype(_jacobian_matrix)),
#     ηt=zero(eltype(_jacobian_matrix)),
#     J=det(_jacobian_matrix),
#   )
# end

@inline function metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, (vx, vy))
  static = metrics(m, (i, j))
  @unpack ξx, ξy, ηx, ηy = static

  return merge(static, (
    ξt=-(vx * ξx + vy * ξy), # dynamic / moving mesh terms
    ηt=-(vx * ηx + vy * ηy), # dynamic / moving mesh terms
  ))
end

function jacobian_matrix(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  return checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo))
end

function jacobian(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  return det(jacobian_matrix(m, (i, j)))
end

"""
    coords(mesh::CurvilinearGrid2D, T=Float64) -> Array{Real}

Return the array of coordinate points, indexed as `[xy,i,j]`.
This does _not_ include halo regions since the geometry can be undefined.
"""
function coords!(xy::Array{T,3}, m::CurvilinearGrid2D) where {T}
  @inbounds for j in axes(xy, 3)
    for i in axes(xy, 2)
      xy[1, i, j] = m.x_func(i, j)
      xy[2, i, j] = m.y_func(i, j)
    end
  end

  return nothing
end

function coords(m::CurvilinearGrid2D, T=Float64)
  xy = zeros(T, 2, m.nnodes...)
  coords!(xy, m)
  return xy
end

function coords_withhalo(m::CurvilinearGrid2D, T=Float64)
  xy = zeros(T, 2, (m.nnodes .+ 2nhalo)...)

  @inbounds for j in axes(xy, 3)
    for i in axes(xy, 2)
      xy[1, i, j] = m.x_func(i - m.nhalo, j - m.nhalo)
      xy[2, i, j] = m.y_func(i - m.nhalo, j - m.nhalo)
    end
  end
  return xy
end

"""
    centroids(m::CurvilinearGrid2D, T=Float64) -> Array{Real}

Return the array of coordinate points, indexed as `[xy,i,j]`.
This does _not_ include halo regions since the geometry can be undefined.
"""
function centroids(m::CurvilinearGrid2D, T=Float64)
  xy = zeros(T, 2, (m.nnodes .- 1)...)

  @inbounds for j in axes(xy, 3)
    for i in axes(xy, 2)
      xy[1, i, j] = m.x_func(i + 0.5, j + 0.5)
      xy[2, i, j] = m.y_func(i + 0.5, j + 0.5)
    end
  end

  return xy
end
