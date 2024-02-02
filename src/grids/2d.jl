
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
    cell_center_metrics = (
      J=zeros(T, _iterators.cell.full.indices),
      ξ=[Metric2D(T) for _ in _iterators.cell.full],
      η=[Metric2D(T) for _ in _iterators.cell.full],
    )

    edge_metrics = (
      i₊½=(
        ξ̂=[Metric2D(T) for _ in _iterators.cell.full],
        η̂=[Metric2D(T) for _ in _iterators.cell.full],
      ),
      j₊½=(
        ξ̂=[Metric2D(T) for _ in _iterators.cell.full],
        η̂=[Metric2D(T) for _ in _iterators.cell.full],
      ),
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
  # update_metrics_old!(m)
  update_metrics!(m)
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

function update_metrics!(m::CurvilinearGrid2D, t=0)

  # cell metrics
  for idx in m.iterators.cell.full
    cell_idx = idx.I .+ 0.5
    @unpack J, ξ, η = metrics(m, cell_idx, t)

    m.cell_center_metrics.ξ[idx] = ξ
    m.cell_center_metrics.η[idx] = η
    m.cell_center_metrics.J[idx] = J
  end

  # i₊½ conserved metrics
  for idx in m.iterators.cell.full
    i, j = idx.I .+ 0.5 # centroid index

    # get the conserved metrics at (i₊½, j)
    @unpack ξ̂, η̂ = conservative_metrics(m, (i + 1 / 2, j), t)

    m.edge_metrics.i₊½.ξ̂[idx] = ξ̂
    m.edge_metrics.i₊½.η̂[idx] = η̂
  end

  # j₊½ conserved metrics
  for idx in m.iterators.cell.full
    i, j = idx.I .+ 0.5 # centroid index

    # get the conserved metrics at (i, j₊½)
    @unpack ξ̂, η̂ = conservative_metrics(m, (i, j + 1 / 2), t)

    m.edge_metrics.j₊½.ξ̂[idx] = ξ̂
    m.edge_metrics.j₊½.η̂[idx] = η̂
  end

  return nothing
end

# Get the conservative metrics, e.g.  ξ̂x = ξx * J

@inline conservative_metrics(m::CurvilinearGrid2D, idx) = conservative_metrics(m, idx, 0)

@inline function conservative_metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t::Real)
  _jacobian_matrix = checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo))
  J = det(_jacobian_matrix)
  inv_jacobian_matrix = inv(_jacobian_matrix)
  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]

  ξt = zero(J)
  ηt = zero(J)
  ξ̂ = Metric2D(ξx * J, ξy * J, ξt * J)
  η̂ = Metric2D(ηx * J, ηy * J, ηt * J)

  return (; ξ̂, η̂)
end

@inline function conservative_metrics(
  m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t::Real, v⃗::SVector{2,Real}
)
  _jacobian_matrix = checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo))
  J = det(_jacobian_matrix)
  inv_jacobian_matrix = inv(_jacobian_matrix)
  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]

  # dynamic / moving mesh terms
  ξt = -(v⃗.x * ξ̂x + v⃗.y * ξ̂y)
  ηt = -(v⃗.x * η̂x + v⃗.y * η̂y)

  ξ̂ = Metric2D(ξx * J, ξy * J, ξt * J)
  η̂ = Metric2D(ηx * J, ηy * J, ηt * J)

  return (; ξ̂, η̂)
end

@inline function metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t::Real)
  _jacobian_matrix = m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo)
  inv_jacobian_matrix = inv(_jacobian_matrix)

  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]

  J = det(_jacobian_matrix)

  ξt = zero(eltype(_jacobian_matrix))
  ηt = zero(eltype(_jacobian_matrix))

  ξ = Metric2D(ξx, ξy, ξt)
  η = Metric2D(ηx, ηy, ηt)

  return (; ξ, η, J)
end

# @inline function metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, (vx, vy))
#   static = metrics(m, (i, j))
#   @unpack ξx, ξy, ηx, ηy = static

#   return merge(static, (
#     ξt=-(vx * ξx + vy * ξy), # dynamic / moving mesh terms
#     ηt=-(vx * ηx + vy * ηy), # dynamic / moving mesh terms
#   ))
# end

function grid_velocity(m::CurvilinearGrid2D, (i, j), t)
  x = m.x_func(i, j, t)
  xₜ = 0
  yₜ = 0

  # xy(i, j) = @SVector [x(i, j), y(i, j)]
  # ForwardDiff.derivative(f, @SVector [i, j, t])

  return xₜ, yₜ
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
