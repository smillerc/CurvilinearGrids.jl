
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
struct CurvilinearGrid3D{F1,F2,F3,EM,CCM,F4,F5,L,CI} <: AbstractCurvilinearGrid
  x::F1 # x(ijk)
  y::F2 # y(ijk)
  z::F3 # z(ijk)
  edge_metrics::EM
  cell_center_metrics::CCM
  jacobian_matrix_func::F4 # jacobian_matrix(ξ,η,ζ)
  conserv_metric_func::F5 # f(ξ,η,ζ) to get ξ̂x, ξ̂y, ...
  nhalo::Int # number of halo cells (for all dimensions)
  nnodes::NTuple{3,Int}
  domain_limits::L
  iterators::CI
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
  cons_metric_func = _setup_conservative_metrics_func(x, y, z)
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

  if cache
    _metric = (
      J=zero(T),
      ξx=zero(T),
      ξ̂x=zero(T),
      ξy=zero(T),
      ξ̂y=zero(T),
      ξz=zero(T),
      ξ̂z=zero(T),
      ηx=zero(T),
      η̂x=zero(T),
      ηy=zero(T),
      η̂y=zero(T),
      ηz=zero(T),
      η̂z=zero(T),
      ζx=zero(T),
      ζ̂x=zero(T),
      ζy=zero(T),
      ζ̂y=zero(T),
      ζz=zero(T),
      ζ̂z=zero(T),
      ξt=zero(T),
      ηt=zero(T),
      ζt=zero(T),
    )

    M = typeof(_metric)
    cell_center_metrics = Array{M,3}(undef, ncells .+ 2nhalo)

    edge_metrics = (
      i₊½=Array{M,3}(undef, ncells .+ 2nhalo),
      j₊½=Array{M,3}(undef, ncells .+ 2nhalo),
      k₊½=Array{M,3}(undef, ncells .+ 2nhalo),
    )
  else
    edge_metrics = nothing
    cell_center_metrics = nothing
  end

  m = CurvilinearGrid3D(
    x,
    y,
    z,
    edge_metrics,
    cell_center_metrics,
    jacobian_matrix_func,
    cons_metric_func,
    nhalo,
    nnodes,
    limits,
    _iterators,
  )

  update_metrics!(m)

  return m
end

function update_metrics!(m::CurvilinearGrid3D)
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

@inline function conservative_metrics(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real})
  _jacobian_matrix = checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo))
  inv_jacobian_matrix = inv(_jacobian_matrix)
  J = jacobian(m, (i, j, k))
  return (
    ξ̂x=inv_jacobian_matrix[1, 1] / J,
    ξ̂y=inv_jacobian_matrix[1, 2] / J,
    ξ̂z=inv_jacobian_matrix[1, 3] / J,
    η̂x=inv_jacobian_matrix[2, 1] / J,
    η̂y=inv_jacobian_matrix[2, 2] / J,
    η̂z=inv_jacobian_matrix[2, 3] / J,
    ζ̂x=inv_jacobian_matrix[3, 1] / J,
    ζ̂y=inv_jacobian_matrix[3, 2] / J,
    ζ̂z=inv_jacobian_matrix[3, 3] / J,
    ξt=zero(eltype(_jacobian_matrix)),
    ηt=zero(eltype(_jacobian_matrix)),
    ζt=zero(eltype(_jacobian_matrix)),
  )

  # M1, M2, M3 = m.conserv_metric_func(i - m.nhalo, j - m.nhalo, k - m.nhalo) # get the matrices

  # yξzη = M1[1, 2] # ∂(z ∂y/∂ξ)/∂η
  # yξzζ = M1[1, 3] # ∂(z ∂y/∂ξ)/∂ζ
  # yηzξ = M1[2, 1] # ∂(z ∂y/∂η)/∂ξ
  # yηzζ = M1[2, 3] # ∂(z ∂y/∂η)/∂ζ
  # yζzξ = M1[3, 1] # ∂(z ∂y/∂ζ)/∂ξ
  # yζzη = M1[3, 2] # ∂(z ∂y/∂ζ)/∂η

  # zξxη = M2[1, 2] # ∂(x ∂z/∂ξ)/∂η
  # zξxζ = M2[1, 3] # ∂(x ∂z/∂ξ)/∂ζ
  # zηxξ = M2[2, 1] # ∂(x ∂z/∂η)/∂ξ
  # zηxζ = M2[2, 3] # ∂(x ∂z/∂η)/∂ζ
  # zζxξ = M2[3, 1] # ∂(x ∂z/∂ζ)/∂ξ
  # zζxη = M2[3, 2] # ∂(x ∂z/∂ζ)/∂η

  # xξyη = M3[1, 2] # ∂(y ∂x/∂ξ)/∂η
  # xξyζ = M3[1, 3] # ∂(y ∂x/∂ξ)/∂ζ
  # xηyξ = M3[2, 1] # ∂(y ∂x/∂η)/∂ξ
  # xηyζ = M3[2, 3] # ∂(y ∂x/∂η)/∂ζ
  # xζyξ = M3[3, 1] # ∂(y ∂x/∂ζ)/∂ξ
  # xζyη = M3[3, 2] # ∂(y ∂x/∂ζ)/∂η

  # return (
  #   ξ̂x=yηzζ - yζzη,
  #   ξ̂y=zηxζ - zζxη,
  #   ξ̂z=xηyζ - xζyη,
  #   η̂x=yζzξ - yξzζ,
  #   η̂y=zζxξ - zξxζ,
  #   η̂z=xζyξ - xξyζ,
  #   ζ̂x=yξzη - yηzξ,
  #   ζ̂y=zξxη - zηxξ,
  #   ζ̂z=xξyη - xηyξ,
  #   ξt=zero(eltype(M1)),
  #   ηt=zero(eltype(M1)),
  #   ζt=zero(eltype(M1)),
  # )
end

@inline function conservative_metrics(
  m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, (vx, vy, vz)
)
  static = conservative_metrics(m, (i, j, k))
  @unpack ξ̂x, ξ̂y, ξ̂z, η̂x, η̂y, η̂z, ζ̂x, ζ̂y, ζ̂z = static

  return merge(
    static,
    (
      ξt=-(vx * ξ̂x + vy * ξ̂y + vz * ξ̂z), # dynamic / moving mesh terms
      ηt=-(vx * η̂x + vy * η̂y + vz * η̂z), # dynamic / moving mesh terms
      ζt=-(vx * ζ̂x + vy * ζ̂y + vz * ζ̂z), # dynamic / moving mesh terms
    ),
  )
end

# Get the grid metrics for a static grid
@inline function metrics(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real})
  _jacobian_matrix = checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo))
  inv_jacobian_matrix = inv(_jacobian_matrix)

  return (
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
    J=det(_jacobian_matrix),
  )
end

@inline function metrics(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, (vx, vy, vz))
  static = metrics(m, (i, j, k))
  @unpack ξx, ξy, ξz, ηx, ηy, ηz, ζx, ζy, ζz = static

  return merge(
    static,
    (
      ξt=-(vx * ξx + vy * ξy + vz * ξz), # dynamic / moving mesh terms
      ηt=-(vx * ηx + vy * ηy + vz * ηz), # dynamic / moving mesh terms
      ζt=-(vx * ζx + vy * ζy + vz * ζz), # dynamic / moving mesh terms
    ),
  )
end

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
    xyz[1, i, j, k] = m.x(i, j, k)
    xyz[2, i, j, k] = m.y(i, j, k)
    xyz[3, i, j, k] = m.z(i, j, k)
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
