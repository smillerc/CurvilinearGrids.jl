
struct RZAxisymmetricGrid2D{T1,T2,T3,T4} <: AbstractCurvilinearGrid
  r::T1
  θ::T2
  z::T3
  jacobian_matrix_func::T4
  nhalo::Int
  nnodes::NTuple{2,Int}
  limits::NamedTuple{(:ilo, :ihi, :jlo, :jhi),NTuple{4,Int}}
end

function RZAxisymmetricGrid2D(r::Function, z::Function, (ni, nj), nhalo)

  # Ensure that the r and z functions are set up properly, i.e., 
  # they are defined as r(i,j) = ... and z(i,j) = ...
  dim = 2
  check_nargs(r, dim, :r)
  check_nargs(z, dim, :z)
  test_coord_func(r, dim, :r)
  test_coord_func(z, dim, :z)

  # Make a full 3d grid with only 1 cell in θ.
  # This is cheap and very useful for certain applications

  θ1 = 2#π # leave the π off so we can use the more accurate cospi/sinpi functions
  θ(j) = θ1 * (j - 1)

  R3d(i, j, k) = r(i, k) * cospi(θ(j))
  Θ3d(i, j, k) = r(i, k) * sinpi(θ(j))
  Z3d(i, j, k) = z(i, k)

  RΘZ(i, j, k) = @SVector [R3d(i, j, k), Θ3d(i, j, k), Z3d(i, j, k)]
  function jacobian_matrix_func(i, j, k)
    return ForwardDiff.jacobian(x -> RΘZ(x[1], x[2], x[3]), @SVector [i, j, k])
  end

  # jacobian_matrix_func = _setup_jacobian_func(x, y)
  nnodes = (ni, nj)
  ni_cells = ni - 1
  nj_cells = nj - 1
  lo = nhalo + 1
  limits = (ilo=lo, ihi=ni_cells + nhalo, jlo=lo, jhi=nj_cells + nhalo)

  return RZAxisymmetricGrid2D(R3d, Θ3d, Z3d, jacobian_matrix_func, nhalo, nnodes, limits)
end

@inline function metrics(m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real})

  # Get the full 3d jacobian matrix. The 2nd coordinate doesn't matter
  # since it's symmetric about θ
  _jacobian_matrix = jacobian_matrix(m, (i, 1, j))
  T = eltype(_jacobian_matrix)

  inv_jacobian_matrix = inv(_jacobian_matrix)

  # Only extract the ∂(r,z) terms
  ξr = inv_jacobian_matrix[1, 1]
  ξz = inv_jacobian_matrix[1, 3]
  ηr = inv_jacobian_matrix[3, 1]
  ηz = inv_jacobian_matrix[3, 3]

  # In this scenario, J is the volume of the node/cell at (i,j),
  # and it includes the revolution term. This is important!
  _metrics = (
    ξx¹=ξr, # re-name these so the 2D API is consistent
    ξx²=ξz, # 
    ηx¹=ηr, # 
    ηx²=ηz, # 
    ξt=zero(T),
    ηt=zero(T),
    J=det(_jacobian_matrix),
  )

  return _metrics
end

@inline function planar_metrics(m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real})

  # Get the full 3d jacobian matrix. The 2nd coordinate doesn't matter
  # since it's symmetric about θ
  _jacobian_matrix = jacobian_matrix(m, (i, 1, j))
  T = eltype(_jacobian_matrix)

  inv_jacobian_matrix = inv(_jacobian_matrix)

  # Only extract the ∂(r,z) terms
  ξr = inv_jacobian_matrix[1, 1]
  ξz = inv_jacobian_matrix[1, 3]
  ηr = inv_jacobian_matrix[3, 1]
  ηz = inv_jacobian_matrix[3, 3]
  rξ = _jacobian_matrix[1, 1]
  zξ = _jacobian_matrix[1, 3]
  rη = _jacobian_matrix[3, 1]
  zη = _jacobian_matrix[3, 3]
  J = rξ * zη - zξ * rη

  # In this scenario, J is the AREA of the node/cell at (i,j),
  # and it DOES NOT include the revolution term. This is important!
  _metrics = (
    ξx¹=ξr, # re-name these so the 2D API is consistent
    ξx²=ξz, # 
    ηx¹=ηr, # 
    ηx²=ηz, # 
    ξt=zero(T),
    ηt=zero(T),
    J=J,
  )

  return _metrics
end

@inline function metrics(m::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real})
  _jacobian_matrix = jacobian_matrix(m, (i, j, k))
  T = eltype(_jacobian_matrix)
  inv_jacobian_matrix = inv(_jacobian_matrix)

  # In this scenario, J is the true VOLUME of the node/cell at (i,j,k).
  # This is an important distinction!
  # One may ask how a node has a volume, but then you may be accused of
  # being a pain in the ass... 😆 Just think of it as the the volume 
  # occupied by the polyhedral element defined by the centroids 
  # of the surrounding cells.
  _metrics = (
    ξx¹=inv_jacobian_matrix[1, 1],
    ξx²=inv_jacobian_matrix[1, 2],
    ξx³=inv_jacobian_matrix[1, 3],
    ηx¹=inv_jacobian_matrix[2, 1],
    ηx²=inv_jacobian_matrix[2, 2],
    ηx³=inv_jacobian_matrix[2, 3],
    ζx¹=inv_jacobian_matrix[3, 1],
    ζx²=inv_jacobian_matrix[3, 2],
    ζx³=inv_jacobian_matrix[3, 3],
    ξt=zero(T),
    ηt=zero(T),
    ζt=zero(T),
    J=det(_jacobian_matrix),
  )

  return _metrics
end

function jacobian_matrix(m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real})
  _jacobian_matrix = jacobian_matrix(m, (i, 1, j))
  T = eltype(_jacobian_matrix)

  # extract only the 2d portion
  _jacobian_matrix_2d = SMatrix{2,2}(
    _jacobian_matrix[1, 1],
    _jacobian_matrix[3, 1],
    _jacobian_matrix[1, 3],
    _jacobian_matrix[3, 3],
  )

  return checkeps(_jacobian_matrix_2d)
end

function jacobian_matrix(m::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real})
  return checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo, k - m.nhalo))
end

function jacobian(m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real})
  _jacobian_matrix = checkeps(jacobian_matrix(m, (i, j)))
  return det(_jacobian_matrix)
end

function jacobian(m::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real})
  _jacobian_matrix = checkeps(jacobian_matrix(m, (i, j, k)))
  return det(_jacobian_matrix)
end

area(m::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real}) = jacobian(m, (i, k))
function area(::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real})
  return error(
    """
    You're trying to get the area of a 3d index in a RZAxisymmetricGrid2D, 
    which doesn't make physical sense! Use one of the following:
    1. `area(grid, (i,j))`, which will get you the non-rotate area of the node/cell
    2. `volume(grid, (i,j,k))` to get the true "rotate" volume of the node/cell
    """
  )
end

volume(m::RZAxisymmetricGrid2D, (i, j)::NTuple{3,Real}) = jacobian(m, (i, j))
volume(m::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real}) = jacobian(m, (i, k))
