
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
struct CurvilinearGrid3D{F1,F2,F3,F4,F5} <: AbstractCurvilinearGrid
  x::F1 # x(ξ,η,ζ)
  y::F2 # y(ξ,η,ζ)
  z::F3 # z(ξ,η,ζ)
  jacobian_matrix_func::F4 # jacobian_matrix(ξ,η,ζ)
  conserv_metric_func::F5 # f(ξ,η,ζ) to get ξ̂x, ξ̂y, ...
  nhalo::Int # number of halo cells (for all dimensions)
  nnodes::NTuple{3,Int}
  limits::NamedTuple{(:ilo, :ihi, :jlo, :jhi, :klo, :khi),NTuple{6,Int}}
end

function CurvilinearGrid3D(x::Function, y::Function, z::Function, (n_ξ, n_η, n_ζ), nhalo)
  jacobian_matrix_func = _setup_jacobian_func(x, y, z)
  cons_metric_func = _setup_conservative_metrics_func(x, y, z)
  nnodes = (n_ξ, n_η, n_ζ)
  ni_cells = n_ξ - 1
  nj_cells = n_η - 1
  nk_cells = n_ζ - 1
  lo = nhalo + 1

  limits = (
    ilo=lo, ihi=ni_cells - nhalo, jlo=lo, jhi=nj_cells - nhalo, klo=lo, khi=nk_cells - nhalo
  )

  return CurvilinearGrid3D(
    x, y, z, jacobian_matrix_func, cons_metric_func, nhalo, nnodes, limits
  )
end

@inline function conservative_metrics(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real})
  M1, M2, M3 = m.conserv_metric_func(i - m.nhalo, j - m.nhalo, k - m.nhalo) # get the matrices 

  yξzη = M1[1, 2] # ∂(z ∂y/∂ξ)/∂η
  yξzζ = M1[1, 3] # ∂(z ∂y/∂ξ)/∂ζ
  yηzξ = M1[2, 1] # ∂(z ∂y/∂η)/∂ξ
  yηzζ = M1[2, 3] # ∂(z ∂y/∂η)/∂ζ
  yζzξ = M1[3, 1] # ∂(z ∂y/∂ζ)/∂ξ
  yζzη = M1[3, 2] # ∂(z ∂y/∂ζ)/∂η

  zξxη = M2[1, 2] # ∂(x ∂z/∂ξ)/∂η
  zξxζ = M2[1, 3] # ∂(x ∂z/∂ξ)/∂ζ
  zηxξ = M2[2, 1] # ∂(x ∂z/∂η)/∂ξ
  zηxζ = M2[2, 3] # ∂(x ∂z/∂η)/∂ζ
  zζxξ = M2[3, 1] # ∂(x ∂z/∂ζ)/∂ξ
  zζxη = M2[3, 2] # ∂(x ∂z/∂ζ)/∂η

  xξyη = M3[1, 2] # ∂(y ∂x/∂ξ)/∂η
  xξyζ = M3[1, 3] # ∂(y ∂x/∂ξ)/∂ζ
  xηyξ = M3[2, 1] # ∂(y ∂x/∂η)/∂ξ
  xηyζ = M3[2, 3] # ∂(y ∂x/∂η)/∂ζ
  xζyξ = M3[3, 1] # ∂(y ∂x/∂ζ)/∂ξ
  xζyη = M3[3, 2] # ∂(y ∂x/∂ζ)/∂η

  return (
    ξ̂x=yηzζ - yζzη,
    ξ̂y=zηxζ - zζxη,
    ξ̂z=xηyζ - xζyη,
    η̂x=yζzξ - yξzζ,
    η̂y=zζxξ - zξxζ,
    η̂z=xζyξ - xξyζ,
    ζ̂x=yξzη - yηzξ,
    ζ̂y=zξxη - zηxξ,
    ζ̂z=xξyη - xηyξ,
    ξt=zero(eltype(M1)),
    ηt=zero(eltype(M1)),
    ζt=zero(eltype(M1)),
  )
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
  )
end

@inline function metrics(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, (vx, vy, vz))
  static = metrics(m, (i, j))
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

function _setup_jacobian_func(x, y, z)

  # for the gradient calls to work, the x,y,z functions
  # need to work with a single argument, so we make these
  # functions below:
  error_found = false
  for func in (x, y, z)
    try
      func((1, 1, 1))
    catch
      @warn(
        "The $(Symbol(func)) function needs to have definitions for both $(Symbol(func))(i,j,k) and $(Symbol(func))((i,j,k)), or else the gradient calculations will fail",
      )
      error_found = true
    end
  end

  if error_found
    error("Grid function definitions aren't complete... see the warning messages")
  end

  ∇x(ξ, η, ζ) = ForwardDiff.gradient(x, @SVector [ξ, η, ζ]) # ∂x∂ξ, ∂x∂η, ∂x∂ζ
  ∇y(ξ, η, ζ) = ForwardDiff.gradient(y, @SVector [ξ, η, ζ]) # ∂z∂ξ, ∂z∂η, ∂z∂ζ
  ∇z(ξ, η, ζ) = ForwardDiff.gradient(z, @SVector [ξ, η, ζ]) # ∂y∂ξ, ∂y∂η, ∂y∂ζ

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
  gradx_y((ξ, η, ζ)) = ForwardDiff.gradient(x, @SVector [ξ, η, ζ]) * y(ξ, η, ζ)
  grady_z((ξ, η, ζ)) = ForwardDiff.gradient(y, @SVector [ξ, η, ζ]) * z(ξ, η, ζ)
  gradz_x((ξ, η, ζ)) = ForwardDiff.gradient(z, @SVector [ξ, η, ζ]) * x(ξ, η, ζ)

  grady_z_jacobian(ξ, η, ζ) = ForwardDiff.jacobian(grady_z, @SVector [ξ, η, ζ])
  gradz_x_jacobian(ξ, η, ζ) = ForwardDiff.jacobian(gradz_x, @SVector [ξ, η, ζ])
  gradx_y_jacobian(ξ, η, ζ) = ForwardDiff.jacobian(gradx_y, @SVector [ξ, η, ζ])

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