
struct CurvilinearMesh3D{T,T1,T2,T3,T4,T5} <: AbstractCurvilinearMesh
  x::T1 # x(ξ,η,ζ)
  y::T2 # y(ξ,η,ζ)
  z::T3 # z(ξ,η,ζ)
  jacobian_matrix_func::T4 # jacobian_matrix(ξ,η,ζ)
  conserv_metric_func::T5 # f(ξ,η,ζ) to get ξ̂x, ξ̂y, ...
  nhalo::Int # number of halo cells (for all dimensions)
  nnodes::NTuple{3,Int}
  limits::NamedTuple{(:ilo, :ihi, :jlo, :jhi, :klo, :khi),NTuple{6,Int}}
  cell_center_metrics::Array{Metrics3D{T},3}
  edge_metrics::NamedTuple{
    (:ξ, :η, :ζ), # one for each edge direction
    NTuple{
      3, # (ξ,η,ζ)
      Array{Metrics3D{T},3}, # an array to hold the 3d metrics
    },
  }

  J::Array{T,3} # cell-centered Jacobian J
  # J⁻¹::Array{T,3} # cell-centered inverse Jacobian J⁻¹
end

function CurvilinearMesh3D(x::Function, y::Function, z::Function, (n_ξ, n_η, n_ζ), nhalo)
  jacobian_matrix_func = _setup_jacobian_func(x, y, z)
  cons_metric_func = _setup_conservative_metrics_func(x, y, z)
  nnodes = (n_ξ, n_η, n_ζ)
  ni_cells = n_ξ - 1
  nj_cells = n_η - 1
  nk_cells = n_ζ - 1
  lo = nhalo + 1
  cell_dims = (ni_cells, nj_cells, nk_cells)

  limits = (
    ilo=lo, ihi=ni_cells - nhalo, jlo=lo, jhi=nj_cells - nhalo, klo=lo, khi=nk_cells - nhalo
  )

  cell_center_metrics = empty_metrics(cell_dims)

  # metric terms for the {i,j,k}±1/2 edges
  edge_metrics = (
    ξ=empty_metrics(cell_dims), η=empty_metrics(cell_dims), ζ=empty_metrics(cell_dims)
  )

  J = zeros(cell_dims)
  # J⁻¹ = zeros(cell_dims)
  m = CurvilinearMesh3D(
    x,
    y,
    z,
    jacobian_matrix_func,
    cons_metric_func,
    nhalo,
    nnodes,
    limits,
    cell_center_metrics,
    edge_metrics,
    J,
    # J⁻¹,
  )

  # initialize the metric terms
  update(m)

  return m
end

function coords(m::CurvilinearMesh3D)
  dims = m.nnodes
  xyz = zeros(3, dims...)
  @inbounds for I in CartesianIndices(dims)
    i, j, k = Tuple(I)
    xyz[1, i, j, k] = m.x(i, j, k)
    xyz[2, i, j, k] = m.y(i, j, k)
    xyz[3, i, j, k] = m.z(i, j, k)
  end

  return xyz
end

function centroids(m::CurvilinearMesh3D)
  dims = (m.nnodes .- 1)
  xyz = zeros(3, dims...)
  @inbounds for I in CartesianIndices(dims)
    i, j, k = Tuple(I)
    xyz[1, i, j, k] = m.x(i + 0.5, j + 0.5, k + 0.5)
    xyz[2, i, j, k] = m.y(i + 0.5, j + 0.5, k + 0.5)
    xyz[3, i, j, k] = m.z(i + 0.5, j + 0.5, k + 0.5)
  end

  return xyz
end

"""
    update_metrics(m::CurvilinearMesh3D)

Update the conservative grid metrics (ξ̂x, ξ̂y, ...). This only needs to be called once for a 
static mesh. In a dynamic mesh, this needs to be called each time the mesh moves.
"""
function update_metrics(m::CurvilinearMesh3D)
  # cell-centered metrics
  @inline for I in CartesianIndices(m.cell_center_metrics)
    ξ, η, ζ = Tuple(I)
    m.cell_center_metrics[I] = _get_metric_terms(m, (ξ, η, ζ))
  end

  # edge metric terms
  @inline for I in CartesianIndices(m.cell_center_metrics)
    ξ, η, ζ = Tuple(I)
    m.edge_metrics.ξ[I] = _get_metric_terms(m, (ξ + 0.5, η, ζ))
    m.edge_metrics.η[I] = _get_metric_terms(m, (ξ, η + 0.5, ζ))
    m.edge_metrics.ζ[I] = _get_metric_terms(m, (ξ, η, ζ + 0.5))
  end

  return nothing
end

@inline function _get_metric_terms(m::CurvilinearMesh3D, (ξ, η, ζ))
  M1, M2, M3 = m.conserv_metric_func(ξ, η, ζ) # get the matrices 

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
    η̂x=yζzξ - yξzζ,
    ζ̂x=yξzη - yηzξ,
    ξ̂y=zηxζ - zζxη,
    η̂y=zζxξ - zξxζ,
    ζ̂y=zξxη - zηxξ,
    ξ̂z=xηyζ - xζyη,
    η̂z=xζyξ - xξyζ,
    ζ̂z=xξyη - xηyξ,
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

# Get the grid metrics for a static grid
function _get_metrics(_jacobian_matrix::SMatrix{3,3,T}) where {T}
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

  return (
    ξ̂x=ξx * J⁻¹,
    ξ̂y=ξy * J⁻¹,
    ξ̂z=ξz * J⁻¹,
    η̂x=ηx * J⁻¹,
    η̂y=ηy * J⁻¹,
    η̂z=ηz * J⁻¹,
    ζ̂x=ζx * J⁻¹,
    ζ̂y=ζy * J⁻¹,
    ζ̂z=ζz * J⁻¹,
  )
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
  @inline gradx(ξ, η, ζ) = ForwardDiff.gradient(x, @SVector [ξ, η, ζ]) # ∂x∂ξ, ∂x∂η, ∂x∂ζ
  @inline grady(ξ, η, ζ) = ForwardDiff.gradient(y, @SVector [ξ, η, ζ]) # ∂z∂ξ, ∂z∂η, ∂z∂ζ
  @inline gradz(ξ, η, ζ) = ForwardDiff.gradient(z, @SVector [ξ, η, ζ]) # ∂y∂ξ, ∂y∂η, ∂y∂ζ

  @inline gradx_y(ξ, η, ζ) = gradx(ξ, η, ζ) * y(ξ, η, ζ)
  @inline grady_z(ξ, η, ζ) = grady(ξ, η, ζ) * z(ξ, η, ζ)
  @inline gradz_x(ξ, η, ζ) = gradz(ξ, η, ζ) * x(ξ, η, ζ)

  # these are needed to evaluate as "f(x)", with a single argument
  @inline gradx_y((ξ, η, ζ)) = gradx_y(ξ, η, ζ)
  @inline grady_z((ξ, η, ζ)) = grady_z(ξ, η, ζ)
  @inline gradz_x((ξ, η, ζ)) = gradz_x(ξ, η, ζ)

  @inline grady_z_jacobian(ξ, η, ζ) = ForwardDiff.jacobian(grady_z, @SVector [ξ, η, ζ])
  @inline gradz_x_jacobian(ξ, η, ζ) = ForwardDiff.jacobian(gradz_x, @SVector [ξ, η, ζ])
  @inline gradx_y_jacobian(ξ, η, ζ) = ForwardDiff.jacobian(gradx_y, @SVector [ξ, η, ζ])

  # Get all the matrices at once
  @inline function conserv_metric_matricies(ξ, η, ζ)
    return (grady_z_jacobian(ξ, η, ζ), gradz_x_jacobian(ξ, η, ζ), gradx_y_jacobian(ξ, η, ζ))
  end
  conserv_metric_matricies((ξ, η, ζ)) = conserv_metric_matricies(ξ, η, ζ)

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
