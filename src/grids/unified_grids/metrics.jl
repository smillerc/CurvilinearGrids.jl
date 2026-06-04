# Unified grid metric payloads and metric computation.

#
# Unified metric AoS type
#

"""
    Metric{N,T,M}

Unified AoS metric payload used by unified-grid metric caches.

The payload stores the metric tensor and associated Jacobian determinant for
forward, inverse, or conserved inverse metric forms.

# Fields
  - `jacobian_matrix`: Metric tensor as a static matrix.
  - `J`: Jacobian determinant associated with `jacobian_matrix`.
"""
struct Metric{N,T,M<:StaticMatrix{N,N,T}}
  jacobian_matrix::M
  J::T
end

"""
    ConservedMetric{N,T,M}

Tensor-only metric payload used for conserved face metrics.

# Fields
  - `jacobian_matrix`: Conserved metric tensor as a static matrix.
"""
struct ConservedMetric{N,T,M<:StaticMatrix{N,N,T}}
  jacobian_matrix::M
end

"""
    Metric(jacobian_matrix::StaticMatrix{N,N,T}) where {N,T}

Construct a metric payload from a Jacobian matrix.

# Arguments
  - `jacobian_matrix`: Jacobian matrix for a cell or face metric.

# Returns
`Metric` with `J = det(jacobian_matrix)`.
"""
@inline function Metric(jacobian_matrix::StaticMatrix{N,N,T}) where {N,T}
  Metric{N,T,typeof(jacobian_matrix)}(jacobian_matrix, det(jacobian_matrix))
end

@inline Base.getindex(metric::Metric, i::Int, j::Int) = metric.jacobian_matrix[i, j]
@inline Base.getindex(metric::ConservedMetric, i::Int, j::Int) = metric.jacobian_matrix[
  i, j
]
@inline function Base.inv(metric::Metric{N,T,M}) where {N,T,M<:StaticMatrix{N,N,T}}
  jinv = inv(metric.jacobian_matrix)
  Metric(jinv, det(jinv))
end

@inline function _metric_eltype(::Val{N}, ::Type{T}) where {N,T}
  Metric{N,T,SMatrix{N,N,T,N * N}}
end

@inline function _conserved_metric_eltype(::Val{N}, ::Type{T}) where {N,T}
  ConservedMetric{N,T,SMatrix{N,N,T,N * N}}
end

function Base.zero(::Type{Metric{N,T,M}}) where {N,T,M<:StaticMatrix{N,N,T}}
  Metric(zero(M), zero(T))
end
Base.zero(metric::Metric{N,T,M}) where {N,T,M<:StaticMatrix{N,N,T}} = zero(Metric{N,T,M})

function Base.zero(::Type{ConservedMetric{N,T,M}}) where {N,T,M<:StaticMatrix{N,N,T}}
  ConservedMetric(zero(M))
end
function Base.zero(metric::ConservedMetric{N,T,M}) where {N,T,M<:StaticMatrix{N,N,T}}
  zero(ConservedMetric{N,T,M})
end

@inline function _as_smatrix(::Val{1}, x)
  if x isa Number
    return @SMatrix [x]
  end
  return SMatrix{1,1}(Tuple(x))
end

@inline function _as_smatrix(::Val{2}, x)
  if x isa StaticMatrix
    return SMatrix{2,2}(Tuple(x))
  end
  return SMatrix{2,2}(Tuple(x))
end

@inline function _as_smatrix(::Val{3}, x)
  if x isa StaticMatrix
    return SMatrix{3,3}(Tuple(x))
  end
  return SMatrix{3,3}(Tuple(x))
end


# Continuous metric cache construction.

struct MetricCache{FM,IM,EM}
  forward::FM
  inverse::IM
  edge::EM
end

"""
Trait for selecting how conserved metric terms are reconstructed at faces.
"""
abstract type EdgeInterpolationSchemeTrait end

"""
Face reconstruction from the average of endpoint values.
"""
struct EndpointAverageReconstruction <: EdgeInterpolationSchemeTrait end

"""
Face reconstruction using endpoint values and first derivatives.
"""
struct GradientCorrectedReconstruction <: EdgeInterpolationSchemeTrait end

"""
Face reconstruction using endpoint values plus first- and second-derivative corrections.
"""
struct CurvatureCorrectedReconstruction <: EdgeInterpolationSchemeTrait end

"""
AD Thomas-Lombard conservative metric path.

This preserves the Thomas-Lombard conservative metric algebra using AD-evaluated
metric potentials and optimized pair-table assembly for the active conserved
face rows.
"""
struct ADThomasLombardMetric <: EdgeInterpolationSchemeTrait end

@inline _coefficient_type(x::Number) = typeof(x)
@inline _coefficient_type(x::StaticArray) = eltype(x)
@inline _edge_coefficient_type(x, Δξ) = promote_type(_coefficient_type(x), typeof(Δξ))

@inline function _edge_reconstruct(ϕᵢ, ϕᵢ₊₁, ::EndpointAverageReconstruction)
  T = _edge_coefficient_type(ϕᵢ, 1)
  half = one(T) / T(2)
  return half * (ϕᵢ + ϕᵢ₊₁)
end

@inline function _edge_reconstruct(
  ϕᵢ, ∂ϕ_∂ξᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ::GradientCorrectedReconstruction, Δξ::Real=1
)
  T = _edge_coefficient_type(ϕᵢ, Δξ)
  half = one(T) / T(2)
  Δ = T(Δξ)
  h = half * Δ
  ϕᴸᵢ₊½ = ϕᵢ + h * ∂ϕ_∂ξᵢ
  ϕᴿᵢ₊½ = ϕᵢ₊₁ - h * ∂ϕ_∂ξᵢ₊₁
  return half * (ϕᴸᵢ₊½ + ϕᴿᵢ₊½)
end

@inline function _edge_reconstruct(
  ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁, ::CurvatureCorrectedReconstruction, Δξ::Real=1
)
  T = _edge_coefficient_type(ϕᵢ, Δξ)
  half = one(T) / T(2)
  twelfth = one(T) / T(12)
  Δ = T(Δξ)
  h = half * Δ
  curvature_weight = Δ * Δ * twelfth
  ϕᴸᵢ₊½ = ϕᵢ + h * ∂ϕ_∂ξᵢ + curvature_weight * ∂²ϕ_∂ξ²ᵢ
  ϕᴿᵢ₊½ = ϕᵢ₊₁ - h * ∂ϕ_∂ξᵢ₊₁ + curvature_weight * ∂²ϕ_∂ξ²ᵢ₊₁
  return half * (ϕᴸᵢ₊½ + ϕᴿᵢ₊½)
end

"""
3D metric cache
"""
function MetricCache(
  x::Function,
  y::Function,
  z::Function,
  backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=CurvatureCorrectedReconstruction(),
)
  xξ(t, i, j, k, p) = derivative(ξ -> x(t, ξ, j, k, p), backend, i)
  xη(t, i, j, k, p) = derivative(η -> x(t, i, η, k, p), backend, j)
  xζ(t, i, j, k, p) = derivative(ζ -> x(t, i, j, ζ, p), backend, k)
  yξ(t, i, j, k, p) = derivative(ξ -> y(t, ξ, j, k, p), backend, i)
  yη(t, i, j, k, p) = derivative(η -> y(t, i, η, k, p), backend, j)
  yζ(t, i, j, k, p) = derivative(ζ -> y(t, i, j, ζ, p), backend, k)
  zξ(t, i, j, k, p) = derivative(ξ -> z(t, ξ, j, k, p), backend, i)
  zη(t, i, j, k, p) = derivative(η -> z(t, i, η, k, p), backend, j)
  zζ(t, i, j, k, p) = derivative(ζ -> z(t, i, j, ζ, p), backend, k)

  function jacobian_matrix(t, i, j, k, p)

    # compute the jacobian matrix w/o any extra logic
    function compute_jacobian_matrix(t, i, j, k, p)
      jac = DifferentiationInterface.jacobian(
        u -> SVector(x(t, u..., p), y(t, u..., p), z(t, u..., p)),
        backend,
        @SVector [i, j, k]
      )

      return jac
    end

    jac_matrix = compute_jacobian_matrix(t, i, j, k, p)
    return jac_matrix
    # J = det(jac_matrix)

  end

  forward_metrics = (;
    jacobian=jacobian_matrix,
    J=jacobian,
    xξ=xξ,
    xη=xη,
    xζ=xζ,
    yξ=yξ,
    yη=yη,
    yζ=yζ,
    zξ=zξ,
    zη=zη,
    zζ=zζ,
  )

  inverse_metrics, edge_metrics = get_inverse_metric_terms(
    x, y, z, backend; edge_interpolation_scheme=edge_interpolation_scheme
  )

  # edge_metrics = get_edge_functions_3d(forward_metrics, inverse_metrics, backend)

  return MetricCache(forward_metrics, inverse_metrics, edge_metrics)
end

"""
2D metric cache
"""
function MetricCache(
  x::Function,
  y::Function,
  backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=CurvatureCorrectedReconstruction(),
)
  xξ(t, i, j, p) = derivative(ξ -> x(t, ξ, j, p), backend, i)
  xη(t, i, j, p) = derivative(η -> x(t, i, η, p), backend, j)
  xτ(t, i, j, p) = derivative(τ -> x(τ, i, j, p), backend, t)

  yξ(t, i, j, p) = derivative(ξ -> y(t, ξ, j, p), backend, i)
  yη(t, i, j, p) = derivative(η -> y(t, i, η, p), backend, j)
  yτ(t, i, j, p) = derivative(τ -> y(τ, i, j, p), backend, t)

  function jacobian_matrix(t, i, j, p)

    # compute the jacobian matrix w/o any extra logic
    # function compute_jacobian_matrix(t, i, j, p)
    jac = DifferentiationInterface.jacobian(
      u -> SVector(x(t, u..., p), y(t, u..., p)), backend, @SVector [i, j]
    )

    return jac
  end

  jacobian(t, i, j, p) = det(jacobian_matrix(t, i, j, p))
  function jinv(t, i, j, p)
    J_inv = inv(jacobian_matrix(t, i, j, p))
    return J_inv
  end

  function normalized_jinv(t, i, j, p)
    jac = jacobian_matrix(t, i, j, p)
    J = det(jac)
    J_inv = inv(jac)

    return J_inv .* J
  end

  forward_metrics = (; jacobian=jacobian_matrix, J=jacobian, xξ=xξ, xη=xη, yξ=yξ, yη=yη)

  inverse_metrics = (;
    # ξ̂x=ξ̂x, ξ̂y=ξ̂y, η̂x=η̂x, η̂y=η̂y,
    Jinv=jinv,
    Jinv_norm=normalized_jinv,
  )

  edge_metrics = get_edge_functions_2d(
    forward_metrics,
    inverse_metrics,
    backend;
    edge_interpolation_scheme=edge_interpolation_scheme,
  )

  return MetricCache(forward_metrics, inverse_metrics, edge_metrics)
end

"""
1D metric cache
"""
function MetricCache(
  x::Function,
  backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=CurvatureCorrectedReconstruction(),
)
  xξ(t, i, p) = derivative(ξ -> x(t, ξ, p), backend, i)

  jacobian_matrix(t, i, p) = @SMatrix [xξ(t, i, p)]

  jacobian(t, i, p) = det(jacobian_matrix(t, i, p))
  jinv(t, i, p) = inv(jacobian_matrix(t, i, p))

  function normalized_jinv(t, i, p)
    jac = jacobian_matrix(t, i, p)
    J = det(jac)
    J_inv = inv(jac)

    return J_inv .* J
  end

  forward_metrics = (; jacobian=jacobian_matrix, J=jacobian, xξ=xξ)

  inverse_metrics = (;
    # ξ̂x=ξ̂x, ξ̂y=ξ̂y, η̂x=η̂x, η̂y=η̂y,
    Jinv=jinv,
    Jinv_norm=normalized_jinv,
  )

  edge_metrics = get_edge_functions_1d(
    forward_metrics,
    inverse_metrics,
    backend;
    edge_interpolation_scheme=edge_interpolation_scheme,
  )

  return MetricCache(forward_metrics, inverse_metrics, edge_metrics)
end

edge_functions_3d(ϕ, backend) = edge_functions_3d(ϕ, backend, CurvatureCorrectedReconstruction())

function edge_functions_3d(ϕ, backend, ::EndpointAverageReconstruction)
  function ϕᵢ₊½(t, i, j, k, p)
    _edge_reconstruct(ϕ(t, i, j, k, p), ϕ(t, i + 1, j, k, p), EndpointAverageReconstruction())
  end
  function ϕⱼ₊½(t, i, j, k, p)
    _edge_reconstruct(ϕ(t, i, j, k, p), ϕ(t, i, j + 1, k, p), EndpointAverageReconstruction())
  end
  function ϕₖ₊½(t, i, j, k, p)
    _edge_reconstruct(ϕ(t, i, j, k, p), ϕ(t, i, j, k + 1, p), EndpointAverageReconstruction())
  end
  return (; ϕᵢ₊½, ϕⱼ₊½, ϕₖ₊½)
end

function edge_functions_3d(ϕ, backend, ::GradientCorrectedReconstruction)
  ξ_derivs(t, i, j, k, p) = value_and_derivative(ξ -> ϕ(t, ξ, j, k, p), backend, i)
  η_derivs(t, i, j, k, p) = value_and_derivative(η -> ϕ(t, i, η, k, p), backend, j)
  ζ_derivs(t, i, j, k, p) = value_and_derivative(ζ -> ϕ(t, i, j, ζ, p), backend, k)

  function ϕᵢ₊½(t, i, j, k, p)
    ϕᵢ, ∂ϕ_∂ξᵢ = ξ_derivs(t, i, j, k, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁ = ξ_derivs(t, i + 1, j, k, p)
    _edge_reconstruct(ϕᵢ, ∂ϕ_∂ξᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, GradientCorrectedReconstruction())
  end
  function ϕⱼ₊½(t, i, j, k, p)
    ϕⱼ, ∂ϕ_∂ηⱼ = η_derivs(t, i, j, k, p)
    ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁ = η_derivs(t, i, j + 1, k, p)
    _edge_reconstruct(ϕⱼ, ∂ϕ_∂ηⱼ, ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, GradientCorrectedReconstruction())
  end
  function ϕₖ₊½(t, i, j, k, p)
    ϕₖ, ∂ϕ_∂ζₖ = ζ_derivs(t, i, j, k, p)
    ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁ = ζ_derivs(t, i, j, k + 1, p)
    _edge_reconstruct(ϕₖ, ∂ϕ_∂ζₖ, ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁, GradientCorrectedReconstruction())
  end
  return (; ϕᵢ₊½, ϕⱼ₊½, ϕₖ₊½)
end

function edge_functions_3d(ϕ, backend, ::CurvatureCorrectedReconstruction)
  function ξ_derivs(t, i, j, k, p)
    value_derivative_and_second_derivative(ξ -> ϕ(t, ξ, j, k, p), backend, i)
  end
  function η_derivs(t, i, j, k, p)
    value_derivative_and_second_derivative(η -> ϕ(t, i, η, k, p), backend, j)
  end
  function ζ_derivs(t, i, j, k, p)
    value_derivative_and_second_derivative(ζ -> ϕ(t, i, j, ζ, p), backend, k)
  end

  function ϕᵢ₊½(t, i, j, k, p)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(t, i, j, k, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(t, i + 1, j, k, p)
    _edge_reconstruct(
      ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁, CurvatureCorrectedReconstruction()
    )
  end
  function ϕⱼ₊½(t, i, j, k, p)
    ϕⱼ, ∂ϕ_∂ηⱼ, ∂²ϕ_∂η²ⱼ = η_derivs(t, i, j, k, p)
    ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, ∂²ϕ_∂η²ⱼ₊₁ = η_derivs(t, i, j + 1, k, p)
    _edge_reconstruct(
      ϕⱼ, ∂ϕ_∂ηⱼ, ∂²ϕ_∂η²ⱼ, ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, ∂²ϕ_∂η²ⱼ₊₁, CurvatureCorrectedReconstruction()
    )
  end
  function ϕₖ₊½(t, i, j, k, p)
    ϕₖ, ∂ϕ_∂ζₖ, ∂²ϕ_∂ζ²ₖ = ζ_derivs(t, i, j, k, p)
    ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁, ∂²ϕ_∂ζ²ₖ₊₁ = ζ_derivs(t, i, j, k + 1, p)
    _edge_reconstruct(
      ϕₖ, ∂ϕ_∂ζₖ, ∂²ϕ_∂ζ²ₖ, ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁, ∂²ϕ_∂ζ²ₖ₊₁, CurvatureCorrectedReconstruction()
    )
  end
  return (; ϕᵢ₊½, ϕⱼ₊½, ϕₖ₊½)
end

edge_functions_2d(ϕ, backend) = edge_functions_2d(ϕ, backend, CurvatureCorrectedReconstruction())

function edge_functions_2d(ϕ, backend, ::EndpointAverageReconstruction)
  function ϕᵢ₊½(t, i, j, p)
    _edge_reconstruct(ϕ(t, i, j, p), ϕ(t, i + 1, j, p), EndpointAverageReconstruction())
  end
  function ϕⱼ₊½(t, i, j, p)
    _edge_reconstruct(ϕ(t, i, j, p), ϕ(t, i, j + 1, p), EndpointAverageReconstruction())
  end
  return (; ϕᵢ₊½, ϕⱼ₊½)
end

function edge_functions_2d(ϕ, backend, ::GradientCorrectedReconstruction)
  ξ_derivs(t, i, j, p) = value_and_derivative(ξ -> ϕ(t, ξ, j, p), backend, i)
  η_derivs(t, i, j, p) = value_and_derivative(η -> ϕ(t, i, η, p), backend, j)

  function ϕᵢ₊½(t, i, j, p)
    ϕᵢ, ∂ϕ_∂ξᵢ = ξ_derivs(t, i, j, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁ = ξ_derivs(t, i + 1, j, p)
    _edge_reconstruct(ϕᵢ, ∂ϕ_∂ξᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, GradientCorrectedReconstruction())
  end
  function ϕⱼ₊½(t, i, j, p)
    ϕⱼ, ∂ϕ_∂ηⱼ = η_derivs(t, i, j, p)
    ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁ = η_derivs(t, i, j + 1, p)
    _edge_reconstruct(ϕⱼ, ∂ϕ_∂ηⱼ, ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, GradientCorrectedReconstruction())
  end
  return (; ϕᵢ₊½, ϕⱼ₊½)
end

function edge_functions_2d(ϕ, backend, ::CurvatureCorrectedReconstruction)
  function ξ_derivs(t, i, j, p)
    value_derivative_and_second_derivative(ξ -> ϕ(t, ξ, j, p), backend, i)
  end
  function η_derivs(t, i, j, p)
    value_derivative_and_second_derivative(η -> ϕ(t, i, η, p), backend, j)
  end

  function ϕᵢ₊½(t, i, j, p)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(t, i, j, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(t, i + 1, j, p)
    _edge_reconstruct(
      ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁, CurvatureCorrectedReconstruction()
    )
  end
  function ϕⱼ₊½(t, i, j, p)
    ϕⱼ, ∂ϕ_∂ηⱼ, ∂²ϕ_∂η²ⱼ = η_derivs(t, i, j, p)
    ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, ∂²ϕ_∂η²ⱼ₊₁ = η_derivs(t, i, j + 1, p)
    _edge_reconstruct(
      ϕⱼ, ∂ϕ_∂ηⱼ, ∂²ϕ_∂η²ⱼ, ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, ∂²ϕ_∂η²ⱼ₊₁, CurvatureCorrectedReconstruction()
    )
  end
  return (; ϕᵢ₊½, ϕⱼ₊½)
end

edge_functions_1d(ϕ, backend) = edge_functions_1d(ϕ, backend, CurvatureCorrectedReconstruction())

function edge_functions_1d(ϕ, backend, ::EndpointAverageReconstruction)
  ϕᵢ₊½(t, i, p) = _edge_reconstruct(ϕ(t, i, p), ϕ(t, i + 1, p), EndpointAverageReconstruction())
  return (; ϕᵢ₊½)
end

function edge_functions_1d(ϕ, backend, ::GradientCorrectedReconstruction)
  ξ_derivs(t, i, p) = value_and_derivative(ξ -> ϕ(t, ξ, p), backend, i)
  function ϕᵢ₊½(t, i, p)
    ϕᵢ, ∂ϕ_∂ξᵢ = ξ_derivs(t, i, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁ = ξ_derivs(t, i + 1, p)
    _edge_reconstruct(ϕᵢ, ∂ϕ_∂ξᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, GradientCorrectedReconstruction())
  end
  return (; ϕᵢ₊½)
end

function edge_functions_1d(ϕ, backend, ::CurvatureCorrectedReconstruction)
  ξ_derivs(t, i, p) = value_derivative_and_second_derivative(ξ -> ϕ(t, ξ, p), backend, i)
  function ϕᵢ₊½(t, i, p)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(t, i, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(t, i + 1, p)
    _edge_reconstruct(
      ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁, CurvatureCorrectedReconstruction()
    )
  end
  return (; ϕᵢ₊½)
end

function get_edge_functions_3d(
  forward_metrics,
  inverse_metrics,
  diff_backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=CurvatureCorrectedReconstruction(),
)
  ξ̂xᵢ₊½, ξ̂xⱼ₊½, ξ̂xₖ₊½ = edge_functions_3d(
    inverse_metrics.ξ̂x, diff_backend, edge_interpolation_scheme
  )
  η̂xᵢ₊½, η̂xⱼ₊½, η̂xₖ₊½ = edge_functions_3d(
    inverse_metrics.η̂x, diff_backend, edge_interpolation_scheme
  )
  ζ̂xᵢ₊½, ζ̂xⱼ₊½, ζ̂xₖ₊½ = edge_functions_3d(
    inverse_metrics.ζ̂x, diff_backend, edge_interpolation_scheme
  )
  ξ̂yᵢ₊½, ξ̂yⱼ₊½, ξ̂yₖ₊½ = edge_functions_3d(
    inverse_metrics.ξ̂y, diff_backend, edge_interpolation_scheme
  )
  η̂yᵢ₊½, η̂yⱼ₊½, η̂yₖ₊½ = edge_functions_3d(
    inverse_metrics.η̂y, diff_backend, edge_interpolation_scheme
  )
  ζ̂yᵢ₊½, ζ̂yⱼ₊½, ζ̂yₖ₊½ = edge_functions_3d(
    inverse_metrics.ζ̂y, diff_backend, edge_interpolation_scheme
  )
  ξ̂zᵢ₊½, ξ̂zⱼ₊½, ξ̂zₖ₊½ = edge_functions_3d(
    inverse_metrics.ξ̂z, diff_backend, edge_interpolation_scheme
  )
  η̂zᵢ₊½, η̂zⱼ₊½, η̂zₖ₊½ = edge_functions_3d(
    inverse_metrics.η̂z, diff_backend, edge_interpolation_scheme
  )
  ζ̂zᵢ₊½, ζ̂zⱼ₊½, ζ̂zₖ₊½ = edge_functions_3d(
    inverse_metrics.ζ̂z, diff_backend, edge_interpolation_scheme
  )
  Jᵢ₊½, Jⱼ₊½, Jₖ₊½ = edge_functions_3d(
    forward_metrics.J, diff_backend, edge_interpolation_scheme
  )
  Jinv_ᵢ₊½, Jinv_ⱼ₊½, Jinv_ₖ₊½ = edge_functions_3d(
    inverse_metrics.Jinv, diff_backend, edge_interpolation_scheme
  )

  #! format: off
  edge_funcs = (;
    ξ̂xᵢ₊½, ξ̂xⱼ₊½, ξ̂xₖ₊½,
    η̂xᵢ₊½, η̂xⱼ₊½, η̂xₖ₊½,
    ζ̂xᵢ₊½, ζ̂xⱼ₊½, ζ̂xₖ₊½,
    ξ̂yᵢ₊½, ξ̂yⱼ₊½, ξ̂yₖ₊½,
    η̂yᵢ₊½, η̂yⱼ₊½, η̂yₖ₊½,
    ζ̂yᵢ₊½, ζ̂yⱼ₊½, ζ̂yₖ₊½,
    ξ̂zᵢ₊½, ξ̂zⱼ₊½, ξ̂zₖ₊½,
    η̂zᵢ₊½, η̂zⱼ₊½, η̂zₖ₊½,
    ζ̂zᵢ₊½, ζ̂zⱼ₊½, ζ̂zₖ₊½,
    Jinv_ᵢ₊½, Jinv_ⱼ₊½, Jinv_ₖ₊½,
    Jᵢ₊½, Jⱼ₊½, Jₖ₊½,
  )
   #! format: on

  return edge_funcs
end

function get_edge_functions_2d(
  forward_metrics,
  inverse_metrics,
  diff_backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=CurvatureCorrectedReconstruction(),
)
  Jinv_ᵢ₊½, Jinv_ⱼ₊½ = edge_functions_2d(
    inverse_metrics.Jinv, diff_backend, edge_interpolation_scheme
  )
  norm_Jinv_ᵢ₊½, norm_Jinv_ⱼ₊½ = edge_functions_2d(
    inverse_metrics.Jinv_norm, diff_backend, edge_interpolation_scheme
  )

  edge_funcs = (; Jinv_ᵢ₊½, Jinv_ⱼ₊½, norm_Jinv_ᵢ₊½, norm_Jinv_ⱼ₊½)

  return edge_funcs
end

function get_edge_functions_1d(
  forward_metrics,
  inverse_metrics,
  diff_backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=CurvatureCorrectedReconstruction(),
)
  Jinv_ᵢ₊½ = edge_functions_1d(
    inverse_metrics.Jinv, diff_backend, edge_interpolation_scheme
  )
  norm_Jinv_ᵢ₊½ = edge_functions_1d(
    inverse_metrics.Jinv_norm, diff_backend, edge_interpolation_scheme
  )

  edge_funcs = (; Jinv_ᵢ₊½, norm_Jinv_ᵢ₊½)

  return edge_funcs
end


function get_ad_thomas_lombard_inverse_metric_terms(x, y, z, backend)
  function jacobian_matrix(t, i, j, k, p)
    DifferentiationInterface.jacobian(
      u -> SVector(
        x(t, u[1], u[2], u[3], p),
        y(t, u[1], u[2], u[3], p),
        z(t, u[1], u[2], u[3], p),
      ),
      backend,
      SVector(i, j, k),
    )
  end

  inverse_jacobian_matrix(t, i, j, k, p) = inv(jacobian_matrix(t, i, j, k, p))

  zero_face_component(t, i, j, k, p) = zero(i + j + k)
  Jinvᵢ₊½(t, i, j, k, p) = inverse_jacobian_matrix(t, i + 0.5, j, k, p)
  Jinvⱼ₊½(t, i, j, k, p) = inverse_jacobian_matrix(t, i, j + 0.5, k, p)
  Jinvₖ₊½(t, i, j, k, p) = inverse_jacobian_matrix(t, i, j, k + 0.5, p)

  inverse_metrics = (; Jinv=inverse_jacobian_matrix)
  edge_metrics = (;
    ad_thomas_lombard_metric=true,
    diff_backend=backend,
    x,
    y,
    z,
    jacobian=jacobian_matrix,
    ξ̂xᵢ₊½=zero_face_component,
    η̂xᵢ₊½=zero_face_component,
    ζ̂xᵢ₊½=zero_face_component,
    ξ̂yᵢ₊½=zero_face_component,
    η̂yᵢ₊½=zero_face_component,
    ζ̂yᵢ₊½=zero_face_component,
    ξ̂zᵢ₊½=zero_face_component,
    η̂zᵢ₊½=zero_face_component,
    ζ̂zᵢ₊½=zero_face_component,
    Jinvᵢ₊½,
    #
    ξ̂xⱼ₊½=zero_face_component,
    η̂xⱼ₊½=zero_face_component,
    ζ̂xⱼ₊½=zero_face_component,
    ξ̂yⱼ₊½=zero_face_component,
    η̂yⱼ₊½=zero_face_component,
    ζ̂yⱼ₊½=zero_face_component,
    ξ̂zⱼ₊½=zero_face_component,
    η̂zⱼ₊½=zero_face_component,
    ζ̂zⱼ₊½=zero_face_component,
    Jinvⱼ₊½,
    #
    ξ̂xₖ₊½=zero_face_component,
    η̂xₖ₊½=zero_face_component,
    ζ̂xₖ₊½=zero_face_component,
    ξ̂yₖ₊½=zero_face_component,
    η̂yₖ₊½=zero_face_component,
    ζ̂yₖ₊½=zero_face_component,
    ξ̂zₖ₊½=zero_face_component,
    η̂zₖ₊½=zero_face_component,
    ζ̂zₖ₊½=zero_face_component,
    Jinvₖ₊½,
  )

  return inverse_metrics, edge_metrics
end


function get_inverse_metric_terms(
  x,
  y,
  z,
  backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=CurvatureCorrectedReconstruction(),
)
  edge_scheme = edge_interpolation_scheme

  if edge_scheme isa ADThomasLombardMetric
    return get_ad_thomas_lombard_inverse_metric_terms(x, y, z, backend)
  end

  #

  function get_jacobian_matrix(x, y, z, backend)
    function jacobian_matrix(t, i, j, k, p)
      DifferentiationInterface.jacobian(
        u -> SVector(x(t, u..., p), y(t, u..., p), z(t, u..., p)),
        backend,
        @SVector [i, j, k]
      )
    end

    return jacobian_matrix
  end

  function get_inverse_jacobian_matrix(x, y, z, backend)
    function inverse_jacobian_matrix(t, i, j, k, p)
      return inv(
        DifferentiationInterface.jacobian(
          u -> SVector(x(t, u..., p), y(t, u..., p), z(t, u..., p)),
          backend,
          @SVector [i, j, k]
        ),
      )
    end
  end

  function get_normalized_inverse_jacobian_matrix(x, y, z, backend)
    function normalized_inverse_jacobian_matrix(t, i, j, k, p)
      jac = DifferentiationInterface.jacobian(
        u -> SVector(x(t, u..., p), y(t, u..., p), z(t, u..., p)),
        backend,
        @SVector [i, j, k]
      )
      return inv(jac) .* det(jac)
    end
  end

  jacobian_matrix = get_jacobian_matrix(x, y, z, backend)
  inverse_jacobian_matrix = get_inverse_jacobian_matrix(x, y, z, backend)
  normalized_inverse_jacobian_matrix = get_normalized_inverse_jacobian_matrix(
    x, y, z, backend
  )

  xξ(t, i, j, k, p) = derivative(ξ -> x(t, ξ, j, k, p), backend, i)
  xη(t, i, j, k, p) = derivative(η -> x(t, i, η, k, p), backend, j)
  xζ(t, i, j, k, p) = derivative(ζ -> x(t, i, j, ζ, p), backend, k)
  yξ(t, i, j, k, p) = derivative(ξ -> y(t, ξ, j, k, p), backend, i)
  yη(t, i, j, k, p) = derivative(η -> y(t, i, η, k, p), backend, j)
  yζ(t, i, j, k, p) = derivative(ζ -> y(t, i, j, ζ, p), backend, k)
  zξ(t, i, j, k, p) = derivative(ξ -> z(t, ξ, j, k, p), backend, i)
  zη(t, i, j, k, p) = derivative(η -> z(t, i, η, k, p), backend, j)
  zζ(t, i, j, k, p) = derivative(ζ -> z(t, i, j, ζ, p), backend, k)

  x_ξ_y(t, i, j, k, p) = xξ(t, i, j, k, p) * y(t, i, j, k, p)
  x_η_y(t, i, j, k, p) = xη(t, i, j, k, p) * y(t, i, j, k, p)
  x_ζ_y(t, i, j, k, p) = xζ(t, i, j, k, p) * y(t, i, j, k, p)
  y_ξ_z(t, i, j, k, p) = yξ(t, i, j, k, p) * z(t, i, j, k, p)
  y_η_z(t, i, j, k, p) = yη(t, i, j, k, p) * z(t, i, j, k, p)
  y_ζ_z(t, i, j, k, p) = yζ(t, i, j, k, p) * z(t, i, j, k, p)
  z_ξ_x(t, i, j, k, p) = zξ(t, i, j, k, p) * x(t, i, j, k, p)
  z_η_x(t, i, j, k, p) = zη(t, i, j, k, p) * x(t, i, j, k, p)
  z_ζ_x(t, i, j, k, p) = zζ(t, i, j, k, p) * x(t, i, j, k, p)

  y_η_z_ζ = ∂ϕ_∂ζ_3d(y_η_z, backend, edge_scheme)
  y_ζ_z_η = ∂ϕ_∂η_3d(y_ζ_z, backend, edge_scheme)

  y_ζ_z_ξ = ∂ϕ_∂ξ_3d(y_ζ_z, backend, edge_scheme)
  y_ξ_z_ζ = ∂ϕ_∂ζ_3d(y_ξ_z, backend, edge_scheme)

  y_ξ_z_η = ∂ϕ_∂η_3d(y_ξ_z, backend, edge_scheme)
  y_η_z_ξ = ∂ϕ_∂ξ_3d(y_η_z, backend, edge_scheme)

  z_η_x_ζ = ∂ϕ_∂ζ_3d(z_η_x, backend, edge_scheme)
  z_ζ_x_η = ∂ϕ_∂η_3d(z_ζ_x, backend, edge_scheme)

  z_ζ_x_ξ = ∂ϕ_∂ξ_3d(z_ζ_x, backend, edge_scheme)
  z_ξ_x_ζ = ∂ϕ_∂ζ_3d(z_ξ_x, backend, edge_scheme)

  z_ξ_x_η = ∂ϕ_∂η_3d(z_ξ_x, backend, edge_scheme)
  z_η_x_ξ = ∂ϕ_∂ξ_3d(z_η_x, backend, edge_scheme)

  x_η_y_ζ = ∂ϕ_∂ζ_3d(x_η_y, backend, edge_scheme)
  x_ζ_y_η = ∂ϕ_∂η_3d(x_ζ_y, backend, edge_scheme)

  x_ζ_y_ξ = ∂ϕ_∂ξ_3d(x_ζ_y, backend, edge_scheme)
  x_ξ_y_ζ = ∂ϕ_∂ζ_3d(x_ξ_y, backend, edge_scheme)

  x_ξ_y_η = ∂ϕ_∂η_3d(x_ξ_y, backend, edge_scheme)
  x_η_y_ξ = ∂ϕ_∂ξ_3d(x_η_y, backend, edge_scheme)

  # Do NOT put eps() tolerance checks on these! It will create GCL-related errors
  ξ̂x(t, i, j, k, p) = y_η_z_ζ(t, i, j, k, p) − y_ζ_z_η(t, i, j, k, p)
  η̂x(t, i, j, k, p) = y_ζ_z_ξ(t, i, j, k, p) − y_ξ_z_ζ(t, i, j, k, p)
  ζ̂x(t, i, j, k, p) = y_ξ_z_η(t, i, j, k, p) − y_η_z_ξ(t, i, j, k, p)
  ξ̂y(t, i, j, k, p) = z_η_x_ζ(t, i, j, k, p) − z_ζ_x_η(t, i, j, k, p)
  η̂y(t, i, j, k, p) = z_ζ_x_ξ(t, i, j, k, p) − z_ξ_x_ζ(t, i, j, k, p)
  ζ̂y(t, i, j, k, p) = z_ξ_x_η(t, i, j, k, p) − z_η_x_ξ(t, i, j, k, p)
  ξ̂z(t, i, j, k, p) = x_η_y_ζ(t, i, j, k, p) − x_ζ_y_η(t, i, j, k, p)
  η̂z(t, i, j, k, p) = x_ζ_y_ξ(t, i, j, k, p) − x_ξ_y_ζ(t, i, j, k, p)
  ζ̂z(t, i, j, k, p) = x_ξ_y_η(t, i, j, k, p) − x_η_y_ξ(t, i, j, k, p)

  ξ̂x_val_and_ξderivs = ξ_derivs(ξ̂x, backend, edge_scheme)
  η̂x_val_and_ξderivs = ξ_derivs(η̂x, backend, edge_scheme)
  ζ̂x_val_and_ξderivs = ξ_derivs(ζ̂x, backend, edge_scheme)
  ξ̂y_val_and_ξderivs = ξ_derivs(ξ̂y, backend, edge_scheme)
  η̂y_val_and_ξderivs = ξ_derivs(η̂y, backend, edge_scheme)
  ζ̂y_val_and_ξderivs = ξ_derivs(ζ̂y, backend, edge_scheme)
  ξ̂z_val_and_ξderivs = ξ_derivs(ξ̂z, backend, edge_scheme)
  η̂z_val_and_ξderivs = ξ_derivs(η̂z, backend, edge_scheme)
  ζ̂z_val_and_ξderivs = ξ_derivs(ζ̂z, backend, edge_scheme)

  ξ̂x_val_and_ηderivs = η_derivs(ξ̂x, backend, edge_scheme)
  η̂x_val_and_ηderivs = η_derivs(η̂x, backend, edge_scheme)
  ζ̂x_val_and_ηderivs = η_derivs(ζ̂x, backend, edge_scheme)
  ξ̂y_val_and_ηderivs = η_derivs(ξ̂y, backend, edge_scheme)
  η̂y_val_and_ηderivs = η_derivs(η̂y, backend, edge_scheme)
  ζ̂y_val_and_ηderivs = η_derivs(ζ̂y, backend, edge_scheme)
  ξ̂z_val_and_ηderivs = η_derivs(ξ̂z, backend, edge_scheme)
  η̂z_val_and_ηderivs = η_derivs(η̂z, backend, edge_scheme)
  ζ̂z_val_and_ηderivs = η_derivs(ζ̂z, backend, edge_scheme)

  ξ̂x_val_and_ζderivs = ζ_derivs(ξ̂x, backend, edge_scheme)
  η̂x_val_and_ζderivs = ζ_derivs(η̂x, backend, edge_scheme)
  ζ̂x_val_and_ζderivs = ζ_derivs(ζ̂x, backend, edge_scheme)
  ξ̂y_val_and_ζderivs = ζ_derivs(ξ̂y, backend, edge_scheme)
  η̂y_val_and_ζderivs = ζ_derivs(η̂y, backend, edge_scheme)
  ζ̂y_val_and_ζderivs = ζ_derivs(ζ̂y, backend, edge_scheme)
  ξ̂z_val_and_ζderivs = ζ_derivs(ξ̂z, backend, edge_scheme)
  η̂z_val_and_ζderivs = ζ_derivs(η̂z, backend, edge_scheme)
  ζ̂z_val_and_ζderivs = ζ_derivs(ζ̂z, backend, edge_scheme)

  Jinv_val_and_ξderivs = ξ_derivs(inverse_jacobian_matrix, backend, edge_scheme)
  Jinv_val_and_ηderivs = η_derivs(inverse_jacobian_matrix, backend, edge_scheme)
  Jinv_val_and_ζderivs = ζ_derivs(inverse_jacobian_matrix, backend, edge_scheme)

  ξ̂xᵢ₊½(t, i, j, k, p) = ϕ_iedge(ξ̂x_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  η̂xᵢ₊½(t, i, j, k, p) = ϕ_iedge(η̂x_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  ζ̂xᵢ₊½(t, i, j, k, p) = ϕ_iedge(ζ̂x_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  ξ̂yᵢ₊½(t, i, j, k, p) = ϕ_iedge(ξ̂y_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  η̂yᵢ₊½(t, i, j, k, p) = ϕ_iedge(η̂y_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  ζ̂yᵢ₊½(t, i, j, k, p) = ϕ_iedge(ζ̂y_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  ξ̂zᵢ₊½(t, i, j, k, p) = ϕ_iedge(ξ̂z_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  η̂zᵢ₊½(t, i, j, k, p) = ϕ_iedge(η̂z_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  ζ̂zᵢ₊½(t, i, j, k, p) = ϕ_iedge(ζ̂z_val_and_ξderivs, t, i, j, k, p, edge_scheme)

  ξ̂xⱼ₊½(t, i, j, k, p) = ϕ_jedge(ξ̂x_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  η̂xⱼ₊½(t, i, j, k, p) = ϕ_jedge(η̂x_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  ζ̂xⱼ₊½(t, i, j, k, p) = ϕ_jedge(ζ̂x_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  ξ̂yⱼ₊½(t, i, j, k, p) = ϕ_jedge(ξ̂y_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  η̂yⱼ₊½(t, i, j, k, p) = ϕ_jedge(η̂y_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  ζ̂yⱼ₊½(t, i, j, k, p) = ϕ_jedge(ζ̂y_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  ξ̂zⱼ₊½(t, i, j, k, p) = ϕ_jedge(ξ̂z_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  η̂zⱼ₊½(t, i, j, k, p) = ϕ_jedge(η̂z_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  ζ̂zⱼ₊½(t, i, j, k, p) = ϕ_jedge(ζ̂z_val_and_ηderivs, t, i, j, k, p, edge_scheme)

  ξ̂xₖ₊½(t, i, j, k, p) = ϕ_kedge(ξ̂x_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  η̂xₖ₊½(t, i, j, k, p) = ϕ_kedge(η̂x_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  ζ̂xₖ₊½(t, i, j, k, p) = ϕ_kedge(ζ̂x_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  ξ̂yₖ₊½(t, i, j, k, p) = ϕ_kedge(ξ̂y_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  η̂yₖ₊½(t, i, j, k, p) = ϕ_kedge(η̂y_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  ζ̂yₖ₊½(t, i, j, k, p) = ϕ_kedge(ζ̂y_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  ξ̂zₖ₊½(t, i, j, k, p) = ϕ_kedge(ξ̂z_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  η̂zₖ₊½(t, i, j, k, p) = ϕ_kedge(η̂z_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  ζ̂zₖ₊½(t, i, j, k, p) = ϕ_kedge(ζ̂z_val_and_ζderivs, t, i, j, k, p, edge_scheme)

  Jinvᵢ₊½(t, i, j, k, p) = ϕ_iedge(Jinv_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  Jinvⱼ₊½(t, i, j, k, p) = ϕ_jedge(Jinv_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  Jinvₖ₊½(t, i, j, k, p) = ϕ_kedge(Jinv_val_and_ζderivs, t, i, j, k, p, edge_scheme)

  return (
    (;
      ξ̂x,
      η̂x,
      ζ̂x,
      ξ̂y,
      η̂y,
      ζ̂y,
      ξ̂z,
      η̂z,
      ζ̂z,
      # J=jacobian_matrix,
      Jinv=inverse_jacobian_matrix,
      # normJinv=normalized_inverse_jacobian_matrix,
    ),
    (;
      ξ̂xᵢ₊½,
      η̂xᵢ₊½,
      ζ̂xᵢ₊½,
      ξ̂yᵢ₊½,
      η̂yᵢ₊½,
      ζ̂yᵢ₊½,
      ξ̂zᵢ₊½,
      η̂zᵢ₊½,
      ζ̂zᵢ₊½,
      Jinvᵢ₊½,
      #
      ξ̂xⱼ₊½,
      η̂xⱼ₊½,
      ζ̂xⱼ₊½,
      ξ̂yⱼ₊½,
      η̂yⱼ₊½,
      ζ̂yⱼ₊½,
      ξ̂zⱼ₊½,
      η̂zⱼ₊½,
      ζ̂zⱼ₊½,
      Jinvⱼ₊½,
      #
      ξ̂xₖ₊½,
      η̂xₖ₊½,
      ζ̂xₖ₊½,
      ξ̂yₖ₊½,
      η̂yₖ₊½,
      ζ̂yₖ₊½,
      ξ̂zₖ₊½,
      η̂zₖ₊½,
      ζ̂zₖ₊½,
      Jinvₖ₊½,
    ),
  )
end
#-------------------------------------------------------------
#-------------------------------------------------------------

@inline _ξ_eval_3d_metriccache(ξ, ϕ, t, j, k, p) = ϕ(t, ξ, j, k, p)
@inline _η_eval_3d_metriccache(η, ϕ, t, i, k, p) = ϕ(t, i, η, k, p)
@inline _ζ_eval_3d_metriccache(ζ, ϕ, t, i, j, p) = ϕ(t, i, j, ζ, p)

ξ_derivs(ϕ, backend) = ξ_derivs(ϕ, backend, CurvatureCorrectedReconstruction())
η_derivs(ϕ, backend) = η_derivs(ϕ, backend, CurvatureCorrectedReconstruction())
ζ_derivs(ϕ, backend) = ζ_derivs(ϕ, backend, CurvatureCorrectedReconstruction())

function ξ_derivs(ϕ, backend, ::EndpointAverageReconstruction)
  ϕval(t, i, j, k, p) = _ξ_eval_3d_metriccache(i, ϕ, t, j, k, p)
  return ϕval
end

function ξ_derivs(ϕ, backend, ::GradientCorrectedReconstruction)
  cϕ = DifferentiationInterface.Constant(ϕ)
  prep = prepare_derivative(
    _ξ_eval_3d_metriccache,
    backend,
    0.0,
    cϕ,
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(nothing);
    strict=Val(false),
  )

  function ϕall(t, i, j, k, p)
    value_and_derivative(
      _ξ_eval_3d_metriccache,
      prep,
      backend,
      i,
      cϕ,
      DifferentiationInterface.Constant(t),
      DifferentiationInterface.Constant(j),
      DifferentiationInterface.Constant(k),
      DifferentiationInterface.Constant(p),
    )
  end

  return ϕall
end

function ξ_derivs(ϕ, backend, ::CurvatureCorrectedReconstruction)
  cϕ = DifferentiationInterface.Constant(ϕ)
  prep = prepare_second_derivative(
    _ξ_eval_3d_metriccache,
    backend,
    0.0,
    cϕ,
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(nothing);
    strict=Val(false),
  )

  function ϕall(t, i, j, k, p)
    value_derivative_and_second_derivative(
      _ξ_eval_3d_metriccache,
      prep,
      backend,
      i,
      cϕ,
      DifferentiationInterface.Constant(t),
      DifferentiationInterface.Constant(j),
      DifferentiationInterface.Constant(k),
      DifferentiationInterface.Constant(p),
    )
  end

  return ϕall
end

function η_derivs(ϕ, backend, ::EndpointAverageReconstruction)
  ϕval(t, i, j, k, p) = _η_eval_3d_metriccache(j, ϕ, t, i, k, p)
  return ϕval
end

function η_derivs(ϕ, backend, ::GradientCorrectedReconstruction)
  cϕ = DifferentiationInterface.Constant(ϕ)
  prep = prepare_derivative(
    _η_eval_3d_metriccache,
    backend,
    0.0,
    cϕ,
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(nothing);
    strict=Val(false),
  )

  function ϕall(t, i, j, k, p)
    value_and_derivative(
      _η_eval_3d_metriccache,
      prep,
      backend,
      j,
      cϕ,
      DifferentiationInterface.Constant(t),
      DifferentiationInterface.Constant(i),
      DifferentiationInterface.Constant(k),
      DifferentiationInterface.Constant(p),
    )
  end

  return ϕall
end

function η_derivs(ϕ, backend, ::CurvatureCorrectedReconstruction)
  cϕ = DifferentiationInterface.Constant(ϕ)
  prep = prepare_second_derivative(
    _η_eval_3d_metriccache,
    backend,
    0.0,
    cϕ,
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(nothing);
    strict=Val(false),
  )

  function ϕall(t, i, j, k, p)
    value_derivative_and_second_derivative(
      _η_eval_3d_metriccache,
      prep,
      backend,
      j,
      cϕ,
      DifferentiationInterface.Constant(t),
      DifferentiationInterface.Constant(i),
      DifferentiationInterface.Constant(k),
      DifferentiationInterface.Constant(p),
    )
  end

  return ϕall
end

function ζ_derivs(ϕ, backend, ::EndpointAverageReconstruction)
  ϕval(t, i, j, k, p) = _ζ_eval_3d_metriccache(k, ϕ, t, i, j, p)
  return ϕval
end

function ζ_derivs(ϕ, backend, ::GradientCorrectedReconstruction)
  cϕ = DifferentiationInterface.Constant(ϕ)
  prep = prepare_derivative(
    _ζ_eval_3d_metriccache,
    backend,
    0.0,
    cϕ,
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(nothing);
    strict=Val(false),
  )

  function ϕall(t, i, j, k, p)
    value_and_derivative(
      _ζ_eval_3d_metriccache,
      prep,
      backend,
      k,
      cϕ,
      DifferentiationInterface.Constant(t),
      DifferentiationInterface.Constant(i),
      DifferentiationInterface.Constant(j),
      DifferentiationInterface.Constant(p),
    )
  end

  return ϕall
end

function ζ_derivs(ϕ, backend, ::CurvatureCorrectedReconstruction)
  cϕ = DifferentiationInterface.Constant(ϕ)
  prep = prepare_second_derivative(
    _ζ_eval_3d_metriccache,
    backend,
    0.0,
    cϕ,
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(nothing);
    strict=Val(false),
  )

  function ϕall(t, i, j, k, p)
    value_derivative_and_second_derivative(
      _ζ_eval_3d_metriccache,
      prep,
      backend,
      k,
      cϕ,
      DifferentiationInterface.Constant(t),
      DifferentiationInterface.Constant(i),
      DifferentiationInterface.Constant(j),
      DifferentiationInterface.Constant(p),
    )
  end

  return ϕall
end


# Edge reconstruction helpers used by metric caches.

#-------------------------------------------------------------
# 1D edge reconstruction
#-------------------------------------------------------------
function ϕ_iedge(ϕ_eval, t, i::Real, p::NamedTuple)
  ϕ_iedge(ϕ_eval, t, i, p, CurvatureCorrectedReconstruction(), 1)
end

function ϕ_iedge(
  ϕ_eval, t, i::Real, p::NamedTuple, edge_interpolation_scheme::EdgeInterpolationSchemeTrait
)
  ϕ_iedge(ϕ_eval, t, i, p, edge_interpolation_scheme, 1)
end

function ϕ_iedge(
  ϕ_eval,
  t,
  i::Real,
  p::NamedTuple,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
  Δξ::Real,
)
  ϕ_iedge(ϕ_eval, t, i, p, edge_interpolation_scheme, Δξ)
end

function ϕ_iedge(ϕ_values, t, i::Real, p::NamedTuple, ::EndpointAverageReconstruction, Δξ::Real)
  return _edge_reconstruct(
    ϕ_values(t, i, p), ϕ_values(t, i + 1, p), EndpointAverageReconstruction()
  )
end

function ϕ_iedge(
  ϕ_val_and_derivs, t, i::Real, p::NamedTuple, ::GradientCorrectedReconstruction, Δξ::Real
)
  ϕᵢ, ϕξᵢ = ϕ_val_and_derivs(t, i, p)
  ϕᵢ₊₁, ϕξᵢ₊₁ = ϕ_val_and_derivs(t, i + 1, p)
  return _edge_reconstruct(ϕᵢ, ϕξᵢ, ϕᵢ₊₁, ϕξᵢ₊₁, GradientCorrectedReconstruction(), Δξ)
end

function ϕ_iedge(
  ϕ_val_and_derivs, t, i::Real, p::NamedTuple, ::CurvatureCorrectedReconstruction, Δξ::Real
)
  ϕᵢ, ϕξᵢ, ϕξξᵢ = ϕ_val_and_derivs(t, i, p)
  ϕᵢ₊₁, ϕξᵢ₊₁, ϕξξᵢ₊₁ = ϕ_val_and_derivs(t, i + 1, p)
  return _edge_reconstruct(
    ϕᵢ, ϕξᵢ, ϕξξᵢ, ϕᵢ₊₁, ϕξᵢ₊₁, ϕξξᵢ₊₁, CurvatureCorrectedReconstruction(), Δξ
  )
end

#-------------------------------------------------------------
# 2D edge reconstruction
#-------------------------------------------------------------
function ϕ_iedge(ϕ_eval, t, i::Real, j::Real, p::NamedTuple)
  ϕ_iedge(ϕ_eval, t, i, j, p, CurvatureCorrectedReconstruction(), 1)
end

function ϕ_iedge(
  ϕ_eval,
  t,
  i::Real,
  j::Real,
  p::NamedTuple,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
)
  ϕ_iedge(ϕ_eval, t, i, j, p, edge_interpolation_scheme, 1)
end

function ϕ_iedge(
  ϕ_eval,
  t,
  i::Real,
  j::Real,
  p::NamedTuple,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
  Δξ::Real,
)
  ϕ_iedge(ϕ_eval, t, i, j, p, edge_interpolation_scheme, Δξ)
end

function ϕ_iedge(
  ϕ_values, t, i::Real, j::Real, p::NamedTuple, ::EndpointAverageReconstruction, Δξ::Real
)
  return _edge_reconstruct(
    ϕ_values(t, i, j, p), ϕ_values(t, i + 1, j, p), EndpointAverageReconstruction()
  )
end

function ϕ_iedge(
  ϕ_val_and_derivs, t, i::Real, j::Real, p::NamedTuple, ::GradientCorrectedReconstruction, Δξ::Real
)
  ϕᵢ, ϕξᵢ = ϕ_val_and_derivs(t, i, j, p)
  ϕᵢ₊₁, ϕξᵢ₊₁ = ϕ_val_and_derivs(t, i + 1, j, p)
  return _edge_reconstruct(ϕᵢ, ϕξᵢ, ϕᵢ₊₁, ϕξᵢ₊₁, GradientCorrectedReconstruction(), Δξ)
end

function ϕ_iedge(
  ϕ_val_and_derivs, t, i::Real, j::Real, p::NamedTuple, ::CurvatureCorrectedReconstruction, Δξ::Real
)
  ϕᵢ, ϕξᵢ, ϕξξᵢ = ϕ_val_and_derivs(t, i, j, p)
  ϕᵢ₊₁, ϕξᵢ₊₁, ϕξξᵢ₊₁ = ϕ_val_and_derivs(t, i + 1, j, p)
  return _edge_reconstruct(
    ϕᵢ, ϕξᵢ, ϕξξᵢ, ϕᵢ₊₁, ϕξᵢ₊₁, ϕξξᵢ₊₁, CurvatureCorrectedReconstruction(), Δξ
  )
end

function ϕ_jedge(ϕ_eval, t, i::Real, j::Real, p::NamedTuple)
  ϕ_jedge(ϕ_eval, t, i, j, p, CurvatureCorrectedReconstruction(), 1)
end

function ϕ_jedge(
  ϕ_eval,
  t,
  i::Real,
  j::Real,
  p::NamedTuple,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
)
  ϕ_jedge(ϕ_eval, t, i, j, p, edge_interpolation_scheme, 1)
end

function ϕ_jedge(
  ϕ_eval,
  t,
  i::Real,
  j::Real,
  p::NamedTuple,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
  Δξ::Real,
)
  ϕ_jedge(ϕ_eval, t, i, j, p, edge_interpolation_scheme, Δξ)
end

function ϕ_jedge(
  ϕ_values, t, i::Real, j::Real, p::NamedTuple, ::EndpointAverageReconstruction, Δξ::Real
)
  return _edge_reconstruct(
    ϕ_values(t, i, j, p), ϕ_values(t, i, j + 1, p), EndpointAverageReconstruction()
  )
end

function ϕ_jedge(
  ϕ_val_and_derivs, t, i::Real, j::Real, p::NamedTuple, ::GradientCorrectedReconstruction, Δξ::Real
)
  ϕⱼ, ϕηⱼ = ϕ_val_and_derivs(t, i, j, p)
  ϕⱼ₊₁, ϕηⱼ₊₁ = ϕ_val_and_derivs(t, i, j + 1, p)
  return _edge_reconstruct(ϕⱼ, ϕηⱼ, ϕⱼ₊₁, ϕηⱼ₊₁, GradientCorrectedReconstruction(), Δξ)
end

function ϕ_jedge(
  ϕ_val_and_derivs, t, i::Real, j::Real, p::NamedTuple, ::CurvatureCorrectedReconstruction, Δξ::Real
)
  ϕⱼ, ϕηⱼ, ϕηηⱼ = ϕ_val_and_derivs(t, i, j, p)
  ϕⱼ₊₁, ϕηⱼ₊₁, ϕηηⱼ₊₁ = ϕ_val_and_derivs(t, i, j + 1, p)
  return _edge_reconstruct(
    ϕⱼ, ϕηⱼ, ϕηηⱼ, ϕⱼ₊₁, ϕηⱼ₊₁, ϕηηⱼ₊₁, CurvatureCorrectedReconstruction(), Δξ
  )
end

#-------------------------------------------------------------
# 3D edge reconstruction
#-------------------------------------------------------------
function ϕ_iedge(ϕ_eval, t, i::Real, j::Real, k::Real, p::NamedTuple)
  ϕ_iedge(ϕ_eval, t, i, j, k, p, CurvatureCorrectedReconstruction(), 1)
end

function ϕ_iedge(
  ϕ_eval,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
)
  ϕ_iedge(ϕ_eval, t, i, j, k, p, edge_interpolation_scheme, 1)
end

function ϕ_iedge(
  ϕ_eval,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
  Δξ::Real,
)
  ϕ_iedge(ϕ_eval, t, i, j, k, p, edge_interpolation_scheme, Δξ)
end

function ϕ_iedge(
  ϕ_values, t, i::Real, j::Real, k::Real, p::NamedTuple, ::EndpointAverageReconstruction, Δξ::Real
)
  return _edge_reconstruct(
    ϕ_values(t, i, j, k, p), ϕ_values(t, i + 1, j, k, p), EndpointAverageReconstruction()
  )
end

function ϕ_iedge(
  ϕ_val_and_derivs,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  ::GradientCorrectedReconstruction,
  Δξ::Real,
)
  ϕᵢ, ϕξᵢ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕᵢ₊₁, ϕξᵢ₊₁ = ϕ_val_and_derivs(t, i + 1, j, k, p)
  return _edge_reconstruct(ϕᵢ, ϕξᵢ, ϕᵢ₊₁, ϕξᵢ₊₁, GradientCorrectedReconstruction(), Δξ)
end

function ϕ_iedge(
  ϕ_val_and_derivs,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  ::CurvatureCorrectedReconstruction,
  Δξ::Real,
)
  ϕᵢ, ϕξᵢ, ϕξξᵢ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕᵢ₊₁, ϕξᵢ₊₁, ϕξξᵢ₊₁ = ϕ_val_and_derivs(t, i + 1, j, k, p)
  return _edge_reconstruct(
    ϕᵢ, ϕξᵢ, ϕξξᵢ, ϕᵢ₊₁, ϕξᵢ₊₁, ϕξξᵢ₊₁, CurvatureCorrectedReconstruction(), Δξ
  )
end

function ϕ_jedge(ϕ_eval, t, i::Real, j::Real, k::Real, p::NamedTuple)
  ϕ_jedge(ϕ_eval, t, i, j, k, p, CurvatureCorrectedReconstruction(), 1)
end

function ϕ_jedge(
  ϕ_eval,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
)
  ϕ_jedge(ϕ_eval, t, i, j, k, p, edge_interpolation_scheme, 1)
end

function ϕ_jedge(
  ϕ_eval,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
  Δξ::Real,
)
  ϕ_jedge(ϕ_eval, t, i, j, k, p, edge_interpolation_scheme, Δξ)
end

function ϕ_jedge(
  ϕ_values, t, i::Real, j::Real, k::Real, p::NamedTuple, ::EndpointAverageReconstruction, Δξ::Real
)
  return _edge_reconstruct(
    ϕ_values(t, i, j, k, p), ϕ_values(t, i, j + 1, k, p), EndpointAverageReconstruction()
  )
end

function ϕ_jedge(
  ϕ_val_and_derivs,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  ::GradientCorrectedReconstruction,
  Δξ::Real,
)
  ϕⱼ, ϕηⱼ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕⱼ₊₁, ϕηⱼ₊₁ = ϕ_val_and_derivs(t, i, j + 1, k, p)
  return _edge_reconstruct(ϕⱼ, ϕηⱼ, ϕⱼ₊₁, ϕηⱼ₊₁, GradientCorrectedReconstruction(), Δξ)
end

function ϕ_jedge(
  ϕ_val_and_derivs,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  ::CurvatureCorrectedReconstruction,
  Δξ::Real,
)
  ϕⱼ, ϕηⱼ, ϕηηⱼ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕⱼ₊₁, ϕηⱼ₊₁, ϕηηⱼ₊₁ = ϕ_val_and_derivs(t, i, j + 1, k, p)
  return _edge_reconstruct(
    ϕⱼ, ϕηⱼ, ϕηηⱼ, ϕⱼ₊₁, ϕηⱼ₊₁, ϕηηⱼ₊₁, CurvatureCorrectedReconstruction(), Δξ
  )
end

function ϕ_kedge(ϕ_eval, t, i::Real, j::Real, k::Real, p::NamedTuple)
  ϕ_kedge(ϕ_eval, t, i, j, k, p, CurvatureCorrectedReconstruction(), 1)
end

function ϕ_kedge(
  ϕ_eval,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
)
  ϕ_kedge(ϕ_eval, t, i, j, k, p, edge_interpolation_scheme, 1)
end

function ϕ_kedge(
  ϕ_eval,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
  Δξ::Real,
)
  ϕ_kedge(ϕ_eval, t, i, j, k, p, edge_interpolation_scheme, Δξ)
end

function ϕ_kedge(
  ϕ_values, t, i::Real, j::Real, k::Real, p::NamedTuple, ::EndpointAverageReconstruction, Δξ::Real
)
  return _edge_reconstruct(
    ϕ_values(t, i, j, k, p), ϕ_values(t, i, j, k + 1, p), EndpointAverageReconstruction()
  )
end

function ϕ_kedge(
  ϕ_val_and_derivs,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  ::GradientCorrectedReconstruction,
  Δξ::Real,
)
  ϕₖ, ϕζₖ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕₖ₊₁, ϕζₖ₊₁ = ϕ_val_and_derivs(t, i, j, k + 1, p)
  return _edge_reconstruct(ϕₖ, ϕζₖ, ϕₖ₊₁, ϕζₖ₊₁, GradientCorrectedReconstruction(), Δξ)
end

function ϕ_kedge(
  ϕ_val_and_derivs,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  ::CurvatureCorrectedReconstruction,
  Δξ::Real,
)
  ϕₖ, ϕζₖ, ϕζζₖ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕₖ₊₁, ϕζₖ₊₁, ϕζζₖ₊₁ = ϕ_val_and_derivs(t, i, j, k + 1, p)
  return _edge_reconstruct(
    ϕₖ, ϕζₖ, ϕζζₖ, ϕₖ₊₁, ϕζₖ₊₁, ϕζζₖ₊₁, CurvatureCorrectedReconstruction(), Δξ
  )
end


# Cell-center derivative helpers.

#-------------------------------------------------------------
# 3D cell-center derivatives
#-------------------------------------------------------------
∂ϕ_∂ξ_3d(ϕ, backend) = ∂ϕ_∂ξ_3d(ϕ, backend, CurvatureCorrectedReconstruction())

function ∂ϕ_∂ξ_3d(ϕ, backend, ::EndpointAverageReconstruction)
  ϕᵢ₊½(t, i, j, k, p) = _edge_reconstruct(
    ϕ(t, i, j, k, p), ϕ(t, i + 1, j, k, p), EndpointAverageReconstruction()
  )
  ∂ϕ_∂ξ(t, i, j, k, p) = ϕᵢ₊½(t, i, j, k, p) - ϕᵢ₊½(t, i - 1, j, k, p)
  return ∂ϕ_∂ξ
end

function ∂ϕ_∂ξ_3d(ϕ, backend, ::GradientCorrectedReconstruction)
  ξ_derivs(t, i, j, k, p) = value_and_derivative(ξ -> ϕ(t, ξ, j, k, p), backend, i)

  function ϕᵢ₊½(t, i, j, k, p)
    ϕᵢ, ∂ϕ_∂ξᵢ = ξ_derivs(t, i, j, k, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁ = ξ_derivs(t, i + 1, j, k, p)
    return _edge_reconstruct(ϕᵢ, ∂ϕ_∂ξᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, GradientCorrectedReconstruction())
  end

  ∂ϕ_∂ξ(t, i, j, k, p) = ϕᵢ₊½(t, i, j, k, p) - ϕᵢ₊½(t, i - 1, j, k, p)
  return ∂ϕ_∂ξ
end

function ∂ϕ_∂ξ_3d(ϕ, backend, ::CurvatureCorrectedReconstruction)
  ξ_derivs(t, i, j, k, p) = value_derivative_and_second_derivative(
    ξ -> ϕ(t, ξ, j, k, p), backend, i
  )

  function ϕᵢ₊½(t, i, j, k, p)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(t, i, j, k, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(t, i + 1, j, k, p)
    return _edge_reconstruct(
      ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁, CurvatureCorrectedReconstruction()
    )
  end

  ∂ϕ_∂ξ(t, i, j, k, p) = ϕᵢ₊½(t, i, j, k, p) - ϕᵢ₊½(t, i - 1, j, k, p)
  return ∂ϕ_∂ξ
end

∂ϕ_∂η_3d(ϕ, backend) = ∂ϕ_∂η_3d(ϕ, backend, CurvatureCorrectedReconstruction())

function ∂ϕ_∂η_3d(ϕ, backend, ::EndpointAverageReconstruction)
  ϕⱼ₊½(t, i, j, k, p) = _edge_reconstruct(
    ϕ(t, i, j, k, p), ϕ(t, i, j + 1, k, p), EndpointAverageReconstruction()
  )
  ∂ϕ_∂η(t, i, j, k, p) = ϕⱼ₊½(t, i, j, k, p) - ϕⱼ₊½(t, i, j - 1, k, p)
  return ∂ϕ_∂η
end

function ∂ϕ_∂η_3d(ϕ, backend, ::GradientCorrectedReconstruction)
  η_derivs(t, i, j, k, p) = value_and_derivative(η -> ϕ(t, i, η, k, p), backend, j)

  function ϕⱼ₊½(t, i, j, k, p)
    ϕⱼ, ∂ϕ_∂ηⱼ = η_derivs(t, i, j, k, p)
    ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁ = η_derivs(t, i, j + 1, k, p)
    return _edge_reconstruct(ϕⱼ, ∂ϕ_∂ηⱼ, ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, GradientCorrectedReconstruction())
  end

  ∂ϕ_∂η(t, i, j, k, p) = ϕⱼ₊½(t, i, j, k, p) - ϕⱼ₊½(t, i, j - 1, k, p)
  return ∂ϕ_∂η
end

function ∂ϕ_∂η_3d(ϕ, backend, ::CurvatureCorrectedReconstruction)
  η_derivs(t, i, j, k, p) = value_derivative_and_second_derivative(
    η -> ϕ(t, i, η, k, p), backend, j
  )

  function ϕⱼ₊½(t, i, j, k, p)
    ϕⱼ, ∂ϕ_∂ηⱼ, ∂²ϕ_∂η²ⱼ = η_derivs(t, i, j, k, p)
    ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, ∂²ϕ_∂η²ⱼ₊₁ = η_derivs(t, i, j + 1, k, p)
    return _edge_reconstruct(
      ϕⱼ, ∂ϕ_∂ηⱼ, ∂²ϕ_∂η²ⱼ, ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, ∂²ϕ_∂η²ⱼ₊₁, CurvatureCorrectedReconstruction()
    )
  end

  ∂ϕ_∂η(t, i, j, k, p) = ϕⱼ₊½(t, i, j, k, p) - ϕⱼ₊½(t, i, j - 1, k, p)
  return ∂ϕ_∂η
end

∂ϕ_∂ζ_3d(ϕ, backend) = ∂ϕ_∂ζ_3d(ϕ, backend, CurvatureCorrectedReconstruction())

function ∂ϕ_∂ζ_3d(ϕ, backend, ::EndpointAverageReconstruction)
  ϕₖ₊½(t, i, j, k, p) = _edge_reconstruct(
    ϕ(t, i, j, k, p), ϕ(t, i, j, k + 1, p), EndpointAverageReconstruction()
  )
  ∂ϕ_∂ζ(t, i, j, k, p) = ϕₖ₊½(t, i, j, k, p) - ϕₖ₊½(t, i, j, k - 1, p)
  return ∂ϕ_∂ζ
end

function ∂ϕ_∂ζ_3d(ϕ, backend, ::GradientCorrectedReconstruction)
  ζ_derivs(t, i, j, k, p) = value_and_derivative(ζ -> ϕ(t, i, j, ζ, p), backend, k)

  function ϕₖ₊½(t, i, j, k, p)
    ϕₖ, ∂ϕ_∂ζₖ = ζ_derivs(t, i, j, k, p)
    ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁ = ζ_derivs(t, i, j, k + 1, p)
    return _edge_reconstruct(ϕₖ, ∂ϕ_∂ζₖ, ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁, GradientCorrectedReconstruction())
  end

  ∂ϕ_∂ζ(t, i, j, k, p) = ϕₖ₊½(t, i, j, k, p) - ϕₖ₊½(t, i, j, k - 1, p)
  return ∂ϕ_∂ζ
end

function ∂ϕ_∂ζ_3d(ϕ, backend, ::CurvatureCorrectedReconstruction)
  ζ_derivs(t, i, j, k, p) = value_derivative_and_second_derivative(
    ζ -> ϕ(t, i, j, ζ, p), backend, k
  )

  function ϕₖ₊½(t, i, j, k, p)
    ϕₖ, ∂ϕ_∂ζₖ, ∂²ϕ_∂ζ²ₖ = ζ_derivs(t, i, j, k, p)
    ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁, ∂²ϕ_∂ζ²ₖ₊₁ = ζ_derivs(t, i, j, k + 1, p)
    return _edge_reconstruct(
      ϕₖ, ∂ϕ_∂ζₖ, ∂²ϕ_∂ζ²ₖ, ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁, ∂²ϕ_∂ζ²ₖ₊₁, CurvatureCorrectedReconstruction()
    )
  end

  ∂ϕ_∂ζ(t, i, j, k, p) = ϕₖ₊½(t, i, j, k, p) - ϕₖ₊½(t, i, j, k - 1, p)
  return ∂ϕ_∂ζ
end


# Metric cache construction for mapping callbacks.

@inline _default_direct_metric_scheme() = CurvatureCorrectedReconstruction()

@inline function _normalize_conserved_metric_scheme(
  dim::Val, scheme::EdgeInterpolationSchemeTrait
)
  return scheme
end

@inline function _normalize_conserved_metric_scheme(
  dim::Val{3}, scheme::ADThomasLombardMetric
)
  return scheme
end

function _normalize_conserved_metric_scheme(
  dim::Val{N}, ::ADThomasLombardMetric
) where {N}
  @warn "ADThomasLombardMetric simplifies to the direct metric form in 1-D/2-D." dimension = N
  return _default_direct_metric_scheme()
end

@inline function _metric_cache_for_mapping(
  dim::Val{1},
  mapping_functions,
  diff_backend,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
)
  scheme = _normalize_conserved_metric_scheme(dim, edge_interpolation_scheme)
  MetricCache(
    mapping_functions.x1, diff_backend; edge_interpolation_scheme=scheme
  )
end
@inline function _metric_cache_for_mapping(
  dim::Val{2},
  mapping_functions,
  diff_backend,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
)
  scheme = _normalize_conserved_metric_scheme(dim, edge_interpolation_scheme)
  MetricCache(
    mapping_functions.x1,
    mapping_functions.x2,
    diff_backend;
    edge_interpolation_scheme=scheme,
  )
end
@inline function _metric_cache_for_mapping(
  dim::Val{3},
  mapping_functions,
  diff_backend,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
)
  scheme = _normalize_conserved_metric_scheme(dim, edge_interpolation_scheme)
  MetricCache(
    mapping_functions.x1,
    mapping_functions.x2,
    mapping_functions.x3,
    diff_backend;
    edge_interpolation_scheme=scheme,
  )
end



# Metric storage kernels and launchers.

@inline function _inverse_and_normalized_edge_metrics(
  ::Val{1}, edge, edge_axis::Int, t, ξηζ, params
)
  if edge_axis != 1
    throw(ArgumentError("Invalid 1D face axis: $edge_axis"))
  end

  jinv_edge = edge.Jinv_ᵢ₊½
  norm_edge = edge.norm_Jinv_ᵢ₊½

  jinv_fun = jinv_edge isa NamedTuple ? jinv_edge.ϕᵢ₊½ : jinv_edge
  norm_fun = norm_edge isa NamedTuple ? norm_edge.ϕᵢ₊½ : norm_edge

  G = _as_smatrix(Val(1), jinv_fun(t, ξηζ..., params))
  Ghat = _as_smatrix(Val(1), norm_fun(t, ξηζ..., params))
  Jinv = det(G)
  return G, Ghat, Jinv
end

@inline function _inverse_and_normalized_edge_metrics(
  ::Val{2}, edge, edge_axis::Int, t, ξηζ, params
)
  if edge_axis == 1
    G = _as_smatrix(Val(2), edge.Jinv_ᵢ₊½(t, ξηζ..., params))
    Ghat = _as_smatrix(Val(2), edge.norm_Jinv_ᵢ₊½(t, ξηζ..., params))
  elseif edge_axis == 2
    G = _as_smatrix(Val(2), edge.Jinv_ⱼ₊½(t, ξηζ..., params))
    Ghat = _as_smatrix(Val(2), edge.norm_Jinv_ⱼ₊½(t, ξηζ..., params))
  else
    throw(ArgumentError("Invalid 2D face axis: $edge_axis"))
  end
  Jinv = det(G)
  return G, Ghat, Jinv
end

@inline function _inverse_and_normalized_edge_metrics(
  ::Val{3}, edge, edge_axis::Int, t, ξηζ, params
)
  if edge_axis == 1
    G = _as_smatrix(Val(3), edge.Jinvᵢ₊½(t, ξηζ..., params))
    Ghat = @SMatrix [
      edge.ξ̂xᵢ₊½(t, ξηζ..., params) edge.ξ̂yᵢ₊½(t, ξηζ..., params) edge.ξ̂zᵢ₊½(t, ξηζ..., params)
      edge.η̂xᵢ₊½(t, ξηζ..., params) edge.η̂yᵢ₊½(t, ξηζ..., params) edge.η̂zᵢ₊½(t, ξηζ..., params)
      edge.ζ̂xᵢ₊½(t, ξηζ..., params) edge.ζ̂yᵢ₊½(t, ξηζ..., params) edge.ζ̂zᵢ₊½(t, ξηζ..., params)
    ]
  elseif edge_axis == 2
    G = _as_smatrix(Val(3), edge.Jinvⱼ₊½(t, ξηζ..., params))
    Ghat = @SMatrix [
      edge.ξ̂xⱼ₊½(t, ξηζ..., params) edge.ξ̂yⱼ₊½(t, ξηζ..., params) edge.ξ̂zⱼ₊½(t, ξηζ..., params)
      edge.η̂xⱼ₊½(t, ξηζ..., params) edge.η̂yⱼ₊½(t, ξηζ..., params) edge.η̂zⱼ₊½(t, ξηζ..., params)
      edge.ζ̂xⱼ₊½(t, ξηζ..., params) edge.ζ̂yⱼ₊½(t, ξηζ..., params) edge.ζ̂zⱼ₊½(t, ξηζ..., params)
    ]
  elseif edge_axis == 3
    G = _as_smatrix(Val(3), edge.Jinvₖ₊½(t, ξηζ..., params))
    Ghat = @SMatrix [
      edge.ξ̂xₖ₊½(t, ξηζ..., params) edge.ξ̂yₖ₊½(t, ξηζ..., params) edge.ξ̂zₖ₊½(t, ξηζ..., params)
      edge.η̂xₖ₊½(t, ξηζ..., params) edge.η̂yₖ₊½(t, ξηζ..., params) edge.η̂zₖ₊½(t, ξηζ..., params)
      edge.ζ̂xₖ₊½(t, ξηζ..., params) edge.ζ̂yₖ₊½(t, ξηζ..., params) edge.ζ̂zₖ₊½(t, ξηζ..., params)
    ]
  else
    throw(ArgumentError("Invalid 3D face axis: $edge_axis"))
  end

  Jinv = det(G)
  return G, Ghat, Jinv
end

@kernel function _fill_cell_metric_storage_1d_kernel!(
  forward,
  inverse,
  jacobian,
  inverse_jacobian,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half

  F = _as_smatrix(Val(1), jacobian(t, ξ, params))
  G = _as_smatrix(Val(1), inverse_jacobian(t, ξ, params))
  J = det(F)
  Jinv = det(G)

  forward[I] = Metric(SMatrix{1,1,T,1}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{1,1,T,1}(Tuple(G)), T(J))
end

@kernel function _fill_cell_metric_storage_2d_kernel!(
  forward,
  inverse,
  jacobian,
  inverse_jacobian,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half

  F = _as_smatrix(Val(2), jacobian(t, ξ, η, params))
  G = _as_smatrix(Val(2), inverse_jacobian(t, ξ, η, params))
  J = det(F)
  Jinv = det(G)

  forward[I] = Metric(SMatrix{2,2,T,4}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{2,2,T,4}(Tuple(G)), T(J))
end

@kernel function _fill_cell_metric_storage_3d_kernel!(
  forward,
  inverse,
  jacobian,
  inverse_jacobian,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half
  ζ = Iglobal.I[3] - nhalo + half

  F = _as_smatrix(Val(3), jacobian(t, ξ, η, ζ, params))
  G = _as_smatrix(Val(3), inverse_jacobian(t, ξ, η, ζ, params))
  J = det(F)
  Jinv = det(G)

  forward[I] = Metric(SMatrix{3,3,T,9}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{3,3,T,9}(Tuple(G)), T(J))
end

@kernel function _fill_face_metric_storage_1d_axis1_kernel!(
  forward,
  inverse,
  conserved,
  edge,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half

  G, Ghat, Jinv = _inverse_and_normalized_edge_metrics(Val(1), edge, 1, t, (ξ,), params)
  F = inv(G)
  J = det(F)

  forward[I] = Metric(SMatrix{1,1,T,1}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{1,1,T,1}(Tuple(G)), T(J))
  conserved[I] = ConservedMetric(SMatrix{1,1,T,1}(Tuple(Ghat)))
end

@kernel function _fill_face_metric_storage_2d_axis1_kernel!(
  forward,
  inverse,
  conserved,
  edge,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half

  G, Ghat, Jinv = _inverse_and_normalized_edge_metrics(Val(2), edge, 1, t, (ξ, η), params)
  F = inv(G)
  J = det(F)

  forward[I] = Metric(SMatrix{2,2,T,4}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{2,2,T,4}(Tuple(G)), T(J))
  conserved[I] = ConservedMetric(SMatrix{2,2,T,4}(Tuple(Ghat)))
end

@kernel function _fill_face_metric_storage_2d_axis2_kernel!(
  forward,
  inverse,
  conserved,
  edge,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half

  G, Ghat, Jinv = _inverse_and_normalized_edge_metrics(Val(2), edge, 2, t, (ξ, η), params)
  F = inv(G)
  J = det(F)

  forward[I] = Metric(SMatrix{2,2,T,4}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{2,2,T,4}(Tuple(G)), T(J))
  conserved[I] = ConservedMetric(SMatrix{2,2,T,4}(Tuple(Ghat)))
end

@kernel function _fill_face_metric_storage_3d_axis1_kernel!(
  forward,
  inverse,
  conserved,
  edge,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half
  ζ = Iglobal.I[3] - nhalo + half

  G, Ghat, Jinv = _inverse_and_normalized_edge_metrics(
    Val(3), edge, 1, t, (ξ, η, ζ), params
  )
  F = inv(G)
  J = det(F)

  forward[I] = Metric(SMatrix{3,3,T,9}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{3,3,T,9}(Tuple(G)), T(J))
  conserved[I] = ConservedMetric(SMatrix{3,3,T,9}(Tuple(Ghat)))
end

@kernel function _fill_face_metric_storage_3d_axis2_kernel!(
  forward,
  inverse,
  conserved,
  edge,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half
  ζ = Iglobal.I[3] - nhalo + half

  G, Ghat, Jinv = _inverse_and_normalized_edge_metrics(
    Val(3), edge, 2, t, (ξ, η, ζ), params
  )
  F = inv(G)
  J = det(F)

  forward[I] = Metric(SMatrix{3,3,T,9}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{3,3,T,9}(Tuple(G)), T(J))
  conserved[I] = ConservedMetric(SMatrix{3,3,T,9}(Tuple(Ghat)))
end

@kernel function _fill_face_metric_storage_3d_axis3_kernel!(
  forward,
  inverse,
  conserved,
  edge,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half
  ζ = Iglobal.I[3] - nhalo + half

  G, Ghat, Jinv = _inverse_and_normalized_edge_metrics(
    Val(3), edge, 3, t, (ξ, η, ζ), params
  )
  F = inv(G)
  J = det(F)

  forward[I] = Metric(SMatrix{3,3,T,9}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{3,3,T,9}(Tuple(G)), T(J))
  conserved[I] = ConservedMetric(SMatrix{3,3,T,9}(Tuple(Ghat)))
end

function _fill_cell_metric_storage!(
  cell_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
  backend,
) where {N,T}
  _fill_cell_metric_storage!(
    backend,
    cell_metric_storage,
    metric_functions_cache,
    iterators,
    nhalo,
    t,
    params,
    Val(N),
    T,
  )
end

function _fill_cell_metric_storage!(
  backend::KernelAbstractions.Backend,
  cell_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
) where {N,T}
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full

  if N == 1
    _fill_cell_metric_storage_1d_kernel!(backend)(
      cell_metric_storage.forward,
      cell_metric_storage.inverse,
      metric_functions_cache.forward.jacobian,
      metric_functions_cache.inverse.Jinv,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
  elseif N == 2
    _fill_cell_metric_storage_2d_kernel!(backend)(
      cell_metric_storage.forward,
      cell_metric_storage.inverse,
      metric_functions_cache.forward.jacobian,
      metric_functions_cache.inverse.Jinv,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
  elseif N == 3
    _fill_cell_metric_storage_3d_kernel!(backend)(
      cell_metric_storage.forward,
      cell_metric_storage.inverse,
      metric_functions_cache.forward.jacobian,
      metric_functions_cache.inverse.Jinv,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
  else
    throw(ArgumentError("Unsupported metric dimension N=$N"))
  end
  KernelAbstractions.synchronize(backend)
  return nothing
end

function _fill_face_metric_storage!(
  face_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
  backend,
) where {N,T}
  _fill_face_metric_storage!(
    backend,
    face_metric_storage,
    metric_functions_cache,
    iterators,
    nhalo,
    t,
    params,
    Val(N),
    T,
  )
  return nothing
end

function _fill_face_metric_storage!(
  backend::KernelAbstractions.CPU,
  face_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
) where {N,T}
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  edge = metric_functions_cache.edge

  if N == 3 && hasproperty(edge, :ad_thomas_lombard_metric)
    _fill_face_metric_storage_ad_thomas_lombard_3d_cpu!(
      face_metric_storage,
      metric_functions_cache,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T,
    )
    return nothing
  end

  return _fill_face_metric_storage_kernel!(
    backend,
    face_metric_storage,
    metric_functions_cache,
    iterators,
    nhalo,
    t,
    params,
    Val(N),
    T,
  )
end

function _fill_face_metric_storage!(
  backend::KernelAbstractions.Backend,
  face_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
) where {N,T}
  return _fill_face_metric_storage_kernel!(
    backend,
    face_metric_storage,
    metric_functions_cache,
    iterators,
    nhalo,
    t,
    params,
    Val(N),
    T,
  )
end

function _fill_face_metric_storage_kernel!(
  backend::KernelAbstractions.Backend,
  face_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
) where {N,T}
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  edge = metric_functions_cache.edge

  if N == 1
    _fill_face_metric_storage_1d_axis1_kernel!(backend)(
      face_metric_storage[1].forward,
      face_metric_storage[1].inverse,
      face_metric_storage[1].conserved,
      edge,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
  elseif N == 2
    _fill_face_metric_storage_2d_axis1_kernel!(backend)(
      face_metric_storage[1].forward,
      face_metric_storage[1].inverse,
      face_metric_storage[1].conserved,
      edge,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
    _fill_face_metric_storage_2d_axis2_kernel!(backend)(
      face_metric_storage[2].forward,
      face_metric_storage[2].inverse,
      face_metric_storage[2].conserved,
      edge,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
  elseif N == 3
    _fill_face_metric_storage_3d_axis1_kernel!(backend)(
      face_metric_storage[1].forward,
      face_metric_storage[1].inverse,
      face_metric_storage[1].conserved,
      edge,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
    _fill_face_metric_storage_3d_axis2_kernel!(backend)(
      face_metric_storage[2].forward,
      face_metric_storage[2].inverse,
      face_metric_storage[2].conserved,
      edge,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
    _fill_face_metric_storage_3d_axis3_kernel!(backend)(
      face_metric_storage[3].forward,
      face_metric_storage[3].inverse,
      face_metric_storage[3].conserved,
      edge,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
  else
    throw(ArgumentError("Unsupported face metric dimension N=$N"))
  end
  KernelAbstractions.synchronize(backend)
  return nothing
end

@inline function _padded_global_cell_index(
  global_domain, i::Int, j::Int, k::Int, ni::Int, nj::Int, nk::Int
)
  ic = clamp(i, 1, ni)
  jc = clamp(j, 1, nj)
  kc = clamp(k, 1, nk)
  Iglobal = global_domain[CartesianIndex(ic, jc, kc)]
  return (
    Iglobal.I[1] + i - ic,
    Iglobal.I[2] + j - jc,
    Iglobal.I[3] + k - kc,
  )
end

@inline function _shift_array_axis(
  i::Int, j::Int, k::Int, axis::Int, offset::Int
)
  if axis == 1
    return i + offset, j, k
  elseif axis == 2
    return i, j + offset, k
  elseif axis == 3
    return i, j, k + offset
  else
    throw(ArgumentError("Invalid 3D derivative axis: $axis"))
  end
end

@inline _shift_array_axis(i::Int, j::Int, k::Int, ::Val{1}, offset::Int) =
  (i + offset, j, k)
@inline _shift_array_axis(i::Int, j::Int, k::Int, ::Val{2}, offset::Int) =
  (i, j + offset, k)
@inline _shift_array_axis(i::Int, j::Int, k::Int, ::Val{3}, offset::Int) =
  (i, j, k + offset)

@inline function _ad_thomas_lombard_potential(edge, potential_id::Int, t, ξ, η, ζ, params)
  F = edge.jacobian(t, ξ, η, ζ, params)

  return _ad_thomas_lombard_potential_from_values(
    F,
    edge.x(t, ξ, η, ζ, params),
    edge.y(t, ξ, η, ζ, params),
    edge.z(t, ξ, η, ζ, params),
    potential_id,
  )
end

@inline function _ad_thomas_lombard_potential_from_values(F, x, y, z, potential_id::Int)
  if potential_id == 1
    return F[2, 2] * z
  elseif potential_id == 2
    return F[2, 3] * z
  elseif potential_id == 3
    return F[2, 1] * z
  elseif potential_id == 4
    return F[3, 2] * x
  elseif potential_id == 5
    return F[3, 3] * x
  elseif potential_id == 6
    return F[3, 1] * x
  elseif potential_id == 7
    return F[1, 2] * y
  elseif potential_id == 8
    return F[1, 3] * y
  elseif potential_id == 9
    return F[1, 1] * y
  else
    throw(ArgumentError("Invalid AD Thomas-Lombard potential id: $potential_id"))
  end
end

@inline function _ad_thomas_lombard_potential_vector(
  edge, potential_ids::NTuple{3,Int}, t, ξ, η, ζ, params
)
  F = edge.jacobian(t, ξ, η, ζ, params)
  x = edge.x(t, ξ, η, ζ, params)
  y = edge.y(t, ξ, η, ζ, params)
  z = edge.z(t, ξ, η, ζ, params)
  return SVector{3}(
    ntuple(
      n -> _ad_thomas_lombard_potential_from_values(F, x, y, z, potential_ids[n]),
      Val(3),
    ),
  )
end

@inline function _ad_thomas_lombard_potential_vector(
  edge, potential_ids::Val, t, ξ, η, ζ, params
)
  return _ad_thomas_lombard_potential_vector_from_axis_derivative(
    edge, potential_ids, t, ξ, η, ζ, params
  )
end

@generated function _ad_thomas_lombard_seed_axis(::Val{Axis}, ξ, η, ζ) where {Axis}
  Axis in (1, 2, 3) || error("Invalid AD Thomas-Lombard seed axis: $Axis")
  vars = (:ξ, :η, :ζ)
  coord = vars[Axis]
  seeded = collect(vars)
  seeded[Axis] = :seeded_coord
  return quote
    coord_value = $coord + zero(ξ) + zero(η) + zero(ζ)
    tag = ForwardDiff.Tag(_ad_thomas_lombard_seed_axis, typeof(coord_value))
    tag_type = typeof(tag)
    seeded_coord = ForwardDiff.Dual{tag_type}(coord_value, oneunit(coord_value))
    return ($(seeded[1]), $(seeded[2]), $(seeded[3]), tag_type)
  end
end

@inline _ad_thomas_lombard_pair_axis_a_tag() = nothing
@inline _ad_thomas_lombard_pair_axis_b_tag() = nothing

@inline function _ad_thomas_lombard_seed_second_axis(
  ::Val{1}, tag_function, ξ, η, ζ
)
  coord_value = ξ
  tag1 = ForwardDiff.Tag(tag_function, typeof(coord_value))
  tag1_type = typeof(tag1)
  coord1 = ForwardDiff.Dual{tag1_type}(coord_value, oneunit(coord_value))
  tag2 = ForwardDiff.Tag(tag_function, typeof(coord1))
  tag2_type = typeof(tag2)
  coord2 = ForwardDiff.Dual{tag2_type}(coord1, oneunit(coord1))
  return coord2, η, ζ, tag1_type, tag2_type
end

@inline function _ad_thomas_lombard_seed_second_axis(
  ::Val{2}, tag_function, ξ, η, ζ
)
  coord_value = η
  tag1 = ForwardDiff.Tag(tag_function, typeof(coord_value))
  tag1_type = typeof(tag1)
  coord1 = ForwardDiff.Dual{tag1_type}(coord_value, oneunit(coord_value))
  tag2 = ForwardDiff.Tag(tag_function, typeof(coord1))
  tag2_type = typeof(tag2)
  coord2 = ForwardDiff.Dual{tag2_type}(coord1, oneunit(coord1))
  return ξ, coord2, ζ, tag1_type, tag2_type
end

@inline function _ad_thomas_lombard_seed_second_axis(
  ::Val{3}, tag_function, ξ, η, ζ
)
  coord_value = ζ
  tag1 = ForwardDiff.Tag(tag_function, typeof(coord_value))
  tag1_type = typeof(tag1)
  coord1 = ForwardDiff.Dual{tag1_type}(coord_value, oneunit(coord_value))
  tag2 = ForwardDiff.Tag(tag_function, typeof(coord1))
  tag2_type = typeof(tag2)
  coord2 = ForwardDiff.Dual{tag2_type}(coord1, oneunit(coord1))
  return ξ, η, coord2, tag1_type, tag2_type
end

@inline function _ad_thomas_lombard_coordinate_axis_values(
  edge, axis::Val, t, ξ, η, ζ, params
)
  ξd, ηd, ζd, tag_type = _ad_thomas_lombard_seed_axis(axis, ξ, η, ζ)
  xd = edge.x(t, ξd, ηd, ζd, params)
  yd = edge.y(t, ξd, ηd, ζd, params)
  zd = edge.z(t, ξd, ηd, ζd, params)
  x = _dual_value_for_tag(tag_type, xd)
  y = _dual_value_for_tag(tag_type, yd)
  z = _dual_value_for_tag(tag_type, zd)
  dx = _dual_derivative_for_tag(tag_type, xd)
  dy = _dual_derivative_for_tag(tag_type, yd)
  dz = _dual_derivative_for_tag(tag_type, zd)
  return x, y, z, dx, dy, dz
end

@inline function _ad_thomas_lombard_potential_vector_from_axis_derivative(
  edge, ::Val{(1, 4, 7)}, t, ξ, η, ζ, params
)
  x, y, z, dx, dy, dz = _ad_thomas_lombard_coordinate_axis_values(
    edge, Val(2), t, ξ, η, ζ, params
  )
  return SVector(dy * z, dz * x, dx * y)
end

@inline function _ad_thomas_lombard_potential_vector_from_axis_derivative(
  edge, ::Val{(2, 5, 8)}, t, ξ, η, ζ, params
)
  x, y, z, dx, dy, dz = _ad_thomas_lombard_coordinate_axis_values(
    edge, Val(3), t, ξ, η, ζ, params
  )
  return SVector(dy * z, dz * x, dx * y)
end

@inline function _ad_thomas_lombard_potential_vector_from_axis_derivative(
  edge, ::Val{(3, 6, 9)}, t, ξ, η, ζ, params
)
  x, y, z, dx, dy, dz = _ad_thomas_lombard_coordinate_axis_values(
    edge, Val(1), t, ξ, η, ζ, params
  )
  return SVector(dy * z, dz * x, dx * y)
end

@inline _dual_value_for_tag(::Type{Tag}, y) where {Tag} = ForwardDiff.value(Tag, y)
@inline _dual_derivative_for_tag(::Type{Tag}, y) where {Tag} =
  ForwardDiff.partials(Tag, y, 1)

@inline function _dual_value_for_tag(::Type{Tag}, y::SVector{N}) where {Tag,N}
  return SVector{N}(ntuple(i -> _dual_value_for_tag(Tag, y[i]), Val(N)))
end

@inline function _dual_derivative_for_tag(::Type{Tag}, y::SVector{N}) where {Tag,N}
  return SVector{N}(ntuple(i -> _dual_derivative_for_tag(Tag, y[i]), Val(N)))
end

@inline function _dual_value_for_tag(::Type{Tag}, y::Tuple) where {Tag}
  return map(v -> _dual_value_for_tag(Tag, v), y)
end

@inline function _dual_derivative_for_tag(::Type{Tag}, y::Tuple) where {Tag}
  return map(v -> _dual_derivative_for_tag(Tag, v), y)
end

@inline function _forwarddiff_extract_value_derivative_and_second_derivative(
  ::Type{T1}, ::Type{T2}, ydual2
) where {T1,T2}
  ydual1 = _dual_value_for_tag(T2, ydual2)
  value = _dual_value_for_tag(T1, ydual1)
  deriv = _dual_derivative_for_tag(T1, ydual1)
  second_deriv = _dual_derivative_for_tag(T1, _dual_derivative_for_tag(T2, ydual2))
  return value, deriv, second_deriv
end

@inline function _forwarddiff_value_derivative_and_second_derivative(f, x)
  tag1 = ForwardDiff.Tag(f, typeof(x))
  T1 = typeof(tag1)
  xdual1 = ForwardDiff.Dual{T1}(x, oneunit(x))
  tag2 = ForwardDiff.Tag(f, typeof(xdual1))
  T2 = typeof(tag2)
  xdual2 = ForwardDiff.Dual{T2}(xdual1, oneunit(xdual1))
  ydual2 = f(xdual2)
  return _forwarddiff_extract_value_derivative_and_second_derivative(T1, T2, ydual2)
end

function _ad_thomas_lombard_potential_axis_derivs(
  edge, potential_id::Int, axis::Int, t, ξ, η, ζ, params
)
  backend = edge.diff_backend
  if axis == 1
    return value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_potential(edge, potential_id, t, s, η, ζ, params),
      backend,
      ξ,
    )
  elseif axis == 2
    return value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_potential(edge, potential_id, t, ξ, s, ζ, params),
      backend,
      η,
    )
  elseif axis == 3
    return value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_potential(edge, potential_id, t, ξ, η, s, params),
      backend,
      ζ,
    )
  else
    throw(ArgumentError("Invalid AD Thomas-Lombard potential derivative axis: $axis"))
  end
end

function _ad_thomas_lombard_potential_axis_derivs_vector(
  edge, potential_ids, axis::Int, t, ξ, η, ζ, params
)
  return _ad_thomas_lombard_potential_axis_derivs_vector(
    edge, potential_ids, Val(axis), t, ξ, η, ζ, params
  )
end

function _ad_thomas_lombard_potential_axis_derivs_vector(
  edge, potential_ids, ::Val{Axis}, t, ξ, η, ζ, params
) where {Axis}
  if Axis == 1
    return _forwarddiff_value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_potential_vector(edge, potential_ids, t, s, η, ζ, params),
      ξ,
    )
  elseif Axis == 2
    return _forwarddiff_value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_potential_vector(edge, potential_ids, t, ξ, s, ζ, params),
      η,
    )
  elseif Axis == 3
    return _forwarddiff_value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_potential_vector(edge, potential_ids, t, ξ, η, s, params),
      ζ,
    )
  else
    throw(ArgumentError("Invalid AD Thomas-Lombard potential derivative axis: $Axis"))
  end
end

function _ad_thomas_lombard_edge_potential(
  edge, potential_id::Int, edge_axis::Int, t, ξ, η, ζ, params
)
  ϕᵢ, ϕaᵢ, ϕaaᵢ = _ad_thomas_lombard_potential_axis_derivs(
    edge, potential_id, edge_axis, t, ξ, η, ζ, params
  )

  if edge_axis == 1
    ϕᵢ₊₁, ϕaᵢ₊₁, ϕaaᵢ₊₁ = _ad_thomas_lombard_potential_axis_derivs(
      edge, potential_id, edge_axis, t, ξ + one(ξ), η, ζ, params
    )
  elseif edge_axis == 2
    ϕᵢ₊₁, ϕaᵢ₊₁, ϕaaᵢ₊₁ = _ad_thomas_lombard_potential_axis_derivs(
      edge, potential_id, edge_axis, t, ξ, η + one(η), ζ, params
    )
  elseif edge_axis == 3
    ϕᵢ₊₁, ϕaᵢ₊₁, ϕaaᵢ₊₁ = _ad_thomas_lombard_potential_axis_derivs(
      edge, potential_id, edge_axis, t, ξ, η, ζ + one(ζ), params
    )
  else
    throw(ArgumentError("Invalid AD Thomas-Lombard edge axis: $edge_axis"))
  end

  return _edge_reconstruct(
    ϕᵢ, ϕaᵢ, ϕaaᵢ, ϕᵢ₊₁, ϕaᵢ₊₁, ϕaaᵢ₊₁, CurvatureCorrectedReconstruction()
  )
end

function _ad_thomas_lombard_edge_potential_derivs(
  edge, potential_id::Int, edge_axis::Int, deriv_axis::Int, t, ξ, η, ζ, params
)
  backend = edge.diff_backend
  if deriv_axis == 1
    return value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_edge_potential(
        edge, potential_id, edge_axis, t, s, η, ζ, params
      ),
      backend,
      ξ,
    )
  elseif deriv_axis == 2
    return value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_edge_potential(
        edge, potential_id, edge_axis, t, ξ, s, ζ, params
      ),
      backend,
      η,
    )
  elseif deriv_axis == 3
    return value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_edge_potential(
        edge, potential_id, edge_axis, t, ξ, η, s, params
      ),
      backend,
      ζ,
    )
  else
    throw(ArgumentError("Invalid AD Thomas-Lombard derivative axis: $deriv_axis"))
  end
end

function _ad_thomas_lombard_potential_pair_derivs_vector(
  edge, potential_ids, axis_a::Int, axis_b::Int, t, ξ, η, ζ, params
)
  return _ad_thomas_lombard_potential_pair_derivs_vector(
    edge, potential_ids, Val(axis_a), Val(axis_b), t, ξ, η, ζ, params
  )
end

function _ad_thomas_lombard_potential_pair_derivs_vector(
  edge, potential_ids, ::Val{AxisA}, ::Val{AxisB}, t, ξ, η, ζ, params
) where {AxisA,AxisB}
  ξb, ηb, ζb, Tb1, Tb2 = _ad_thomas_lombard_seed_second_axis(
    Val(AxisB), _ad_thomas_lombard_pair_axis_b_tag, ξ, η, ζ
  )
  ξab, ηab, ζab, Ta1, Ta2 = _ad_thomas_lombard_seed_second_axis(
    Val(AxisA), _ad_thomas_lombard_pair_axis_a_tag, ξb, ηb, ζb
  )

  Pdual = _ad_thomas_lombard_potential_vector(
    edge, potential_ids, t, ξab, ηab, ζab, params
  )
  P_bdual, Pa_bdual, Paa_bdual =
    _forwarddiff_extract_value_derivative_and_second_derivative(Ta1, Ta2, Pdual)

  P, Pb, Pbb =
    _forwarddiff_extract_value_derivative_and_second_derivative(Tb1, Tb2, P_bdual)
  Pa, Pab, Pabb =
    _forwarddiff_extract_value_derivative_and_second_derivative(Tb1, Tb2, Pa_bdual)
  Paa, Paab, Paabb =
    _forwarddiff_extract_value_derivative_and_second_derivative(Tb1, Tb2, Paa_bdual)
  return P, Pa, Paa, Pb, Pab, Paab, Pbb, Pabb, Paabb
end

function _fill_ad_thomas_lombard_edge_potential_derivs_3d!(
  Q,
  Qb,
  Qbb,
  edge,
  potential_id::Int,
  edge_axis::Int,
  deriv_axis::Int,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  ni, nj, nk = size(local_domain)
  half = T(0.5)
  irange = deriv_axis == 1 ? (1:(ni + 2)) : (2:(ni + 1))
  jrange = deriv_axis == 2 ? (1:(nj + 2)) : (2:(nj + 1))
  krange = deriv_axis == 3 ? (1:(nk + 2)) : (2:(nk + 1))

  Threads.@threads for kp in krange
    @inbounds for jp in jrange, ip in irange
      i = ip - 1
      j = jp - 1
      k = kp - 1
      Iglobal = _padded_global_cell_index(global_domain, i, j, k, ni, nj, nk)
      ξ = T(Iglobal[1] - nhalo) + half
      η = T(Iglobal[2] - nhalo) + half
      ζ = T(Iglobal[3] - nhalo) + half

      q, qb, qbb = _ad_thomas_lombard_edge_potential_derivs(
        edge, potential_id, edge_axis, deriv_axis, t, ξ, η, ζ, params
      )
      Q[ip, jp, kp] = T(q)
      Qb[ip, jp, kp] = T(qb)
      Qbb[ip, jp, kp] = T(qbb)
    end
  end

  return nothing
end

function _fill_ad_thomas_lombard_potential_pair_table_3d!(
  table,
  edge,
  potential_ids,
  axis_a::Int,
  axis_b::Int,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  return _fill_ad_thomas_lombard_potential_pair_table_3d!(
    table,
    edge,
    potential_ids,
    Val(axis_a),
    Val(axis_b),
    local_domain,
    global_domain,
    nhalo,
    t,
    params,
    T,
  )
end

function _fill_ad_thomas_lombard_potential_pair_table_3d!(
  table,
  edge,
  potential_ids,
  ::Val{AxisA},
  ::Val{AxisB},
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T,AxisA,AxisB}
  ni, nj, nk = size(local_domain)
  half = T(0.5)
  irange = AxisA == 1 || AxisB == 1 ? (1:(ni + 2)) : (2:(ni + 1))
  jrange = AxisA == 2 || AxisB == 2 ? (1:(nj + 2)) : (2:(nj + 1))
  krange = AxisA == 3 || AxisB == 3 ? (1:(nk + 2)) : (2:(nk + 1))

  Threads.@threads for kp in krange
    @inbounds for jp in jrange, ip in irange
      i = ip - 1
      j = jp - 1
      k = kp - 1
      Iglobal = _padded_global_cell_index(global_domain, i, j, k, ni, nj, nk)
      ξ = T(Iglobal[1] - nhalo) + half
      η = T(Iglobal[2] - nhalo) + half
      ζ = T(Iglobal[3] - nhalo) + half

      vals = _ad_thomas_lombard_potential_pair_derivs_vector(
        edge, potential_ids, Val(AxisA), Val(AxisB), t, ξ, η, ζ, params
      )
      table[1][ip, jp, kp] = SVector{3,T}(vals[1])
      table[2][ip, jp, kp] = SVector{3,T}(vals[2])
      table[3][ip, jp, kp] = SVector{3,T}(vals[3])
      table[4][ip, jp, kp] = SVector{3,T}(vals[4])
      table[5][ip, jp, kp] = SVector{3,T}(vals[5])
      table[6][ip, jp, kp] = SVector{3,T}(vals[6])
      table[7][ip, jp, kp] = SVector{3,T}(vals[7])
      table[8][ip, jp, kp] = SVector{3,T}(vals[8])
      table[9][ip, jp, kp] = SVector{3,T}(vals[9])
    end
  end

  return nothing
end

@inline function _array_derivative_from_ad_thomas_lombard_edge_potential(
  Q, Qb, Qbb, deriv_axis::Int, i::Int, j::Int, k::Int
)
  ip, jp, kp = _shift_array_axis(i, j, k, deriv_axis, 1)
  im, jm, km = _shift_array_axis(i, j, k, deriv_axis, -1)

  edge_high = _edge_reconstruct(
    Q[i, j, k],
    Qb[i, j, k],
    Qbb[i, j, k],
    Q[ip, jp, kp],
    Qb[ip, jp, kp],
    Qbb[ip, jp, kp],
    CurvatureCorrectedReconstruction(),
  )
  edge_low = _edge_reconstruct(
    Q[im, jm, km],
    Qb[im, jm, km],
    Qbb[im, jm, km],
    Q[i, j, k],
    Qb[i, j, k],
    Qbb[i, j, k],
    CurvatureCorrectedReconstruction(),
  )

  return edge_high - edge_low
end

@inline function _edge_reconstruct_from_pair_table(
  table,
  value_id::Int,
  normal_deriv_id::Int,
  normal_second_deriv_id::Int,
  normal_axis::Int,
  i::Int,
  j::Int,
  k::Int,
)
  return _edge_reconstruct_from_pair_table(
    table,
    value_id,
    normal_deriv_id,
    normal_second_deriv_id,
    Val(normal_axis),
    i,
    j,
    k,
  )
end

@inline function _edge_reconstruct_from_pair_table(
  table,
  value_id::Int,
  normal_deriv_id::Int,
  normal_second_deriv_id::Int,
  normal_axis::Val,
  i::Int,
  j::Int,
  k::Int,
)
  ip, jp, kp = _shift_array_axis(i, j, k, normal_axis, 1)
  return _edge_reconstruct(
    table[value_id][i, j, k],
    table[normal_deriv_id][i, j, k],
    table[normal_second_deriv_id][i, j, k],
    table[value_id][ip, jp, kp],
    table[normal_deriv_id][ip, jp, kp],
    table[normal_second_deriv_id][ip, jp, kp],
    CurvatureCorrectedReconstruction(),
  )
end

@inline function _edge_potential_values_from_pair_table(
  table,
  normal_axis::Int,
  normal_is_axis_a::Bool,
  i::Int,
  j::Int,
  k::Int,
)
  return _edge_potential_values_from_pair_table(
    table, Val(normal_axis), Val(normal_is_axis_a), i, j, k
  )
end

@inline function _edge_potential_values_from_pair_table(
  table,
  normal_axis::Val,
  ::Val{NormalIsAxisA},
  i::Int,
  j::Int,
  k::Int,
) where {NormalIsAxisA}
  if NormalIsAxisA
    Q = _edge_reconstruct_from_pair_table(table, 1, 2, 3, normal_axis, i, j, k)
    Qd = _edge_reconstruct_from_pair_table(table, 4, 5, 6, normal_axis, i, j, k)
    Qdd = _edge_reconstruct_from_pair_table(table, 7, 8, 9, normal_axis, i, j, k)
  else
    Q = _edge_reconstruct_from_pair_table(table, 1, 4, 7, normal_axis, i, j, k)
    Qd = _edge_reconstruct_from_pair_table(table, 2, 5, 8, normal_axis, i, j, k)
    Qdd = _edge_reconstruct_from_pair_table(table, 3, 6, 9, normal_axis, i, j, k)
  end
  return Q, Qd, Qdd
end

@inline function _array_derivative_from_pair_table(
  table,
  normal_axis::Int,
  deriv_axis::Int,
  normal_is_axis_a::Bool,
  i::Int,
  j::Int,
  k::Int,
)
  return _array_derivative_from_pair_table(
    table,
    Val(normal_axis),
    Val(deriv_axis),
    Val(normal_is_axis_a),
    i,
    j,
    k,
  )
end

@inline function _array_derivative_from_pair_table(
  table,
  normal_axis::Val,
  deriv_axis::Val,
  normal_is_axis_a::Val,
  i::Int,
  j::Int,
  k::Int,
)
  ip, jp, kp = _shift_array_axis(i, j, k, deriv_axis, 1)
  im, jm, km = _shift_array_axis(i, j, k, deriv_axis, -1)

  Q, Qd, Qdd = _edge_potential_values_from_pair_table(
    table, normal_axis, normal_is_axis_a, i, j, k
  )
  Qp, Qdp, Qddp = _edge_potential_values_from_pair_table(
    table, normal_axis, normal_is_axis_a, ip, jp, kp
  )
  Qm, Qdm, Qddm = _edge_potential_values_from_pair_table(
    table, normal_axis, normal_is_axis_a, im, jm, km
  )

  edge_high = _edge_reconstruct(Q, Qd, Qdd, Qp, Qdp, Qddp, CurvatureCorrectedReconstruction())
  edge_low = _edge_reconstruct(Qm, Qdm, Qddm, Q, Qd, Qdd, CurvatureCorrectedReconstruction())
  return edge_high - edge_low
end

function _accumulate_ad_thomas_lombard_term_3d!(
  out,
  Q,
  Qb,
  Qbb,
  sign,
  edge,
  potential_id::Int,
  edge_axis::Int,
  deriv_axis::Int,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  _fill_ad_thomas_lombard_edge_potential_derivs_3d!(
    Q,
    Qb,
    Qbb,
    edge,
    potential_id,
    edge_axis,
    deriv_axis,
    local_domain,
    global_domain,
    nhalo,
    t,
    params,
    T,
  )

  Threads.@threads for I in local_domain
    @inbounds begin
      i, j, k = I.I
      out[I] += sign * _array_derivative_from_ad_thomas_lombard_edge_potential(
        Q, Qb, Qbb, deriv_axis, i + 1, j + 1, k + 1
      )
    end
  end

  return nothing
end

function _accumulate_ad_thomas_lombard_pair_orientations_3d!(
  out_a,
  sign_a,
  out_b,
  sign_b,
  table,
  axis_a::Val,
  axis_b::Val,
  local_domain,
)
  Threads.@threads for I in local_domain
    @inbounds begin
      i, j, k = I.I
      ip = i + 1
      jp = j + 1
      kp = k + 1
      dQa = _array_derivative_from_pair_table(
        table, axis_a, axis_b, Val(true), ip, jp, kp
      )
      dQb = _array_derivative_from_pair_table(
        table, axis_b, axis_a, Val(false), ip, jp, kp
      )
      out_a[I] += sign_a * dQa
      out_b[I] += sign_b * dQb
    end
  end

  return nothing
end

function _accumulate_ad_thomas_lombard_pair_3d!(
  out_a,
  sign_a,
  out_b,
  sign_b,
  table,
  edge,
  potential_ids,
  axis_a::Int,
  axis_b::Int,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  return _accumulate_ad_thomas_lombard_pair_3d!(
    out_a,
    sign_a,
    out_b,
    sign_b,
    table,
    edge,
    potential_ids,
    Val(axis_a),
    Val(axis_b),
    local_domain,
    global_domain,
    nhalo,
    t,
    params,
    T,
  )
end

function _accumulate_ad_thomas_lombard_pair_3d!(
  out_a,
  sign_a,
  out_b,
  sign_b,
  table,
  edge,
  potential_ids,
  axis_a::Val,
  axis_b::Val,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  _fill_ad_thomas_lombard_potential_pair_table_3d!(
    table,
    edge,
    potential_ids,
    axis_a,
    axis_b,
    local_domain,
    global_domain,
    nhalo,
    t,
    params,
    T,
  )
  _accumulate_ad_thomas_lombard_pair_orientations_3d!(
    out_a,
    sign_a,
    out_b,
    sign_b,
    table,
    axis_a,
    axis_b,
    local_domain,
  )

  return nothing
end

function _fill_face_metric_storage_ad_thomas_lombard_3d_cpu!(
  face_metric_storage,
  metric_functions_cache,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  edge = metric_functions_cache.edge
  zeroT = zero(T)
  half = T(0.5)
  row_size = size(local_domain)

  ni, nj, nk = size(local_domain)
  pad_size = (ni + 2, nj + 2, nk + 2)
  if edge.diff_backend isa AutoForwardDiff
    ξrow = fill(zero(SVector{3,T}), row_size)
    ηrow = fill(zero(SVector{3,T}), row_size)
    ζrow = fill(zero(SVector{3,T}), row_size)
    table = ntuple(_ -> Array{SVector{3,T}}(undef, pad_size), Val(9))
    _accumulate_ad_thomas_lombard_pair_3d!(
      ξrow,
      one(T),
      ζrow,
      -one(T),
      table,
      edge,
      Val((1, 4, 7)),
      Val(1),
      Val(3),
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T,
    )
    _accumulate_ad_thomas_lombard_pair_3d!(
      ξrow,
      -one(T),
      ηrow,
      one(T),
      table,
      edge,
      Val((2, 5, 8)),
      Val(1),
      Val(2),
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T,
    )
    _accumulate_ad_thomas_lombard_pair_3d!(
      ηrow,
      -one(T),
      ζrow,
      one(T),
      table,
      edge,
      Val((3, 6, 9)),
      Val(2),
      Val(3),
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T,
    )

    Threads.@threads for I in local_domain
      @inbounds begin
        Iglobal = global_domain[I]
        ξc = T(Iglobal.I[1] - nhalo) + half
        ηc = T(Iglobal.I[2] - nhalo) + half
        ζc = T(Iglobal.I[3] - nhalo) + half

        for axis in 1:3
          ξf, ηf, ζf = _face_center_3d(axis, ξc, ηc, ζc)
          F = _as_smatrix(Val(3), edge.jacobian(t, ξf, ηf, ζf, params))
          G = inv(F)
          J = det(F)

          face_metric_storage[axis].forward[I] = Metric(SMatrix{3,3,T,9}(Tuple(F)), T(J))
          face_metric_storage[axis].inverse[I] = Metric(SMatrix{3,3,T,9}(Tuple(G)), T(J))
        end

        ξ = ξrow[I]
        η = ηrow[I]
        ζ = ζrow[I]
        face_metric_storage[1].conserved[I] = ConservedMetric(
          @SMatrix [ξ[1] ξ[2] ξ[3]; zeroT zeroT zeroT; zeroT zeroT zeroT]
        )
        face_metric_storage[2].conserved[I] = ConservedMetric(
          @SMatrix [zeroT zeroT zeroT; η[1] η[2] η[3]; zeroT zeroT zeroT]
        )
        face_metric_storage[3].conserved[I] = ConservedMetric(
          @SMatrix [zeroT zeroT zeroT; zeroT zeroT zeroT; ζ[1] ζ[2] ζ[3]]
        )
      end
    end
  else
    ξx = fill(zeroT, row_size)
    ξy = fill(zeroT, row_size)
    ξz = fill(zeroT, row_size)
    ηx = fill(zeroT, row_size)
    ηy = fill(zeroT, row_size)
    ηz = fill(zeroT, row_size)
    ζx = fill(zeroT, row_size)
    ζy = fill(zeroT, row_size)
    ζz = fill(zeroT, row_size)
    Q = Array{T}(undef, pad_size)
    Qb = Array{T}(undef, pad_size)
    Qbb = Array{T}(undef, pad_size)
    terms = (
      (ξx, one(T), 1, 1, 3),
      (ξx, -one(T), 2, 1, 2),
      (ηx, one(T), 2, 2, 1),
      (ηx, -one(T), 3, 2, 3),
      (ζx, one(T), 3, 3, 2),
      (ζx, -one(T), 1, 3, 1),
      (ξy, one(T), 4, 1, 3),
      (ξy, -one(T), 5, 1, 2),
      (ηy, one(T), 5, 2, 1),
      (ηy, -one(T), 6, 2, 3),
      (ζy, one(T), 6, 3, 2),
      (ζy, -one(T), 4, 3, 1),
      (ξz, one(T), 7, 1, 3),
      (ξz, -one(T), 8, 1, 2),
      (ηz, one(T), 8, 2, 1),
      (ηz, -one(T), 9, 2, 3),
      (ζz, one(T), 9, 3, 2),
      (ζz, -one(T), 7, 3, 1),
    )
    for (out, sign, potential_id, edge_axis, deriv_axis) in terms
      _accumulate_ad_thomas_lombard_term_3d!(
        out,
        Q,
        Qb,
        Qbb,
        sign,
        edge,
        potential_id,
        edge_axis,
        deriv_axis,
        local_domain,
        global_domain,
        nhalo,
        t,
        params,
        T,
      )
    end

    Threads.@threads for I in local_domain
      @inbounds begin
        Iglobal = global_domain[I]
        ξc = T(Iglobal.I[1] - nhalo) + half
        ηc = T(Iglobal.I[2] - nhalo) + half
        ζc = T(Iglobal.I[3] - nhalo) + half

        for axis in 1:3
          ξf, ηf, ζf = _face_center_3d(axis, ξc, ηc, ζc)
          F = _as_smatrix(Val(3), edge.jacobian(t, ξf, ηf, ζf, params))
          G = inv(F)
          J = det(F)

          face_metric_storage[axis].forward[I] = Metric(SMatrix{3,3,T,9}(Tuple(F)), T(J))
          face_metric_storage[axis].inverse[I] = Metric(SMatrix{3,3,T,9}(Tuple(G)), T(J))
        end

        face_metric_storage[1].conserved[I] = ConservedMetric(
          @SMatrix [ξx[I] ξy[I] ξz[I]; zeroT zeroT zeroT; zeroT zeroT zeroT]
        )
        face_metric_storage[2].conserved[I] = ConservedMetric(
          @SMatrix [zeroT zeroT zeroT; ηx[I] ηy[I] ηz[I]; zeroT zeroT zeroT]
        )
        face_metric_storage[3].conserved[I] = ConservedMetric(
          @SMatrix [zeroT zeroT zeroT; zeroT zeroT zeroT; ζx[I] ζy[I] ζz[I]]
        )
      end
    end
  end

  return nothing
end

@inline function _face_center_3d(face_axis::Int, ξc::T, ηc::T, ζc::T) where {T}
  half = T(0.5)
  if face_axis == 1
    return ξc + half, ηc, ζc
  elseif face_axis == 2
    return ξc, ηc + half, ζc
  elseif face_axis == 3
    return ξc, ηc, ζc + half
  else
    throw(ArgumentError("Invalid 3D face axis: $face_axis"))
  end
end
