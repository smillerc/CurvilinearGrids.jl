module MetricTypes

using StaticArrays

export Metric1D, Metric2D, Metric3D
export contravariant_velocity

struct Metric1D{T} <: FieldVector{1,T}
  x₁::T
  t::T
end

struct Metric2D{T} <: FieldVector{2,T}
  x₁::T
  x₂::T
  t::T
end

struct Metric3D{T} <: FieldVector{3,T}
  x₁::T
  x₂::T
  x₃::T
  t::T
end

# Default constructors for Metrics
Metric1D(T::DataType=Float64) = Metric1D(zero(T), zero(T))
Metric2D(T::DataType=Float64) = Metric2D(zero(T), zero(T), zero(T))
Metric3D(T::DataType=Float64) = Metric3D(zero(T), zero(T), zero(T), zero(T))

# the time member defaults to 0, promote all to common datatype
Metric1D(x::Real) = Metric1D(promote(x, 0)...)
Metric2D(x::Real, y::Real) = Metric2D(promote(x, y, 0)...)
Metric3D(x::Real, y::Real, z::Real) = Metric3D(promote(x, y, z, 0)...)

# promote all arguments to a common type
Metric1D(x, t) = Metric1D(promote(x, t)...)
Metric2D(x, y, t) = Metric2D(promote(x, y, t)...)
Metric3D(x, y, z, t) = Metric3D(promote(x, y, z, t)...)

end
