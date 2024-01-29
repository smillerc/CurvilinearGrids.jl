
# 1st derivative operator
@inline function forward_∂_2nd_order(ϕ::OffsetVector{T,SVector{3,T}}) where {T}
  _∂ϕ = -(3 / 2) * ϕ[0] + 2ϕ[+1] - (1 / 2) * ϕ[+2]
  return _∂ϕ
end

@inline function backward_∂_2nd_order(ϕ::OffsetVector{T,SVector{3,T}}) where {T}
  _∂ϕ = +(3 / 2) * ϕ[0] - 2ϕ[-1] + (1 / 2) * ϕ[-2]
  return _∂ϕ
end

@inline function forward_mixed_∂_4th_order(ϕ::OffsetVector{T,SVector{5,T}}) where {T}
  # This is different from the Chandravamsi2023 paper, but these give the correct
  # derivative
  _∂ϕ = (1 / 12) * (-25ϕ[0] + 48ϕ[+1] - 36ϕ[+2] + 16ϕ[+3] - 3ϕ[+4])

  # original version, but this gave odd results in testing...
  # _∂ϕ = (
  #   -(1 / 4) * ϕ[0] - (5 / 6) * ϕ[+1] + (3 / 2) * ϕ[+2] - (1 / 2) * ϕ[+3] + (1 / 12) * ϕ[+4]
  # )
  return _∂ϕ
end

@inline function backward_mixed_∂_4th_order(ϕ::OffsetVector{T,SVector{5,T}}) where {T}
  # This is different from the Chandravamsi2023 paper, but these give the correct
  # derivative
  _∂ϕ = (1 / 12) * (+25ϕ[0] - 48ϕ[-1] + 36ϕ[-2] - 16ϕ[-3] + 3ϕ[-4])

  # original version, but this gave odd results in testing...
  # _∂ϕ = (
  #   -(1 / 4) * ϕ[0] - (5 / 6) * ϕ[-1] + (3 / 2) * ϕ[-2] - (1 / 2) * ϕ[-3] + (1 / 12) * ϕ[-4]
  # )
  return _∂ϕ
end

@inline function central_∂_4th_order(ϕ::OffsetVector{T,SVector{5,T}}) where {T}
  _∂ϕ = ((2 / 3) * (ϕ[+1] - ϕ[-1]) - (1 / 12) * (ϕ[+2] - ϕ[-2]))
  return _∂ϕ
end

@inline function ∂(ϕ::OffsetVector{T,SVector{7,T}}) where {T}
  f1 = 3 / 4
  f2 = 3 / 20
  f3 = 1 / 60
  _∂ϕ = f1 * (ϕ[+1] - ϕ[-1]) - f2 * (ϕ[+2] - ϕ[-2]) + f3 * (ϕ[+3] - ϕ[-3])
  return _∂ϕ
end

# 2nd derivative operator
@inline function ∂²(
  ϕ::OffsetVector{T,SVector{3,T}}, ∂ϕ::OffsetVector{T,SVector{3,T}}
) where {T}
  _∂²ϕ = 2(ϕ[+1] - 2ϕ[0] + ϕ[-1]) - (∂ϕ[+1] - ∂ϕ[-1]) / 2
  return _∂²ϕ
end

@inline function forward_∂²(∂ϕ::OffsetVector{T,SVector{3,T}}) where {T}
  _∂²ϕ = -(3 / 2) * ∂ϕ[+0] + 2∂ϕ[+1] - (1 / 2) * ∂ϕ[+2]
  return _∂²ϕ
end

@inline function forward_mixed_∂²(∂ϕ::OffsetVector{T,SVector{5,T}}) where {T}
  _∂²ϕ =
    -(1 / 4) * ∂ϕ[0] - (5 / 6) * ∂ϕ[+1] + (3 / 2) * ∂ϕ[+2] - (1 / 2) * ∂ϕ[+3] +
    (1 / 12) * ∂ϕ[+4]
  return _∂²ϕ
end

@inline function backward_∂²(∂ϕ::OffsetVector{T,SVector{3,T}}) where {T}
  _∂²ϕ = +(3 / 2) * ∂ϕ[+0] - 2∂ϕ[-1] + (1 / 2) * ∂ϕ[-2]
  return _∂²ϕ
end

@inline function backward_mixed_∂²(∂ϕ::OffsetVector{T,SVector{5,T}}) where {T}
  _∂²ϕ =
    +(1 / 4) * ∂ϕ[0] + (5 / 6) * ∂ϕ[-1] - (3 / 2) * ∂ϕ[-2] + (1 / 2) * ∂ϕ[-3] -
    (1 / 12) * ∂ϕ[-4]
  return _∂²ϕ
end