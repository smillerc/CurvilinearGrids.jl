using Interpolations, LinearAlgebra

"""
    compute_coons_patch_interp(C0::Vector{NTuple{2,T}}, C1::Vector{NTuple{2,T}}, D0::Vector{NTuple{2,T}}, D1::Vector{NTuple{2,T}}; tol=1e-6)

Create a 2D Coon's patch using 4 vectors of (x,y) coordinates. 

# Arguments
 - `C0`: Vector of points defining the C0 curve, e.g. "bottom"
 - `C1`: Vector of points defining the C1 curve, e.g. "top"
 - `D0`: Vector of points defining the D0 curve, e.g. "left"
 - `D1`: Vector of points defining the D1 curve, e.g. "right"

 # Returns
 - `x`: Matrix{T} of x coordinates
 - `y`: Matrix{T} of y coordinates

"""
function coons_patch(
  C0::Vector{NTuple{2,T}},
  C1::Vector{NTuple{2,T}},
  D0::Vector{NTuple{2,T}},
  D1::Vector{NTuple{2,T}},
) where {T}
  nu = length(C0)
  nv = length(D0)

  @assert length(C1) == length(C0)
  @assert length(D1) == length(D0)

  # Check corners
  corners_ok = check_corners(C0, C1, D0, D1)
  if !corners_ok
    error("Boundary curves do not meet at corners within tolerance.")
  end

  # Extract x and y from each boundary
  C0x, C0y = getxy(C0)
  C1x, C1y = getxy(C1)
  D0x, D0y = getxy(D0)
  D1x, D1y = getxy(D1)

  # Create interpolants
  ugrid = range(0, 1; length=nu)
  vgrid = range(0, 1; length=nv)

  itp_C0x = LinearInterpolation(ugrid, C0x)
  itp_C0y = LinearInterpolation(ugrid, C0y)
  itp_C1x = LinearInterpolation(ugrid, C1x)
  itp_C1y = LinearInterpolation(ugrid, C1y)

  itp_D0x = LinearInterpolation(vgrid, D0x)
  itp_D0y = LinearInterpolation(vgrid, D0y)
  itp_D1x = LinearInterpolation(vgrid, D1x)
  itp_D1y = LinearInterpolation(vgrid, D1y)

  x = zeros(eltype(first(C0)), nv, nu)
  y = zeros(eltype(first(C0)), nv, nu)

  for i in 1:nv
    v = vgrid[i]
    for j in 1:nu
      u = ugrid[j]

      C0uv = (itp_C0x(u), itp_C0y(u))
      C1uv = (itp_C1x(u), itp_C1y(u))
      D0uv = (itp_D0x(v), itp_D0y(v))
      D1uv = (itp_D1x(v), itp_D1y(v))

      Lc = lerp(C0uv, C1uv, v)
      Ld = lerp(D0uv, D1uv, u)

      P00 = C0[1]
      P10 = C0[end]
      P01 = C1[1]
      P11 = C1[end]

      B = bilinear_interp(P00, P10, P01, P11, u, v)

      coordinate = Lc .+ Ld .- B

      x[i, j] = coordinate[1]
      y[i, j] = coordinate[2]
    end
  end

  return x, y
end

# Helpers
function getxy(points)
  xs = [p[1] for p in points]
  ys = [p[2] for p in points]
  return xs, ys
end

# function lerp(
#   p0::NTuple{N,Unitful.Length}, p1::NTuple{N,Unitful.Length}, t::Unitful.Length
# ) where {N}
#   _one = 1.0 * unit(t)
#   return (
#     (_one - t) * p0[1] + t * p1[1], # x
#     (_one - t) * p0[2] + t * p1[2], # y
#   )
# end

function lerp(p0, p1, t)
  return (
    (1 - t) * p0[1] + t * p1[1], # x
    (1 - t) * p0[2] + t * p1[2], # y
  )
end

function check_corners(C0, C1, D0, D1, rtol=1e-6)
  corners = [
    ("C0[begin]", "D0[begin]", C0[begin], D0[begin]),
    ("C0[end]", "D1[begin]", C0[end], D1[begin]),
    ("C1[begin]", "D0[end]", C1[begin], D0[end]),
    ("C1[end]", "D1[end]", C1[end], D1[end]),
  ]

  for (name1, name2, p1, p2) in corners
    dist = norm(p1 .- p2)
    if !all(isapprox.(p1, p2, rtol=rtol))
      @warn "Corner mismatch between $name1 and $name2: distance = $dist"
      return false
    end
  end
  return true
end

function bilinear_interp(P00, P10, P01, P11, u, v)
  x =
    (1 - u) * (1 - v) * P00[1] +
    u * (1 - v) * P10[1] +
    (1 - u) * v * P01[1] +
    u * v * P11[1]
  y =
    (1 - u) * (1 - v) * P00[2] +
    u * (1 - v) * P10[2] +
    (1 - u) * v * P01[2] +
    u * v * P11[2]
  return (x, y)
end

# function bilinear_interp(P00, P10, P01, P11, u::Unitful.Length, v::Unitful.Length)
#   _one = 1 * unit(u)
#   x =
#     (_one - u) * (_one - v) * P00[1] +
#     u * (_one - v) * P10[1] +
#     (_one - u) * v * P01[1] +
#     u * v * P11[1]
#   y =
#     (_one - u) * (_one - v) * P00[2] +
#     u * (_one - v) * P10[2] +
#     (_one - u) * v * P01[2] +
#     u * v * P11[2]
#   return (x, y)
# end

# Define discrete boundary curves
function linspace2d(p1, p2, n)
  [(1 - t) .* p1 .+ t .* p2 for t in range(0, 1; length=n)]
end

function coons_patch(c0, c1, d0, d1)
  Lc(s, t) = (1 - t) * c0(s) + t * c1(s)
  Ld(s, t) = (1 - s) * d0(t) + s * d1(t)

  function B(s, t)
    return c0(0) * (1 - s) * (1 - t) +
           c0(1) * s * (1 - t) +
           c1(0) * (1 - s) * t +
           c1(1) * s * t
  end

  C(s, t) = Lc(s, t) + Ld(s, t) - B(s, t)

  return C
end
