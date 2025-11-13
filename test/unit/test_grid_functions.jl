using Plots

U(t, a, b, c) = (a / 2) * tanh(b * (t - c))
V(t, d, e, f1, f2) = ((d - 1) / 2e) * log(cosh(e * (t - f1)) / cosh(e * (t - f2)))

function R(t, C0, C1, u_funcs::NTuple{NU,<:Function}) where {NU}
  u_contrib = sum((U(t) - U(0)) for U in u_funcs)

  return (t + u_contrib) * C0 + C1
end

function R(
  t, C0, C1, u_funcs::NTuple{NU,<:Function}, v_funcs::NTuple{NV,<:Function}
) where {NU,NV}
  u_contrib = sum((U(t) - U(0)) for U in u_funcs)
  v_contrib = sum((V(t) - V(0)) for V in v_funcs)

  return (t + u_contrib + v_contrib) * C0 + C1
end

begin
  U1(t) = U(t, 5.0, 20.0, 0.5)

  r = [R(i, 0.5, 0.0, (U1,)) for i in 0:0.01:1]
  plot(r; marker=:square)
  #   plot(r; marker=:square)
end

begin
  V1(t) = V(t, 5.0, 20.0, 0.35, 0.65)

  t = 0:0.025:1
  r = [R(i, 0.5, 1.0, (V1,)) for i in t]

  # plot(r, marker=:square)
  plot(t, r, ; marker=:square)

  # plot(r, zeros(size(r)); marker=:circle)
end

begin
  U1(t) = U(t, 1.0, 20.0, 0.4)
  V1(t) = V(t, 4.0, 20.0, 0.75, 1.25)

  r = [R(i, 0.45, 0.0, (U1,), (V1,)) for i in 0:0.01:1]
  # plot(r, marker=:square)
  plot(r; marker=:square)
  plot(r, zeros(size(r)); marker=:circle)
end
