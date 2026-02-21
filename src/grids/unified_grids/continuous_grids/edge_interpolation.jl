#-------------------------------------------------------------
# 1D edge interpolation
#-------------------------------------------------------------
function ϕ_iedge(ϕ_eval, t, i::Real, p::NamedTuple)
  ϕ_iedge(ϕ_eval, t, i, p, EdgeInterpolationOrder3(), 1)
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

function ϕ_iedge(ϕ_values, t, i::Real, p::NamedTuple, ::EdgeInterpolationOrder1, Δξ::Real)
  return _edge_reconstruct(
    ϕ_values(t, i, p), ϕ_values(t, i + 1, p), EdgeInterpolationOrder1()
  )
end

function ϕ_iedge(
  ϕ_val_and_derivs, t, i::Real, p::NamedTuple, ::EdgeInterpolationOrder2, Δξ::Real
)
  ϕᵢ, ϕξᵢ = ϕ_val_and_derivs(t, i, p)
  ϕᵢ₊₁, ϕξᵢ₊₁ = ϕ_val_and_derivs(t, i + 1, p)
  return _edge_reconstruct(ϕᵢ, ϕξᵢ, ϕᵢ₊₁, ϕξᵢ₊₁, EdgeInterpolationOrder2(), Δξ)
end

function ϕ_iedge(
  ϕ_val_and_derivs, t, i::Real, p::NamedTuple, ::EdgeInterpolationOrder3, Δξ::Real
)
  ϕᵢ, ϕξᵢ, ϕξξᵢ = ϕ_val_and_derivs(t, i, p)
  ϕᵢ₊₁, ϕξᵢ₊₁, ϕξξᵢ₊₁ = ϕ_val_and_derivs(t, i + 1, p)
  return _edge_reconstruct(
    ϕᵢ, ϕξᵢ, ϕξξᵢ, ϕᵢ₊₁, ϕξᵢ₊₁, ϕξξᵢ₊₁, EdgeInterpolationOrder3(), Δξ
  )
end

#-------------------------------------------------------------
# 2D edge interpolation
#-------------------------------------------------------------
function ϕ_iedge(ϕ_eval, t, i::Real, j::Real, p::NamedTuple)
  ϕ_iedge(ϕ_eval, t, i, j, p, EdgeInterpolationOrder3(), 1)
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
  ϕ_values, t, i::Real, j::Real, p::NamedTuple, ::EdgeInterpolationOrder1, Δξ::Real
)
  return _edge_reconstruct(
    ϕ_values(t, i, j, p), ϕ_values(t, i + 1, j, p), EdgeInterpolationOrder1()
  )
end

function ϕ_iedge(
  ϕ_val_and_derivs, t, i::Real, j::Real, p::NamedTuple, ::EdgeInterpolationOrder2, Δξ::Real
)
  ϕᵢ, ϕξᵢ = ϕ_val_and_derivs(t, i, j, p)
  ϕᵢ₊₁, ϕξᵢ₊₁ = ϕ_val_and_derivs(t, i + 1, j, p)
  return _edge_reconstruct(ϕᵢ, ϕξᵢ, ϕᵢ₊₁, ϕξᵢ₊₁, EdgeInterpolationOrder2(), Δξ)
end

function ϕ_iedge(
  ϕ_val_and_derivs, t, i::Real, j::Real, p::NamedTuple, ::EdgeInterpolationOrder3, Δξ::Real
)
  ϕᵢ, ϕξᵢ, ϕξξᵢ = ϕ_val_and_derivs(t, i, j, p)
  ϕᵢ₊₁, ϕξᵢ₊₁, ϕξξᵢ₊₁ = ϕ_val_and_derivs(t, i + 1, j, p)
  return _edge_reconstruct(
    ϕᵢ, ϕξᵢ, ϕξξᵢ, ϕᵢ₊₁, ϕξᵢ₊₁, ϕξξᵢ₊₁, EdgeInterpolationOrder3(), Δξ
  )
end

function ϕ_jedge(ϕ_eval, t, i::Real, j::Real, p::NamedTuple)
  ϕ_jedge(ϕ_eval, t, i, j, p, EdgeInterpolationOrder3(), 1)
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
  ϕ_values, t, i::Real, j::Real, p::NamedTuple, ::EdgeInterpolationOrder1, Δξ::Real
)
  return _edge_reconstruct(
    ϕ_values(t, i, j, p), ϕ_values(t, i, j + 1, p), EdgeInterpolationOrder1()
  )
end

function ϕ_jedge(
  ϕ_val_and_derivs, t, i::Real, j::Real, p::NamedTuple, ::EdgeInterpolationOrder2, Δξ::Real
)
  ϕⱼ, ϕηⱼ = ϕ_val_and_derivs(t, i, j, p)
  ϕⱼ₊₁, ϕηⱼ₊₁ = ϕ_val_and_derivs(t, i, j + 1, p)
  return _edge_reconstruct(ϕⱼ, ϕηⱼ, ϕⱼ₊₁, ϕηⱼ₊₁, EdgeInterpolationOrder2(), Δξ)
end

function ϕ_jedge(
  ϕ_val_and_derivs, t, i::Real, j::Real, p::NamedTuple, ::EdgeInterpolationOrder3, Δξ::Real
)
  ϕⱼ, ϕηⱼ, ϕηηⱼ = ϕ_val_and_derivs(t, i, j, p)
  ϕⱼ₊₁, ϕηⱼ₊₁, ϕηηⱼ₊₁ = ϕ_val_and_derivs(t, i, j + 1, p)
  return _edge_reconstruct(
    ϕⱼ, ϕηⱼ, ϕηηⱼ, ϕⱼ₊₁, ϕηⱼ₊₁, ϕηηⱼ₊₁, EdgeInterpolationOrder3(), Δξ
  )
end

#-------------------------------------------------------------
# 3D edge interpolation
#-------------------------------------------------------------
function ϕ_iedge(ϕ_eval, t, i::Real, j::Real, k::Real, p::NamedTuple)
  ϕ_iedge(ϕ_eval, t, i, j, k, p, EdgeInterpolationOrder3(), 1)
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
  ϕ_values, t, i::Real, j::Real, k::Real, p::NamedTuple, ::EdgeInterpolationOrder1, Δξ::Real
)
  return _edge_reconstruct(
    ϕ_values(t, i, j, k, p), ϕ_values(t, i + 1, j, k, p), EdgeInterpolationOrder1()
  )
end

function ϕ_iedge(
  ϕ_val_and_derivs,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  ::EdgeInterpolationOrder2,
  Δξ::Real,
)
  ϕᵢ, ϕξᵢ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕᵢ₊₁, ϕξᵢ₊₁ = ϕ_val_and_derivs(t, i + 1, j, k, p)
  return _edge_reconstruct(ϕᵢ, ϕξᵢ, ϕᵢ₊₁, ϕξᵢ₊₁, EdgeInterpolationOrder2(), Δξ)
end

function ϕ_iedge(
  ϕ_val_and_derivs,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  ::EdgeInterpolationOrder3,
  Δξ::Real,
)
  ϕᵢ, ϕξᵢ, ϕξξᵢ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕᵢ₊₁, ϕξᵢ₊₁, ϕξξᵢ₊₁ = ϕ_val_and_derivs(t, i + 1, j, k, p)
  return _edge_reconstruct(
    ϕᵢ, ϕξᵢ, ϕξξᵢ, ϕᵢ₊₁, ϕξᵢ₊₁, ϕξξᵢ₊₁, EdgeInterpolationOrder3(), Δξ
  )
end

function ϕ_jedge(ϕ_eval, t, i::Real, j::Real, k::Real, p::NamedTuple)
  ϕ_jedge(ϕ_eval, t, i, j, k, p, EdgeInterpolationOrder3(), 1)
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
  ϕ_values, t, i::Real, j::Real, k::Real, p::NamedTuple, ::EdgeInterpolationOrder1, Δξ::Real
)
  return _edge_reconstruct(
    ϕ_values(t, i, j, k, p), ϕ_values(t, i, j + 1, k, p), EdgeInterpolationOrder1()
  )
end

function ϕ_jedge(
  ϕ_val_and_derivs,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  ::EdgeInterpolationOrder2,
  Δξ::Real,
)
  ϕⱼ, ϕηⱼ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕⱼ₊₁, ϕηⱼ₊₁ = ϕ_val_and_derivs(t, i, j + 1, k, p)
  return _edge_reconstruct(ϕⱼ, ϕηⱼ, ϕⱼ₊₁, ϕηⱼ₊₁, EdgeInterpolationOrder2(), Δξ)
end

function ϕ_jedge(
  ϕ_val_and_derivs,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  ::EdgeInterpolationOrder3,
  Δξ::Real,
)
  ϕⱼ, ϕηⱼ, ϕηηⱼ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕⱼ₊₁, ϕηⱼ₊₁, ϕηηⱼ₊₁ = ϕ_val_and_derivs(t, i, j + 1, k, p)
  return _edge_reconstruct(
    ϕⱼ, ϕηⱼ, ϕηηⱼ, ϕⱼ₊₁, ϕηⱼ₊₁, ϕηηⱼ₊₁, EdgeInterpolationOrder3(), Δξ
  )
end

function ϕ_kedge(ϕ_eval, t, i::Real, j::Real, k::Real, p::NamedTuple)
  ϕ_kedge(ϕ_eval, t, i, j, k, p, EdgeInterpolationOrder3(), 1)
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
  ϕ_values, t, i::Real, j::Real, k::Real, p::NamedTuple, ::EdgeInterpolationOrder1, Δξ::Real
)
  return _edge_reconstruct(
    ϕ_values(t, i, j, k, p), ϕ_values(t, i, j, k + 1, p), EdgeInterpolationOrder1()
  )
end

function ϕ_kedge(
  ϕ_val_and_derivs,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  ::EdgeInterpolationOrder2,
  Δξ::Real,
)
  ϕₖ, ϕζₖ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕₖ₊₁, ϕζₖ₊₁ = ϕ_val_and_derivs(t, i, j, k + 1, p)
  return _edge_reconstruct(ϕₖ, ϕζₖ, ϕₖ₊₁, ϕζₖ₊₁, EdgeInterpolationOrder2(), Δξ)
end

function ϕ_kedge(
  ϕ_val_and_derivs,
  t,
  i::Real,
  j::Real,
  k::Real,
  p::NamedTuple,
  ::EdgeInterpolationOrder3,
  Δξ::Real,
)
  ϕₖ, ϕζₖ, ϕζζₖ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕₖ₊₁, ϕζₖ₊₁, ϕζζₖ₊₁ = ϕ_val_and_derivs(t, i, j, k + 1, p)
  return _edge_reconstruct(
    ϕₖ, ϕζₖ, ϕζζₖ, ϕₖ₊₁, ϕζₖ₊₁, ϕζζₖ₊₁, EdgeInterpolationOrder3(), Δξ
  )
end
