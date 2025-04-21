using Unitful

reynolds_number(ρ, U, L, μ) = (ρ * U * L / μ) |> NoUnits
skin_friction_coefficient(Re) = (2log10(Re) - 0.65)^(-2.3)

function estimate_wall_distance(
  ρ::Unitful.Density,
  U::Unitful.Velocity,
  L::Unitful.Length,
  μ::Unitful.DynamicViscosity,
  y_plus,
)
  Re = reynolds_number(ρ, U, L, μ)
  Cf = skin_friction_coefficient(Re)
  τ_wall = Cf * (ρ * U^2) / 2
  u_star = sqrt(τ_wall / ρ)

  y = (y_plus * μ) / (ρ * u_star)
  return y |> Unitful.Length
end
