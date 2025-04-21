using Unitful

function reynolds_number(
  ρ::Unitful.Density, U::Unitful.Velocity, L::Unitful.Length, μ::Unitful.DynamicViscosity
)
  (ρ * U * L / μ) |> NoUnits
end
skin_friction_coefficient(Re) = (2log10(Re) - 0.65)^(-2.3)

"""
    estimate_wall_distance(
  ρ::Unitful.Density,
  U::Unitful.Velocity,
  L::Unitful.Length,
  μ::Unitful.DynamicViscosity,
  y_plus,
)

Estimate the wall-thickness needed to capture a turbulent boundary layer
based on the input conditions. This returns (reynolds_number, wall_thickness).
"""
function estimate_wall_distance(
  ρ::Unitful.Density,
  U::Unitful.Velocity,
  L::Unitful.Length,
  μ::Unitful.DynamicViscosity,
  y_plus,
)
  Re = reynolds_number(ρ, U, L, μ)

  if Re > 1e9
    @warn """
The Schlichting skin-friction formula used to estimate the local skin-friction
for a turbulent boundary layer on a smooth flat plate is only valid for Reynolds numbers < 10^9,
and the calculated value is $(Re) based on your inputs.
"""
  end

  Cf = skin_friction_coefficient(Re)
  τ_wall = Cf * (ρ * U^2) / 2 # wall shear stress
  u_star = sqrt(τ_wall / ρ)

  y = (y_plus * μ) / (ρ * u_star)
  return (; reynolds_number=Re, wall_thickness=y |> Unitful.Length)
end
