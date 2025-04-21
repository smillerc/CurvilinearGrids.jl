using Unitful
using CurvilinearGrids
using CurvilinearGrids:
  estimate_wall_distance, one_sided_stretch, one_sided_with_initial_spacing
using Test

@testset "Wall spacing" begin
  cs = 340.29u"m/s" # sound speed

  ρ = 1.205u"kg/m^3"
  U = 0.2 * cs # set U to Mach 0.2
  L = 1.0u"ft"
  μ = 1.82e-5u"kg/m/s"
  y_plus = 1.0

  # Re = reynolds_number(ρ, U, L, μ)
  Re, wall_thickness = estimate_wall_distance(ρ, U, L, μ, y_plus)

  @test Re ≈ 1.3734403556043962e6
  @test wall_thickness ≈ 5.2716756477660915u"μm"

  xcoords = one_sided_with_initial_spacing((0, 5), 100, 1e-6)
  @test first(diff(xcoords)) ≈ 1e-6

  xcoords = one_sided_with_initial_spacing((0u"ft", 1u"ft"), 100, wall_thickness)
  first(diff(xcoords)) ≈ wall_thickness
end
