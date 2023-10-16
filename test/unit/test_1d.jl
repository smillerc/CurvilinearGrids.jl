
@testitem "1D Mesh - Rectlinear Mesh" begin
  include("common.jl")

  function LinearSpacing((x0, x1), ni, nhalo)
    _x(i) = x0 + (x1 - x0) * ((i - 1) / (ni - 1))
    x(i) = _x(i - nhalo)
    return x
  end

  function getmesh()
    nhalo = 0
    ni = 5
    x0, x1 = (0.0, 2.0)
    x = LinearSpacing((x0, x1), ni, nhalo)
    return CurvilinearGrid1D(x, ni, nhalo)
  end

  m = getmesh()
  @test metrics(m, 1) == (ξx=2.0, ξt=0.0)
  @test conservative_metrics(m, 1) == (ξ̂x=4.0, ξt=0.0)
  @test jacobian_matrix(m, 1) == @SMatrix [0.5]
  @test jacobian(m, 1) == 0.5
  @test inv(jacobian_matrix(m, 1)) == @SMatrix [2.0]
  @test inv(jacobian(m, 1)) == 2.0
end