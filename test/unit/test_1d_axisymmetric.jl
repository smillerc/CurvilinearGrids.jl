
using CurvilinearGrids
using Test
using StaticArrays
using LinearAlgebra
using StructArrays
using WriteVTK
using UnPack

function save_vtk(mesh)
  fn = "mesh"
  @info "Writing to $fn.vti"

  domain = mesh.iterators.node.domain
  @unpack x, y, z = StructArrays.components(mesh.node_coordinates.xyz)

  @views vtk_grid(fn, x[domain], y[domain], z[domain]) do vtk
  end
end

@testset "1D Spherical Mesh" begin
  function getmesh()
    nhalo = 2
    ni = 6
    x0, x1 = (0.0, 5.0)
    x(i) = x0 + (x1 - x0) * ((i - 1) / (ni - 1))
    return SphericalGrid1D(x, ni, nhalo)
  end

  mesh = getmesh()
  save_vtk(mesh)
  nothing

  @test mesh.node_coordinates.r == [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
  @test mesh.centroid_coordinates.r == [-1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]

  nh = 2
  mesh._jacobian_matrix_func(5, 1, 1) == jacobian_matrix(mesh, (5 + nh, 1 + nh, 1 + nh))

  J = mesh._jacobian_matrix_func(5, 1.5, 1.5)
  Jexact = @SMatrix [
    1.0 0.0 0.0
    0.0 0.0 2pi
    0.0 -2pi 0.0
  ]
  @test J == Jexact

  i0 = 3
  r0 = 2.0
  θ0, θ1 = (1, 2)
  ϕ0, ϕ1 = (1, 2)

  @test mesh._coordinate_funcs.XYZ(i0, θ0, ϕ0) ≈ [1, -1, √2]
  @test mesh._coordinate_funcs.XYZ(i0, θ1, ϕ0) ≈ [1, -1, -√2]
  @test mesh._coordinate_funcs.XYZ(i0, θ1, ϕ1) ≈ [1, 1, -√2]
  @test mesh._coordinate_funcs.XYZ(i0, θ0, ϕ1) ≈ [1, 1, √2]

  @test mesh._coordinate_funcs.XYZ(i0, θ0 + 0.5, ϕ1) ≈ [√2, √2, 0]
  @test mesh._coordinate_funcs.XYZ(i0, θ0 + 0.5, ϕ0 + 0.5) ≈ [2, 0, 0]

  # for idx in mesh.iterators.cell.domain
  #   i, = idx.I .+ 0.5
  #   m_i₊½ = conservative_metrics(mesh, i + 0.5)
  #   m_i₋½ = conservative_metrics(mesh, i - 0.5)

  #   I₁ = (m_i₊½.ξ̂.x - m_i₋½.ξ̂.x)
  #   # @show I₁, m_i₊½.ξ̂.x, m_i₋½.ξ̂.x
  #   I₁_passes = abs(I₁) < eps()
  #   if !(I₁_passes)
  #     break
  #   end
  # end
  # @test I₁_passes

  conserved_metrics_pass = false
  domain = mesh.iterators.cell.domain
  # for idx in expand(domain, -1)
  for idx in domain
    i, j, k = idx.I
    I₁ = (
      (mesh.edge_metrics.i₊½.ξ̂.x[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.x[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.x[i, j, k] - mesh.edge_metrics.j₊½.η̂.x[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.x[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.x[i, j, k - 1])
    )
    I₂ = (
      (mesh.edge_metrics.i₊½.ξ̂.y[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.y[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.y[i, j, k] - mesh.edge_metrics.j₊½.η̂.y[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.y[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.y[i, j, k - 1])
    )
    I₃ = (
      (mesh.edge_metrics.i₊½.ξ̂.z[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.z[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.z[i, j, k] - mesh.edge_metrics.j₊½.η̂.z[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.z[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.z[i, j, k - 1])
    )

    I₁ = I₁ * (abs(I₁) > 5e-15)
    I₂ = I₂ * (abs(I₂) > 5e-15)
    I₃ = I₃ * (abs(I₃) > 5e-15)
    conserved_metrics_pass = iszero(I₁) && iszero(I₂) && iszero(I₃)
    if !conserved_metrics_pass
      break
    end
  end
  @test conserved_metrics_pass
end

@testset "1D Cylindrical Mesh" begin
  function getmesh()
    nhalo = 2
    ni = 6
    x0, x1 = (0.0, 5.0)
    x(i) = x0 + (x1 - x0) * ((i - 1) / (ni - 1))
    return CylindricalGrid1D(x, ni, nhalo)
  end

  mesh = getmesh()
  save_vtk(mesh)
  nothing

  @test mesh.node_coordinates.r == [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
  @test mesh.centroid_coordinates.r == [-1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]

  # nh = 2
  # mesh._jacobian_matrix_func(5, 1, 1) == jacobian_matrix(mesh, (5 + nh, 1 + nh, 1 + nh))

  # J = mesh._jacobian_matrix_func(5, 1.5, 1.5)
  # Jexact = @SMatrix [
  #   1.0 0.0 0.0
  #   0.0 0.0 2pi
  #   0.0 -2pi 0.0
  # ]
  # @test J == Jexact

  # i0 = 3
  # r0 = 2.0
  # θ0, θ1 = (1, 2)
  # ϕ0, ϕ1 = (1, 2)

  # @test mesh._coordinate_funcs.XYZ(i0, θ0, ϕ0) ≈ [1, -1, √2]
  # @test mesh._coordinate_funcs.XYZ(i0, θ1, ϕ0) ≈ [1, -1, -√2]
  # @test mesh._coordinate_funcs.XYZ(i0, θ1, ϕ1) ≈ [1, 1, -√2]
  # @test mesh._coordinate_funcs.XYZ(i0, θ0, ϕ1) ≈ [1, 1, √2]

  # @test mesh._coordinate_funcs.XYZ(i0, θ0 + 0.5, ϕ1) ≈ [√2, √2, 0]
  # @test mesh._coordinate_funcs.XYZ(i0, θ0 + 0.5, ϕ0 + 0.5) ≈ [2, 0, 0]

  # # for idx in mesh.iterators.cell.domain
  # #   i, = idx.I .+ 0.5
  # #   m_i₊½ = conservative_metrics(mesh, i + 0.5)
  # #   m_i₋½ = conservative_metrics(mesh, i - 0.5)

  # #   I₁ = (m_i₊½.ξ̂.x - m_i₋½.ξ̂.x)
  # #   # @show I₁, m_i₊½.ξ̂.x, m_i₋½.ξ̂.x
  # #   I₁_passes = abs(I₁) < eps()
  # #   if !(I₁_passes)
  # #     break
  # #   end
  # # end
  # # @test I₁_passes

  conserved_metrics_pass = false
  domain = mesh.iterators.cell.domain
  # for idx in expand(domain, -1)
  for idx in domain
    i, j, k = idx.I
    I₁ = (
      (mesh.edge_metrics.i₊½.ξ̂.x[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.x[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.x[i, j, k] - mesh.edge_metrics.j₊½.η̂.x[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.x[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.x[i, j, k - 1])
    )
    I₂ = (
      (mesh.edge_metrics.i₊½.ξ̂.y[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.y[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.y[i, j, k] - mesh.edge_metrics.j₊½.η̂.y[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.y[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.y[i, j, k - 1])
    )
    I₃ = (
      (mesh.edge_metrics.i₊½.ξ̂.z[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.z[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.z[i, j, k] - mesh.edge_metrics.j₊½.η̂.z[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.z[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.z[i, j, k - 1])
    )

    I₁ = I₁ * (abs(I₁) > 5e-15)
    I₂ = I₂ * (abs(I₂) > 5e-15)
    I₃ = I₃ * (abs(I₃) > 5e-15)
    conserved_metrics_pass = iszero(I₁) && iszero(I₂) && iszero(I₃)
    if !conserved_metrics_pass
      break
    end
  end
  @test conserved_metrics_pass
end
