function rectilinear_nodes_2d(x0, x1, y0, y1, ni, nj; halo_coords_included=false)
  Δx = (x1 - x0) / ni
  Δy = (y1 - y0) / nj

  if halo_coords_included
    xnodes = collect(range(x0, x1; length=ni + 1))
    ynodes = collect(range(y0, y1; length=nj + 1))
  else
    xnodes = collect(range(x0, x1; length=ni + 1))
    ynodes = collect(range(y0, y1; length=nj + 1))
  end

  x = [xnodes[i] for i in eachindex(xnodes), j in eachindex(ynodes)]
  y = [ynodes[j] for i in eachindex(xnodes), j in eachindex(ynodes)]
  return x, y, Δx, Δy
end

function assert_rectilinear_metrics_2d(mesh, domain, Δx, Δy)
  cm = cell_metrics(mesh)
  fm = face_metrics(mesh)

  J = Δx * Δy
  ξx = inv(Δx)
  ηy = inv(Δy)
  ξ̂x = J * ξx
  η̂y = J * ηy

  @test all(cm.forward[idx].J ≈ J for idx in domain)
  @test all(cm.forward[idx][1, 1] ≈ Δx for idx in domain)
  @test all(cm.forward[idx][1, 2] ≈ 0.0 for idx in domain)
  @test all(cm.forward[idx][2, 1] ≈ 0.0 for idx in domain)
  @test all(cm.forward[idx][2, 2] ≈ Δy for idx in domain)

  @test all(cm.inverse[idx][1, 1] ≈ ξx for idx in domain)
  @test all(cm.inverse[idx][1, 2] ≈ 0.0 for idx in domain)
  @test all(cm.inverse[idx][2, 1] ≈ 0.0 for idx in domain)
  @test all(cm.inverse[idx][2, 2] ≈ ηy for idx in domain)

  iaxis, jaxis = (1, 2)
  i₊½_domain = expand(domain, iaxis, -1)
  j₊½_domain = expand(domain, jaxis, -1)

  @test all(fm[1].forward[idx].J ≈ J for idx in i₊½_domain)
  @test all(fm[2].forward[idx].J ≈ J for idx in j₊½_domain)
  @test all(fm[1].forward[idx][1, 1] ≈ Δx for idx in i₊½_domain)
  @test all(fm[2].forward[idx][1, 1] ≈ Δx for idx in j₊½_domain)
  @test all(fm[1].forward[idx][2, 2] ≈ Δy for idx in i₊½_domain)
  @test all(fm[2].forward[idx][2, 2] ≈ Δy for idx in j₊½_domain)

  @test all(fm[1].inverse[idx][1, 1] ≈ ξx for idx in i₊½_domain)
  @test all(fm[2].inverse[idx][1, 1] ≈ ξx for idx in j₊½_domain)
  @test all(fm[1].inverse[idx][2, 2] ≈ ηy for idx in i₊½_domain)
  @test all(fm[2].inverse[idx][2, 2] ≈ ηy for idx in j₊½_domain)

  @test all(fm[1].conserved[idx][1, 1] ≈ ξ̂x for idx in i₊½_domain)
  @test all(fm[2].conserved[idx][1, 1] ≈ ξ̂x for idx in j₊½_domain)
  @test all(fm[1].conserved[idx][1, 2] ≈ 0.0 for idx in i₊½_domain)
  @test all(fm[2].conserved[idx][1, 2] ≈ 0.0 for idx in j₊½_domain)
  @test all(fm[1].conserved[idx][2, 1] ≈ 0.0 for idx in i₊½_domain)
  @test all(fm[2].conserved[idx][2, 1] ≈ 0.0 for idx in j₊½_domain)
  @test all(fm[1].conserved[idx][2, 2] ≈ η̂y for idx in i₊½_domain)
  @test all(fm[2].conserved[idx][2, 2] ≈ η̂y for idx in j₊½_domain)
end

@testset "2D Uniform Mesh" begin
  x0, x1 = (0, 4)
  y0, y1 = (1, 9)
  ni, nj = (40, 80)
  x, y, Δx, Δy = rectilinear_nodes_2d(x0, x1, y0, y1, ni, nj)

  mesh = DiscreteGrid(x, y, :MEG6)
  domain = mesh.iterators.cell.domain

  @test mesh.iterators.cell.full == CartesianIndices((50, 90))
  @test mesh.iterators.cell.domain == CartesianIndices((6:45, 6:85))
  @test mesh.iterators.node.full == CartesianIndices((51, 91))
  @test mesh.iterators.node.domain == CartesianIndices((6:46, 6:86))

  assert_rectilinear_metrics_2d(mesh, domain, Δx, Δy)

  ilo_c = mesh.nhalo + 1
  jlo_c = mesh.nhalo + 1
  @test coord(mesh, (ilo_c + 1, jlo_c + 1)) == [0.1, 1.1]
  @test centroid(mesh, (ilo_c, jlo_c)) == [0.05, 1.05]

  xn, yn = coords(mesh)
  @test size(xn) == (41, 81)
  @test size(yn) == (41, 81)

  xc, yc = centroids(mesh)
  @test size(xc) == (40, 80)
  @test size(yc) == (40, 80)
end

@testset "2D Rectangular Mesh" begin
  ni, nj = (40, 80)
  x0, x1 = (0, 2)
  y0, y1 = (1, 3)
  x, y, Δx, Δy = rectilinear_nodes_2d(x0, x1, y0, y1, ni, nj)
  mesh = DiscreteGrid(x, y, :MEG6)
  domain = mesh.iterators.cell.domain

  @test mesh.iterators.cell.full == CartesianIndices((50, 90))
  @test mesh.iterators.cell.domain == CartesianIndices((6:45, 6:85))
  @test mesh.iterators.node.full == CartesianIndices((51, 91))
  @test mesh.iterators.node.domain == CartesianIndices((6:46, 6:86))

  assert_rectilinear_metrics_2d(mesh, domain, Δx, Δy)

  ilo_c = mesh.nhalo + 1
  jlo_c = mesh.nhalo + 1
  @test coord(mesh, (ilo_c + 1, jlo_c + 1)) == [0.05, 1.025]
  @test centroid(mesh, (ilo_c, jlo_c)) == [0.025, 1.0125]

  xn, yn = coords(mesh)
  @test size(xn) == (41, 81)
  @test size(yn) == (41, 81)

  xc, yc = centroids(mesh)
  @test size(xc) == (40, 80)
  @test size(yc) == (40, 80)
end

@testset "2D Rectangular Mesh -- Halo Geometry Defined" begin
  ni, nj = (40, 80)
  x0, x1 = (0, 2)
  y0, y1 = (1, 3)
  x, y, Δx, Δy = rectilinear_nodes_2d(x0, x1, y0, y1, ni, nj; halo_coords_included=true)
  mesh = DiscreteGrid(x, y, :MEG6; halo_coords_included=true)
  domain = mesh.iterators.cell.full

  @test mesh.iterators.cell.full == CartesianIndices((40, 80))
  @test mesh.iterators.cell.domain == CartesianIndices((6:35, 6:75))
  @test mesh.iterators.node.full == CartesianIndices((41, 81))
  @test mesh.iterators.node.domain == CartesianIndices((6:36, 6:76))

  assert_rectilinear_metrics_2d(mesh, domain, Δx, Δy)

  ilo_c = mesh.nhalo + 1
  jlo_c = mesh.nhalo + 1
  @test coord(mesh, (ilo_c + 1, jlo_c + 1)) == [0.3, 1.15]
  @test centroid(mesh, (ilo_c, jlo_c)) == [0.275, 1.1375]

  xn, yn = coords(mesh)
  @test size(xn) == (31, 71)
  @test size(yn) == (31, 71)

  xc, yc = centroids(mesh)
  @test size(xc) == (30, 70)
  @test size(yc) == (30, 70)
end

@testset "2D Wavy Mesh GCL" begin
  function wavy_grid(nx, ny)
    x0, x1 = (0, 1)
    y0, y1 = (0, 1)
    a0 = 0.1

    function x(i, j)
      x1d = x0 + (x1 - x0) * ((i - 1) / (nx - 1))
      y1d = y0 + (y1 - y0) * ((j - 1) / (ny - 1))
      return x1d + a0 * sin(2 * pi * x1d) * sin(2 * pi * y1d)
    end

    function y(i, j)
      x1d = x0 + (x1 - x0) * ((i - 1) / (nx - 1))
      y1d = y0 + (y1 - y0) * ((j - 1) / (ny - 1))
      return y1d + a0 * sin(2 * pi * x1d) * sin(2 * pi * y1d)
    end

    X = zeros(nx, ny)
    Y = zeros(nx, ny)
    for j in 1:ny
      for i in 1:nx
        X[i, j] = x(i, j)
        Y[i, j] = y(i, j)
      end
    end

    return X, Y
  end

  ni, nj = (41, 41)
  x, y = wavy_grid(ni, nj)
  mesh = DiscreteGrid(x, y, :MEG6)

  domain = mesh.iterators.cell.domain
  I₁, I₂ = CurvilinearGrids.GridTypes.gcl(face_metrics(mesh), domain)
  gcl_identities = (all(abs.(I₁[domain]) .< 5e-15), all(abs.(I₂[domain]) .< 5e-15))
  max_vals = (maximum(abs, I₁[domain]), maximum(abs, I₂[domain]))
  @test all(gcl_identities)

  save_vtk(coords(mesh), "wavy")
  nothing
end

@testset "2D Wavy Mesh GCL -- Halo Geometry Defined" begin
  using CurvilinearGrids

  function wavy_grid(nx, ny)
    x0, x1 = (0, 1)
    y0, y1 = (0, 1)
    a0 = 0.1

    function x(i, j)
      x1d = x0 + (x1 - x0) * ((i - 1) / (nx - 1))
      y1d = y0 + (y1 - y0) * ((j - 1) / (ny - 1))
      return x1d + a0 * sin(2 * pi * x1d) * sin(2 * pi * y1d)
    end

    function y(i, j)
      x1d = x0 + (x1 - x0) * ((i - 1) / (nx - 1))
      y1d = y0 + (y1 - y0) * ((j - 1) / (ny - 1))
      return y1d + a0 * sin(2 * pi * x1d) * sin(2 * pi * y1d)
    end

    X = zeros(nx, ny)
    Y = zeros(nx, ny)
    for j in 1:ny
      for i in 1:nx
        X[i, j] = x(i, j)
        Y[i, j] = y(i, j)
      end
    end

    return X, Y
  end

  ni, nj = (41, 41)
  x, y = wavy_grid(ni, nj)
  mesh = DiscreteGrid(x, y, :MEG6; halo_coords_included=true)

  domain = mesh.iterators.cell.domain
  I₁, I₂ = CurvilinearGrids.GridTypes.gcl(face_metrics(mesh), domain)
  gcl_identities = (all(abs.(I₁[domain]) .< 5e-15), all(abs.(I₂[domain]) .< 5e-15))
  max_vals = (maximum(abs, I₁[domain]), maximum(abs, I₂[domain]))
  @test all(gcl_identities)
end
