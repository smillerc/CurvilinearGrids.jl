@testset "Uniform MappedGrid1D" begin
  function get_uniform_mapping()
    function x(t, i, p)
      @unpack xmin, Δx = p
      return xmin + (i - 1) * Δx
    end

    return x
  end

  function get_uniform_params()
    celldims = (40,)
    xmin, xmax = (0.0, 2.0)

    ni, = celldims
    Δx = (xmax - xmin) / ni

    params = (; xmin, Δx)
    return params, celldims
  end

  params, celldims = get_uniform_params()
  x = get_uniform_mapping()

  backend = AutoForwardDiff()
  mesh = MappedGrid(x, params, celldims, :meg6; backend=CPU(), diff_backend=backend)

  I1 = CurvilinearGrids.GridTypes.gcl(face_metrics(mesh), mesh.iterators.cell.domain)
  @test all(abs.(extrema(I1)) .< eps())

  cm = cell_metrics(mesh)
  cell_volume = 0.05

  @test all([m.J for m in cm.forward] .≈ cell_volume)
  @test all([m[1, 1] for m in cm.forward] .≈ 0.05)
  @test all([m[1, 1] for m in cm.inverse] .≈ 20.0)
  @test all([m[1, 1] * m.J for m in cm.inverse] .≈ 1.0)

  iaxis = 1
  domain = mesh.iterators.cell.domain
  i₊½_domain = expand(domain, iaxis, -1)
  fm = face_metrics(mesh)

  @test all([m[1, 1] for m in fm[1].conserved[i₊½_domain]] .≈ 1.0)
end
