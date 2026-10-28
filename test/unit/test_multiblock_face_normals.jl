using Test, LinearAlgebra, StaticArrays, CurvilinearGrids

# Intent: evaluate interface normals at physical faces, independent of halo width.
# Oracle: the normal to q_axis = s*(1 + sum(a_d*u_d^2)). At face midpoints,
# linear interpolation of each quadratic has the exact midpoint derivative.
# Budget: 2D/3D grids with at most 6 cells per direction; no time integration.
function curved_normal_blocks(::Val{N}, axis, halos, discrete) where {N}
  dims = ntuple(d -> d + 3, N)
  grids = ntuple(2) do block
    p = (; dims, axis, lo=block == 1 ? 0.0 : 0.4, width=block == 1 ? 0.4 : 0.6)
    functions = ntuple(N) do component
      function mapping(t, args...)
        parameters = args[end]
        u = ntuple(d -> (args[d] - 1) / parameters.dims[d], N)
        component != parameters.axis && return u[component]
        s = parameters.lo + parameters.width * u[component]
        return s * (1 + sum(0.03d * u[d]^2 for d in 1:N if d != component))
      end
    end
    g = MappedGrid(
      functions..., p, dims, halos[block];
      coordinate_system=CartesianCS(), basis=CartesianBasis(),
    )
    discrete || return g
    return DiscreteGrid(
      map(Array, g.node_coordinates)..., halos[block];
      halo_coords_included=true, coordinate_system=CartesianCS(), basis=CartesianBasis(),
    )
  end
  return grids, dims
end

@testset "Curved interface normals use face coordinates" begin
  for N in (2, 3), axis in 1:N, halos in ((1, 1), (2, 3)), discrete in (false, true)
    @testset "N=$N axis=$axis halos=$halos discrete=$discrete" begin
      grids, dims = curved_normal_blocks(Val(N), axis, halos, discrete)
      faces = (BlockFace{N}(1, axis, :max), BlockFace{N}(2, axis, :min))
      interface = BlockInterface(faces..., ntuple(identity, N - 1), ntuple(_ -> false, N - 1))
      # Test each normal directly, even when the public constructor rejects the mesh.
      for block in 1:2
        g, face = grids[block], faces[block]
        ranges = g.iterators.cell.domain.indices
        boundary = block == 1 ? last(ranges[axis]) : first(ranges[axis])
        sample_ranges = ntuple(d -> d == axis ? (boundary:boundary) : ranges[d], N)
        for I in CartesianIndices(sample_ranges)
          u = ntuple(d -> (I[d] - halos[block] - 0.5) / dims[d], N)
          raw = SVector{N}(ntuple(d -> d == axis ? 1.0 : -0.4 * 0.06d * u[d], N))
          expected = (block == 1 ? 1 : -1) * normalize(raw)
          actual = CurvilinearGrids.MultiBlockMeshes._face_outward_normal_cartesian(
            g, I, face, Val(N)
          )
          @test actual ≈ expected rtol=2e-12 atol=2e-12
        end
      end
      @test MultiBlockMesh(grids, (interface,)) isa MultiBlockMesh
    end
  end
end
