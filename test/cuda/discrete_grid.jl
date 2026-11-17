# Run explicitly from an environment with CUDA installed and a working device.
using CurvilinearGrids, CUDA, Adapt, KernelAbstractions, Test
CUDA.functional() || error("a working CUDA device is required; no host fallback")
CUDA.allowscalar(false)

host_arrays(a::AbstractArray) = Array(a)
host_arrays(a::Union{Tuple,NamedTuple}) = map(host_arrays,a)

@testset "DiscreteGrid CUDA storage and lifecycle" begin
  for T in (Float32,Float64), N in 1:3, mode in (:eager,:lazy,:off)
    indices = CartesianIndices(ntuple(d->6-d,N))
    nodes = ntuple(axis -> [T(axis+1)*I[axis]+T(0.01)*sum(T(I[d])^2 for d in 1:N) for I in indices],N)
    cpu = DiscreteGrid(nodes...,2;T,cache_mode=mode)
    gpu = Adapt.adapt(CuArray,cpu)
    @test all(a -> a isa CuArray,gpu.node_coordinates)
    @test gpu.backend isa typeof(CUDA.CUDABackend())
    @test host_arrays(gpu.node_coordinates) == cpu.node_coordinates
    @test host_arrays(gpu.centroid_coordinates) == cpu.centroid_coordinates
    @test host_arrays(gpu.face_coordinates) == cpu.face_coordinates
    @test host_arrays(cell_metrics(gpu)) == cell_metrics(cpu)
    @test host_arrays(face_metrics(gpu)) == face_metrics(cpu)
    direct = DiscreteGrid(nodes...,2;T,backend=CUDA.CUDABackend(),cache_mode=mode)
    @test host_arrays(cell_metrics(direct)) == cell_metrics(cpu)
    update!(gpu,T(2),(;))
    @test cpu.state[].t == zero(T)
    @test gpu.state[].t == T(2)
    @test host_arrays(cell_metrics(gpu)) == cell_metrics(cpu)
    @test host_arrays(face_metrics(gpu)) == face_metrics(cpu)
    returned = Adapt.adapt(Array,gpu)
    @test returned.backend isa CPU
    @test returned.node_coordinates == cpu.node_coordinates
    @test cell_metrics(returned) == cell_metrics(cpu)
  end
end

# Preserve explicitly requested backend options during direct construction.
configured = CUDA.CUDABackend(; always_inline=true)
configured_grid = DiscreteGrid(collect(0.0:4.0), 2; backend=configured)
@test configured_grid.backend === configured
