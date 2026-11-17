using Test
include("discrete_device_test_support.jl")

@testset "Discrete device geometry preserves CPU calculations" begin
  for T in (Float32,Float64), N in 1:3, mode in (:eager,:lazy,:off)
    dims = ntuple(d -> 6-d,N)
    indices = CartesianIndices(dims)
    nodes = ntuple(N) do axis
      [T(axis+1)*T(I[axis]) + T(0.01)*sum(T(I[d])^2 for d in 1:N) for I in indices]
    end
    cpu = DiscreteGrid(nodes...,2;T,cache_mode=mode)
    device = Adapt.adapt(DiscreteTestAdaptor(),cpu)
    @test device.backend isa DiscreteTestBackend
    @test device.state !== cpu.state
    @test device.metric_caches !== cpu.metric_caches
    @test host_data(device.node_coordinates) == cpu.node_coordinates
    @test host_data(device.centroid_coordinates) == cpu.centroid_coordinates
    @test host_data(device.face_coordinates) == cpu.face_coordinates
    @test host_data(cell_metrics(device)) == cell_metrics(cpu)
    @test host_data(face_metrics(device)) == face_metrics(cpu)

    direct = DiscreteGrid(nodes...,2;T,backend=DiscreteTestBackend(),cache_mode=mode)
    @test host_data(direct.node_coordinates) == cpu.node_coordinates
    @test host_data(cell_metrics(direct)) == cell_metrics(cpu)
    @test host_data(face_metrics(direct)) == face_metrics(cpu)
    device_inputs = map(a -> DiscreteTestArray(copy(a)),nodes)
    from_device = DiscreteGrid(device_inputs...,2;T,backend=DiscreteTestBackend(),cache_mode=mode)
    @test host_data(from_device.node_coordinates) == cpu.node_coordinates

    update!(device,T(2),(;))
    @test device.state[].t == T(2)
    @test cpu.state[].t == zero(T)
    @test !device.metric_caches.cell.valid && !device.metric_caches.face.valid
    @test host_data(device.node_coordinates) == cpu.node_coordinates
    @test host_data(cell_metrics(device)) == cell_metrics(cpu)
    @test host_data(face_metrics(device)) == face_metrics(cpu)
    refresh_cell_metrics!(device)
    refresh_face_metrics!(device)
    @test host_data(cell_metrics(device)) == cell_metrics(cpu)
    @test host_data(face_metrics(device)) == face_metrics(cpu)

    roundtrip = Adapt.adapt(Array,device)
    @test roundtrip.backend isa CPU
    @test roundtrip.state !== device.state
    @test roundtrip.node_coordinates == cpu.node_coordinates
    @test cell_metrics(roundtrip) == cell_metrics(cpu)
    @test face_metrics(roundtrip) == face_metrics(cpu)
  end
end

@testset "Disabled discrete metric storage on devices" begin
  cpu = DiscreteGrid(collect(0.0:4.0),2;compute_metrics=false,cache_mode=:off)
  device = Adapt.adapt(DiscreteTestAdaptor(),cpu)
  @test device.metric_caches === nothing
  update!(device,1.0,(;))
  @test Array(first(device.node_coordinates)) == first(cpu.node_coordinates)
  @test_throws ArgumentError cell_metrics(device)
  @test_throws ArgumentError face_metrics(device)
end
