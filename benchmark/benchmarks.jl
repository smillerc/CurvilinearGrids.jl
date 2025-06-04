using BenchmarkTools, Random

include("src/RectlinearArrays.jl")

using .RectlinearArrays

const SUITE = BenchmarkGroup()

A = rand(Float64, (10,10,10))

SUITE["construct"] = BenchmarkGroup()

SUITE["construct"]["3d-1f"] = @benchmarkable begin
    B = RectlinearArray($A, (1,))
end
SUITE["construct"]["3d-2f"] = @benchmarkable begin
    B = RectlinearArray($A, (1,2))
end
SUITE["construct"]["3d-3f"] = @benchmarkable begin
    B = RectlinearArray($A, (1,2,3))
end

B₁ = RectlinearArray(A, (1,))
B₂ = RectlinearArray(A, (1,2))
B₃ = RectlinearArray(A, (1,2,3))

SUITE["getindex"] = BenchmarkGroup()

SUITE["getindex"]["3d-1f"] = @benchmarkable begin
    b = $B₁[1,2,3]
end
SUITE["getindex"]["3d-2f"] = @benchmarkable begin
    b = $B₂[1,2,3]
end
SUITE["getindex"]["3d-3f"] = @benchmarkable begin
    b = $B₃[1,2,3]
end

SUITE["setindex!"] = BenchmarkGroup()

SUITE["setindex!"]["3d-1f"] = @benchmarkable begin
    $B₁[1,2,3] = 1
end
SUITE["setindex!"]["3d-2f"] = @benchmarkable begin
    $B₂[1,2,3] = 1
end
SUITE["setindex!"]["3d-3f"] = @benchmarkable begin
    $B₃[1,2,3] = 1
end
