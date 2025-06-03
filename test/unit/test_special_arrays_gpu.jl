using KernelAbstractions, CUDA, Random, Test

include("src/special_arrays_gpu.jl")

using .SpecialArrays 

@testset "RectlinearArray2D GPU MatMul" begin
    # Simple kernel for matrix multiplication
    @kernel function matmul_kernel!(output, a, b)
        i, j = @index(Global, NTuple)
    
        tmp_sum = zero(eltype(output))
        for k in 1:size(a)[2]
            tmp_sum += a[i, k] * b[k, j]
        end
    
        output[i, j] = tmp_sum
    end
    
    # Creating a wrapper kernel for launching with error checks
    function matmul!(output, a, b)
        if size(a)[2] != size(b)[1]
            println("Matrix size mismatch!")
            return nothing
        end
        kernel! = matmul_kernel!(CUDABackend())
        kernel!(output, a, b, ndrange = size(output))
        return
    end

    A = repeat(rand(Float32, 4), 1, 5)
    Aᵪ = CuArray(A)
    Aᵣ = RectlinearArray(Aᵪ, (2,))

    B = repeat(rand(Float32, 3), 1, 5)
    Bᵪ = CuArray(Array(transpose(B)))
    Bᵣ = RectlinearArray(Bᵪ, (1,))

    output = CUDA.zeros(Float32, (4,3))

    matmul!(output, Aᵣ, Bᵣ)
    KernelAbstractions.synchronize(CUDABackend())

    @test isapprox(output, a*b)

end
