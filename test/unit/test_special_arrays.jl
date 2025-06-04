using Test, Random, KernelAbstractions, CUDA

include("src/special_arrays.jl")

using .SpecialArrays 

@testset "Rectlinear2D 1 fixed axis" begin
    slice = rand(-0.9:0.01:0.9, 4)
    A = reshape(repeat(slice, 5, 1), 4, 5)
    B = RectlinearArray(A, (2,))

    # Test to see if the array looks like it should (this also tests the getindex function)
    @test B[1,1] == B[1,2] == B[1,3] == B[1,4] == B[1,5]
    @test B[2,1] == B[2,2] == B[2,3] == B[2,4] == B[2,5]
    @test B[3,1] == B[3,2] == B[3,3] == B[3,4] == B[3,5]
    @test B[4,1] == B[4,2] == B[4,3] == B[4,4] == B[4,5]

    # Test to see if the array "looks" like it should (in terms of dimensionality)
    @test size(A) == size(B)
    @test ndims(A) == ndims(B)

    # Test to see if the the array stores less data than A 
    @test Base.summarysize(B.data) ≤ Base.summarysize(A)
    @test Base.summarysize(B) ≤ Base.summarysize(A)

    # Test to see if the array can be modified with setindex! (note: the array will never have 3, so this is a valid test)
    B[1,4] = 3
    B[3,5] = 4
    @test B[1,1] == B[1,2] == B[1,3] == B[1,4] == B[1,5] == 3
    @test B[3,1] == B[3,2] == B[3,3] == B[3,4] == B[3,5] == 4

    # Test similar function
    C = similar(B)
    @test typeof(C) == typeof(B)
    @test typeof(C.data) == typeof(B.data)
    @test eltype(C) == eltype(B)

    # Test broadcasting (this also tests the Base.similar implementation necessary for broadcasting)
    C[1,1] = C[2,1] = C[3,1] = C[4,1] = 2
    C .+= 10
    @test C[1,1] == C[1,2] == C[1,3] == C[1,4] == C[1,5] == 12
    @test C[2,1] == C[2,2] == C[2,3] == C[2,4] == C[2,5] == 12
    @test C[3,1] == C[3,2] == C[3,3] == C[3,4] == C[3,5] == 12
    @test C[4,1] == C[4,2] == C[4,3] == C[4,4] == C[4,5] == 12
end 

@testset "Rectlinear3D 1/2 fixed axes" begin
    # --- Begin 1 fixed axis --- #
    slice = rand(-0.9:0.01:0.9, 4, 5)
    A = repeat(slice, 1, 1, 3)
    B = RectlinearArray(A, (3,))

    # Test to see if the array looks like it should (this also tests the getindex function)
    @test B[1,1,1] == B[1,1,2] == B[1,1,3]
    @test B[2,1,1] == B[2,1,2] == B[2,1,3]
    @test B[3,1,1] == B[3,1,2] == B[3,1,3]

    @test B[1,3,1] == B[1,3,2] == B[1,3,3]
    @test B[2,3,1] == B[2,3,2] == B[2,3,3]
    @test B[3,3,1] == B[3,3,2] == B[3,3,3]

    # Test to see if the array "looks" like it should (in terms of dimensionality)
    @test size(A) == size(B)
    @test ndims(A) == ndims(B)

    # Test to see if the the array stores less data than A 
    @test Base.summarysize(B.data) ≤ Base.summarysize(A)
    @test Base.summarysize(B) ≤ Base.summarysize(A)

    # Test to see if the array can be modified with setindex! (note: the array will never have 3, so this is a valid test)
    B[1,4,1] = 3
    B[3,5,1] = 4
    @test B[1,4,1] == B[1,4,2] == B[1,4,3] == 3
    @test B[3,5,1] == B[3,5,2] == B[3,5,3] == 4

    # Test similar function
    B_sim = similar(B)
    @test typeof(B_sim) == typeof(B)
    @test typeof(B_sim.data) == typeof(B.data)
    @test eltype(B_sim) == eltype(B)

    # Test broadcasting (this also tests the Base.similar implementation necessary for broadcasting)
    B_zeros = RectlinearArray(zeros(Float64, (4,5,3)), (3,))
    B_zeros .+= 10
    @test all(i -> i == 10, B_zeros)
    # --- End 1 fixed axis --- #
    
    # --- Begin 2 fixed axes --- #
    slice = rand(-0.9:0.01:0.9, 4)
    C = repeat(slice, 1, 5, 3)
    D = RectlinearArray(C, (2,3))

    # Test to see if the array looks like it should (this also tests the getindex function)
    @test D[1,1,1] == D[1,1,2] == D[1,1,3] == D[1,3,1] == D[1,3,2] == D[1,3,3]
    @test D[2,1,1] == D[2,1,2] == D[2,1,3] == D[2,3,1] == D[2,3,2] == D[2,3,3]
    @test D[3,1,1] == D[3,1,2] == D[3,1,3] == D[3,3,1] == D[3,3,2] == D[3,3,3]

    # Test to see if the array "looks" like it should (in terms of dimensionality)
    @test size(C) == size(D)
    @test ndims(C) == ndims(D)

    # Test to see if the the array stores less data than A 
    @test Base.summarysize(D.data) ≤ Base.summarysize(C)
    @test Base.summarysize(D) ≤ Base.summarysize(C)

    # Test to see if the array can be modified with setindex! (note: the array will never have 3, so this is a valid test)
    D[1,4,1] = 3
    D[3,5,1] = 4
    @test D[1,4,1] == D[1,4,2] == D[1,4,3] == D[1,2,1] == D[1,2,2] == D[1,2,3] == 3
    @test D[3,5,1] == D[3,5,2] == D[3,5,3] == D[3,2,1] == D[3,2,2] == D[3,2,3] == 4

    # Test similar function
    D_sim = similar(D)
    @test typeof(D_sim) == typeof(D)
    @test typeof(D_sim.data) == typeof(D.data)
    @test eltype(D_sim) == eltype(D)

    # Test broadcasting (this also tests the Base.similar implementation necessary for broadcasting)
    D_zeros = RectlinearArray(zeros(Float64, (4,5,3)), (2,3))
    D_zeros .+= 10
    @test all(i -> i == 10, D_zeros)
    # --- End 2 fixed axes --- #
end 

@testset "Uniform1/2/3D" begin
    # --- Begin 1D --- #
    A = fill(5.0, (10,))
    B = RectlinearArray(A, (1,))

    # Test to see if the array looks like it should (this also tests the getindex function)
    @test all(i -> i == 5.0, B)

    # Test to see if the array "looks" like it should (in terms of dimensionality)
    @test size(A) == size(B)
    @test ndims(A) == ndims(B)

    # Test to see if the the array stores less data than A 
    @test Base.summarysize(B.data) ≤ Base.summarysize(A)
    @test Base.summarysize(B) ≤ Base.summarysize(A)

    # Test to see if the array can be modified with setindex!
    B[1] = 2.0
    @test all(i -> i == 2.0, B)

    # Test similar function
    B_sim = similar(B)
    @test typeof(B_sim) == typeof(B)
    @test typeof(B_sim.data) == typeof(B.data)
    @test eltype(B_sim) == eltype(B)

    # Test broadcasting (this also tests the Base.similar implementation necessary for broadcasting)
    B .+= 10
    @test all(i -> i == 12.0, B)
    # --- End 1D --- #
    
    # --- Begin 2D --- #
    C = fill(5.0, (4,5))
    D = RectlinearArray(C, (1,2))

    # Test to see if the array looks like it should (this also tests the getindex function)
    @test all(i -> i == 5.0, D)

    # Test to see if the array "looks" like it should (in terms of dimensionality)
    @test size(C) == size(D)
    @test ndims(C) == ndims(D)

    # Test to see if the the array stores less data than A 
    @test Base.summarysize(D.data) ≤ Base.summarysize(C)
    @test Base.summarysize(D) ≤ Base.summarysize(C)

    # Test to see if the array can be modified with setindex! (note: the array will never have 3, so this is a valid test)
    D[3,4] = 2.0
    @test all(i -> i == 2.0, D)

    # Test similar function
    D_sim = similar(D)
    @test typeof(D_sim) == typeof(D)
    @test typeof(D_sim.data) == typeof(D.data)
    @test eltype(D_sim) == eltype(D)

    # Test broadcasting (this also tests the Base.similar implementation necessary for broadcasting)
    D .+= 10
    @test all(i -> i == 12.0, D)
    # --- End 2D --- #
    
    # --- Begin 3D --- #
    E = fill(5.0, (4,5,6))
    F = RectlinearArray(E, (1,2,3))

    # Test to see if the array looks like it should (this also tests the getindex function)
    @test all(i -> i == 5.0, F)

    # Test to see if the array "looks" like it should (in terms of dimensionality)
    @test size(E) == size(F)
    @test ndims(E) == ndims(F)

    # Test to see if the the array stores less data than A 
    @test Base.summarysize(F.data) ≤ Base.summarysize(E)
    @test Base.summarysize(F) ≤ Base.summarysize(E)

    # Test to see if the array can be modified with setindex! (note: the array will never have 3, so this is a valid test)
    F[3,4,5] = 2.0
    @test all(i -> i == 2.0, F)

    # Test similar function
    F_sim = similar(F)
    @test typeof(F_sim) == typeof(F)
    @test typeof(F_sim.data) == typeof(F.data)
    @test eltype(F_sim) == eltype(F)

    # Test broadcasting (this also tests the Base.similar implementation necessary for broadcasting)
    F .+= 10
    @test all(i -> i == 12.0, F)
    # --- End 3D --- #
end 

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

    @test isapprox(Array(output), A * Array(transpose(B)))

end
