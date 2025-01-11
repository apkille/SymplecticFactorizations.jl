@testitem "symplectic features" begin
    using SymplecticFactorizations
    using SymplecticFactorizations: _rand_orthogonal_symplectic, _rand_unitary
    using LinearAlgebra
    using Base: size, axes, eltype, getindex, setindex!, similar, Matrix, Array, AbstractMatrix, parent, copy, copyto!

    @testset "symplectic type" begin

        n = rand(1:5)
        J = BlockForm(n)
        Omega = PairForm(n)

        data_block = randsymplectic(J)
        data_pair = randsymplectic(Omega)

        S_block = Symplectic(J, data_block)
        S_pair = Symplectic(Omega, data_pair)

        @test randsymplectic(Symplectic, J) isa Symplectic && randsymplectic(Symplectic, Omega) isa Symplectic
        @test Symplectic(S_block) == S_block && Symplectic(S_pair) == S_pair
        @test eltype(S_block) == eltype(data_block) && eltype(S_pair) == eltype(data_pair)

        @test isequal(S_block, copy(S_block)) && isequal(S_pair, copy(S_pair))
        @test isapprox(S_block, copy(S_block)) && isapprox(S_pair, copy(S_pair))

        @test size(S_block) == size(data_block) && size(S_pair) == size(data_pair)
        @test size(S_block, 1) == size(data_block, 1) && size(S_pair, 1) == size(data_pair, 1)
        @test axes(S_block) == axes(data_block) && axes(S_pair) == axes(data_pair)
        @test S_block[n] == data_block[n] && S_pair[n] == data_pair[n]
        @test S_block[n, n] == data_block[n, n] && S_pair[n, n] == data_pair[n, n]
        @test setindex!(S_block, 0.0, n) == setindex!(data_block, 0.0, n) && setindex!(S_pair, 0.0, n) == setindex!(data_pair, 0.0, n)
        @test setindex!(S_block, 0.0, n, n) == setindex!(data_block, 0.0, n, n)
        @test setindex!(S_pair, 0.0, n, n) == setindex!(data_pair, 0.0, n, n)

        sim_S_block, sim_S_pair = similar(S_block), similar(S_pair)
        sim_S_block_float32, sim_S_pair_float32 = similar(S_block, Float32), similar(S_pair, Float32)
        sim_S_block_dim, sim_S_pair_dim = similar(S_block, (2,2)), similar(S_pair, (2,2))
        sim_S_block_dim_float32, sim_S_pair_dim_float32 = similar(S_block, Float32, (2,2)), similar(S_pair, Float32, (2,2))
        sim_S_block .= 0.0; sim_S_pair .= 0.0; sim_S_block_float32 .= 0.0; sim_S_pair_float32 .= 0.0;
        sim_S_block_dim .= 0.0; sim_S_pair_dim .= 0.0; sim_S_block_dim_float32 .= 0.0; sim_S_pair_dim_float32 .= 0.0;

        @test sim_S_block == Symplectic(J, sim_S_block.data) && sim_S_pair == Symplectic(Omega, sim_S_pair.data)
        @test sim_S_block_float32 == Symplectic(J, sim_S_block_float32.data) && sim_S_pair_float32 == Symplectic(Omega, sim_S_pair_float32.data)
        @test sim_S_block_dim == Symplectic(J, sim_S_block_dim.data) && sim_S_pair_dim == Symplectic(Omega, sim_S_pair_dim.data)
        @test sim_S_block_dim_float32 == Symplectic(J, sim_S_block_dim_float32.data) && sim_S_pair_dim_float32 == Symplectic(Omega, sim_S_pair_dim_float32.data)

        @test Matrix(S_block) == Matrix(data_block) && Matrix(S_pair) == Matrix(data_pair)
        @test Array(S_block) == Matrix(data_block) && Array(S_pair) == Matrix(data_pair)
        @test AbstractMatrix{Float32}(S_block) == Symplectic(J, AbstractMatrix{Float32}(data_block)) && AbstractMatrix{Float32}(S_pair) == Symplectic(Omega, AbstractMatrix{Float32}(data_pair))
        @test parent(S_block) == data_block && parent(S_pair) == data_pair

        zero_S_block, zero_S_pair = zero(S_block), zero(S_pair)
        @test copy(S_block) == S_block && copy(S_pair) == S_pair
        @test copyto!(copy(zero_S_block), S_block) == S_block && copyto!(copy(zero_S_pair), S_pair) == S_pair
        @test copyto!(copy(zero_S_block), data_block) == S_block && copyto!(copy(zero_S_pair), data_pair) == S_pair
        @test copyto!(zeros(2n, 2n), S_block) == data_block && copyto!(zeros(2n, 2n), S_pair) == data_pair
    end

    @testset "random objects" begin

        n = rand(1:5)
        J = BlockForm(n)
        Omega = PairForm(n)
        U_block = _rand_unitary(J)
        U_pair = _rand_unitary(Omega)
        @test isapprox(U_pair', inv(U_pair), atol = 1e-5)
        @test isapprox(U_block', inv(U_block), atol = 1e-5)

        O_block = _rand_orthogonal_symplectic(J)
        O_pair = _rand_orthogonal_symplectic(Omega)
        @test isapprox(O_pair', inv(O_pair), atol = 1e-5)
        @test isapprox(O_block', inv(O_block), atol = 1e-5)
        @test issymplectic(J, O_block, atol = 1e-5)
        @test issymplectic(Omega, O_pair, atol = 1e-5)

        S_block = randsymplectic(J)
        S_pair = randsymplectic(Omega)
        @test issymplectic(J, S_block, atol = 1e-5)
        @test issymplectic(Omega, S_pair, atol = 1e-5)

    end

    @testset "linear algebra" begin

        n = rand(1:5)
        J = BlockForm(n)
        Omega = PairForm(n)

        S_block = randsymplectic(Symplectic, J)
        S_pair = randsymplectic(Symplectic, Omega)

        data_block = S_block.data
        data_pair = S_pair.data

        for f in (:svd!, :lu!, :svdvals!, :lq!, :qr!, :eigvals!, :schur!, :eigen!, :hessenberg!)
            @test @eval(LinearAlgebra.$f)(copy(S_block)) == @eval(LinearAlgebra.$f)(copy(data_block))
            @test @eval(LinearAlgebra.$f)(copy(S_pair)) == @eval(LinearAlgebra.$f)(copy(data_pair))
        end
        for f in (:svd, :lu, :svdvals, :lq, :qr, :eigvals, :schur, :eigvecs, :eigen, :hessenberg)
            @test @eval(LinearAlgebra.$f)(S_block) == @eval(LinearAlgebra.$f)(data_block)
            @test @eval(LinearAlgebra.$f)(S_pair) == @eval(LinearAlgebra.$f)(data_pair)
        end
        for f in (:det, :tr, :inv, :pinv, :logdet)
            @test @eval(LinearAlgebra.$f)(S_block) == @eval(LinearAlgebra.$f)(data_block)
            @test @eval(LinearAlgebra.$f)(S_pair) == @eval(LinearAlgebra.$f)(data_pair)
        end
    end
end