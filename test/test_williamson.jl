@testitem "williamson decomposition" begin
    using SymplecticFactorizations
    using LinearAlgebra: adjoint, eigvals, transpose, Diagonal

    @testset "random objects" begin

        rn = rand(2:10)
        n = 2 * rn
        J = BlockForm(rn)
        Omega = PairForm(rn)
        X = rand(n, n)
        V = X' * X

        F_block = williamson(J, V)
        S_block, spectrum_block = williamson(J, V)
        @test F_block.S == S_block && F_block.spectrum == spectrum_block
        F_pair = williamson(Omega, V)
        S_pair, spectrum_pair = williamson(Omega, V)
        @test F_pair.S == S_pair && F_pair.spectrum == spectrum_pair
        @test issymplectic(Omega, S_pair, atol = 1e-5) && issymplectic(J, S_block, atol = 1e-5)
        @test isapprox(S_pair * V * transpose(S_pair), Matrix(Diagonal(repeat(spectrum_pair, inner = 2))), atol = 1e-5)
        @test isapprox(S_block * V * transpose(S_block), Matrix(Diagonal(repeat(spectrum_block, 2))), atol = 1e-5)
    end
end