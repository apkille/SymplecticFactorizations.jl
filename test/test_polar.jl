@testitem "symplectic features" begin
    using SymplecticFactorizations
    using LinearAlgebra: adjoint, isposdef, eigvals

    @testset "random objects" begin

        n = rand(1:20)
        J = BlockForm(n)
        Omega = PairForm(n)
        S_block = randsymplectic(J)
        S_pair = randsymplectic(Omega)

        F_block = polar(S_block)
        P_block, O_block = polar(S_block)
        @test F_block.P == P_block && F_block.O == O_block
        F_pair = polar(S_pair)
        P_pair, O_pair = polar(S_pair)
        @test F_pair.P == P_pair && F_pair.O == O_pair
        @test issymplectic(J, P_block, atol = 1e-5) && issymplectic(J, O_block, atol = 1e-5)
        @test issymplectic(Omega, P_pair, atol = 1e-5) && issymplectic(Omega, O_pair, atol = 1e-5)
        @test isapprox(P_block, transpose(P_block), atol = 1e-5) && all(i > 0 for i in eigvals(P_block))
        @test isapprox(P_pair, transpose(P_pair), atol = 1e-5) && all(i > 0 for i in eigvals(P_pair))
        @test isapprox(inv(O_pair), transpose(O_pair), atol = 1e-5) && isapprox(inv(O_block), transpose(O_block), atol = 1e-5)
        
    end
end