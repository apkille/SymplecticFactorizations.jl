@testitem "polar decomposition" begin
    using SymplecticFactorizations
    using LinearAlgebra: adjoint, isposdef, eigvals

    @testset "random objects" begin

        n = rand(1:5)
        J = BlockForm(n)
        Omega = PairForm(n)
        S_block = randsymplectic(J)
        S_pair = randsymplectic(Omega)

        F_block = polar(S_block)
        O_block, P_block = polar(S_block)
        @test F_block.O == O_block && F_block.P == P_block
        F_pair = polar(S_pair)
        O_pair, P_pair = polar(S_pair)
        @test F_pair.O == O_pair && F_pair.P == P_pair 
        @test issymplectic(J, O_block, atol = 1e-5) && issymplectic(J, P_block, atol = 1e-5) 
        @test issymplectic(Omega, O_pair, atol = 1e-5) && issymplectic(Omega, P_pair, atol = 1e-5)
        @test isapprox(inv(O_pair), transpose(O_pair), atol = 1e-5) && isapprox(inv(O_block), transpose(O_block), atol = 1e-5)
        @test isapprox(P_block, transpose(P_block), atol = 1e-5) && all(i > 0 for i in eigvals(P_block))
        @test isapprox(P_pair, transpose(P_pair), atol = 1e-5) && all(i > 0 for i in eigvals(P_pair))
        @test isapprox(O_block * P_block, S_block, atol = 1e-5) && isapprox(O_pair * P_pair, S_pair, atol = 1e-5)
    end

    @testset "symplectic type" begin
        n = rand(1:5)
        J = BlockForm(n)
        Omega = PairForm(n)
        S_block = randsymplectic(Symplectic, J)
        S_pair = randsymplectic(Symplectic, Omega)

        F_block = polar(S_block)
        O_block, P_block = polar(S_block)
        @test F_block.O == O_block && F_block.P == P_block
        F_pair = polar(S_pair)
        O_pair, P_pair = polar(S_pair)
        @test F_pair.O == O_pair && F_pair.P == P_pair
        @test O_block isa Symplectic && O_pair isa Symplectic
        @test P_block isa Symplectic && P_pair isa Symplectic
        @test issymplectic(O_block, atol = 1e-5) && issymplectic(P_block, atol = 1e-5) 
        @test issymplectic(O_pair, atol = 1e-5) && issymplectic(P_pair, atol = 1e-5)
        @test isapprox(inv(O_pair), transpose(O_pair), atol = 1e-5) && isapprox(inv(O_block), transpose(O_block), atol = 1e-5)
        @test isapprox(P_block, transpose(P_block), atol = 1e-5) && all(i > 0 for i in eigvals(P_block.data))
        @test isapprox(P_pair, transpose(P_pair), atol = 1e-5) && all(i > 0 for i in eigvals(P_pair.data))
        @test isapprox(O_block * P_block, S_block, atol = 1e-5) && isapprox(O_pair * P_pair, S_pair, atol = 1e-5)
    end
end