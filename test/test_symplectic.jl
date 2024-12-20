@testitem "symplectic features" begin
    using SymplecticFactorizations
    using SymplecticFactorizations: _rand_orthogonal_symplectic, _rand_unitary
    using LinearAlgebra: eigvals, adjoint

    @testset "random objects" begin

        n = rand(1:20)
        J = BlockForm(n)
        Omega = PairForm(n)
        U_block = _rand_unitary(J)
        U_pair = _rand_unitary(Omega)
        @test isapprox(U_pair', inv(U_pair), atol = 1e-4)
        @test isapprox(U_block', inv(U_block), atol = 1e-4)

        O_block = _rand_orthogonal_symplectic(J)
        O_pair = _rand_orthogonal_symplectic(Omega)
        @test isapprox(O_pair', inv(O_pair), atol = 1e-4)
        @test isapprox(O_block', inv(O_block), atol = 1e-4)
        @test issymplectic(J, O_block, atol = 1e-4)
        @test issymplectic(Omega, O_pair, atol = 1e-4)

        S_block = randsymplectic(J)
        S_pair = randsymplectic(Omega)
        @test issymplectic(J, S_block, atol = 1e-4)
        @test issymplectic(Omega, S_pair, atol = 1e-4)
        
    end


end