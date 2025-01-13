@testitem "bloch-messiah decomposition" begin
    using SymplecticFactorizations
    using LinearAlgebra: transpose, Diagonal

    @testset "random objects" begin

        rn = rand(2:5)
        n = 2 * rn
        J = BlockForm(rn)
        Omega = PairForm(rn)
        S_block = randsymplectic(J)
        S_pair = randsymplectic(Omega)

        F_block = blochmessiah(J, S_block)
        O_block, vals_block, Q_block = blochmessiah(J, S_block)
        @test F_block.O == O_block && F_block.values == vals_block && F_block.Q == Q_block
        F_pair = blochmessiah(Omega, S_pair)
        O_pair, vals_pair, Q_pair = blochmessiah(Omega, S_pair)
        @test F_pair.O == O_pair && F_pair.values == vals_pair && F_pair.Q == Q_pair
        @test issymplectic(Omega, O_pair, atol = 1e-5) && issymplectic(J, O_block, atol = 1e-5)
        @test issymplectic(Omega, Q_pair, atol = 1e-5) && issymplectic(J, Q_block, atol = 1e-5)
        @test isapprox(inv(O_pair), transpose(O_pair), atol = 1e-5) && isapprox(inv(O_block), transpose(O_block), atol = 1e-5)
        @test isapprox(inv(Q_pair), transpose(Q_pair), atol = 1e-5) && isapprox(inv(Q_block), transpose(Q_block), atol = 1e-5)

        @test all(x -> x > 0.0, vals_block) && all(x -> x > 0.0, vals_pair)
        D_block = Diagonal(vcat(vals_block, vals_block .^ (-1)))
        D_pair = Diagonal(collect(Iterators.flatten(zip(vals_pair, vals_pair .^ (-1)))))
        @test isapprox(O_block * D_block * Q_block, S_block, atol = 1e-5)
        @test isapprox(O_pair * D_pair * Q_pair, S_pair, atol = 1e-5)
    end

    @testset "symplectic type" begin

        rn = rand(2:5)
        n = 2 * rn
        J = BlockForm(rn)
        Omega = PairForm(rn)
        S_block = randsymplectic(Symplectic, J)
        S_pair = randsymplectic(Symplectic, Omega)

        F_block = blochmessiah(S_block)
        O_block, vals_block, Q_block = blochmessiah(S_block)
        @test F_block.O == O_block && F_block.values == vals_block && F_block.Q == Q_block
        F_pair = blochmessiah(S_pair)
        O_pair, vals_pair, Q_pair = blochmessiah(S_pair)
        @test F_pair.O == O_pair && F_pair.values == vals_pair && F_pair.Q == Q_pair
        @test issymplectic(O_pair, atol = 1e-5) && issymplectic(O_block, atol = 1e-5)
        @test issymplectic(Q_pair, atol = 1e-5) && issymplectic(Q_block, atol = 1e-5)
        @test isapprox(inv(O_pair), transpose(O_pair), atol = 1e-5) && isapprox(inv(O_block), transpose(O_block), atol = 1e-5)
        @test isapprox(inv(Q_pair), transpose(Q_pair), atol = 1e-5) && isapprox(inv(Q_block), transpose(Q_block), atol = 1e-5)

        @test all(x -> x > 0.0, vals_block) && all(x -> x > 0.0, vals_pair)
        D_block = Diagonal(vcat(vals_block, vals_block .^ (-1)))
        D_pair = Diagonal(collect(Iterators.flatten(zip(vals_pair, vals_pair .^ (-1)))))
        @test isapprox(O_block * D_block * Q_block, S_block, atol = 1e-5)
        @test isapprox(O_pair * D_pair * Q_pair, S_pair, atol = 1e-5)
    end
end