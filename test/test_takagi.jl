@testitem "symplectic features" begin
    using SymplecticFactorizations
    using LinearAlgebra: Diagonal, adjoint, Symmetric

    @testset "random objects" begin

        n = rand(1:20)
        R_sym = Symmetric(rand(Float64, n, n))
        C_sym = Symmetric(rand(ComplexF64, n, n))

        Fr_sym = takagi(R_sym)
        Qr_sym, Sr_sym = takagi(R_sym)
        @test Fr_sym.Q == Qr_sym && Fr_sym.S == Sr_sym
        @test isapprox(Qr_sym * Diagonal(Sr_sym) * transpose(Qr_sym), Matrix{ComplexF64}(R_sym), atol=1e-5)
        
        #= error on LinearAlgebra.jl end
        @test Fc_sym = takagi(C_sym)
        @test Qc_sym, Sc_sym = takagi(C_sym)
        @test Fc_sym.Q == Qc_sym && Fc_sym.S == Sc_sym
        @testisapprox(Qc_sym * Diagonal(Sc_sym) * transpose(Qc_sym), Matrix(C_sym), atol=1e-5)
        =#
    end
end