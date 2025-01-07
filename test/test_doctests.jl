@testitem "Doctests" tags=[:doctests] begin
    using Documenter
    using SymplecticFactorizations

    DocMeta.setdocmeta!(SymplecticFactorizations, :DocTestSetup, :(using SymplecticFactorizations, LinearAlgebra); recursive=true)
    doctest(SymplecticFactorizations)
end