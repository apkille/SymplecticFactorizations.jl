@testitem "Doctests" tags=[:doctests] begin
    using Documenter
    using SymplecticFactorizations

    DocMeta.setdocmeta!(SymplecticFactorizations, :DocTestSetup, :(using SymplecticFactorizations); recursive=true)
    doctest(SymplecticFactorizations)
end