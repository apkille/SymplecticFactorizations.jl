using Revise # for interactive work on docs
push!(LOAD_PATH,"../src/")

using Documenter
using DocumenterCitations
using SymplecticFactorizations

DocMeta.setdocmeta!(SymplecticFactorizations, :DocTestSetup, :(using SymplecticFactorizations); recursive=true)

function main()
    #bib = CitationBibliography(joinpath(@__DIR__,"src/references.bib"), style=:authoryear)

    makedocs(
    #plugins=[bib],
    doctest = false,
    clean = true,
    sitename = "SymplecticFactorizations.jl",
    format = Documenter.HTML(
        assets=["assets/init.js"],
        canonical = "https://apkille.github.io/SymplecticFactorizations.jl"
    ),
    modules = [SymplecticFactorizations],
    checkdocs = :exports,
    warnonly = false,
    authors = "Andrew Kille",
    pages = [
        "SymplecticFactorizations.jl" => "index.md"
    ]
    )

    deploydocs(
        repo = "github.com/apkille/SymplecticFactorizations.jl.git"
    )
end

main()
