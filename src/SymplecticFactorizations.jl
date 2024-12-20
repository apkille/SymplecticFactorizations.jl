module SymplecticFactorizations

import LinearAlgebra
using LinearAlgebra: Diagonal, qr, eigvals, Factorization, svd, require_one_based_indexing

export 
    # symplectic stuff
    issymplectic, symplecticform, BlockForm, PairForm, randsymplectic,
    # polar decomposition
    polar, Polar


include("symplectic.jl")

include("polar.jl")

end
