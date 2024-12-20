module SymplecticFactorizations

import LinearAlgebra
using LinearAlgebra: mul!, Diagonal, qr, Factorization, svd, require_one_based_indexing

export 
    # symplectic stuff
    issymplectic, symplecticform, BlockForm, PairForm, randsymplectic,
    # polar decomposition
    polar, Polar


include("symplectic.jl")

include("polar.jl")

end
