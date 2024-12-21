module SymplecticFactorizations

import LinearAlgebra
using LinearAlgebra: mul!, Diagonal, qr, Factorization, svd, require_one_based_indexing, Symmetric, eigen

export 
    # symplectic stuff
    issymplectic, symplecticform, BlockForm, PairForm, randsymplectic,
    # polar decomposition
    polar, Polar,
    # takagi/autonne decomposition
    takagi, Takagi

include("symplectic.jl")

include("polar.jl")

include("takagi.jl")

end
