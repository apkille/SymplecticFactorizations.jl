module SymplecticFactorizations

import LinearAlgebra
using LinearAlgebra: mul!, Diagonal, qr, Factorization, svd, require_one_based_indexing, Symmetric, eigen, I, eigvals, adjoint, eigvecs, normalize!

export 
    # symplectic stuff
    issymplectic, symplecticform, BlockForm, PairForm, randsymplectic,
    # polar decomposition
    polar, Polar,
    # takagi/autonne decomposition
    takagi, Takagi,
    # williamson decomposition
    williamson, Williamson

include("symplectic.jl")

include("polar.jl")

include("takagi.jl")

include("euler.jl")

include("williamson.jl")

end
