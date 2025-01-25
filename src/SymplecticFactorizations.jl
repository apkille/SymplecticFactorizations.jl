module SymplecticFactorizations

import LinearAlgebra
using LinearAlgebra: mul!, Diagonal, qr, Factorization, svd, require_one_based_indexing, Symmetric, eigen, eigen!, I, eigvals, adjoint, eigvecs, normalize!

export 
    # symplectic stuff
    Symplectic, issymplectic, symplecticform, BlockForm, PairForm, randsymplectic,
    # polar decomposition
    polar, Polar,
    # takagi/autonne decomposition
    takagi, Takagi,
    # williamson decomposition
    williamson, Williamson,
    # bloch-messiah/euler decomposition
    blochmessiah, BlochMessiah

include("form.jl")

include("symplectic.jl")

include("polar.jl")

include("takagi.jl")

include("williamson.jl")

include("blochmessiah.jl")

end
