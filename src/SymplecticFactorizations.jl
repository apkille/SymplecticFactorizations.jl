module SymplecticFactorizations

import LinearAlgebra
import LinearAlgebra: givens
using LinearAlgebra: mul!, Diagonal, qr, Factorization, svd, require_one_based_indexing, Symmetric, eigen, eigen!, I, eigvals, adjoint, eigvecs, normalize!, AbstractRotation

export 
    # symplectic stuff
    Symplectic, issymplectic, symplecticform, BlockForm, PairForm, randsymplectic,
    # symplectic givens
    givens, SymplecticGivens,
    # symplectic householder
    householder, SymplecticHouseholder,
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

include("givens.jl")

include("householder.jl")

include("polar.jl")

include("takagi.jl")

include("williamson.jl")

include("blochmessiah.jl")

end
