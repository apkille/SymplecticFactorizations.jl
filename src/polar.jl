"""
    Polar <: Factorization

Matrix factorization type of the polar decomposition of a symplectic matrix `S`.
This is the return type of [`polar(_)`](@ref), the corresponding matrix factorization function.

If `F::Polar` is the factorization object, `O` and `P` can be obtained
via `F.O` and `F.P`, such that `S = O * P`.

Iterating the decomposition produces the components `O` and `P`.

# Examples
```jldoctest
julia> S = [1. 1.; 0. 1.]
2×2 Matrix{Float64}:
 1.0  1.0
 0.0  1.0

julia> issymplectic(BlockForm(1), S)
true

julia> F = polar(S)
Polar{Float64, Matrix{Float64}, Matrix{Float64}}
O factor:
2×2 Matrix{Float64}:
  0.894427  0.447214
 -0.447214  0.894427
P factor:
2×2 Matrix{Float64}:
 0.894427  0.447214
 0.447214  1.34164

julia> isapprox(F.O * F.P, S)
true

julia> O, P = F; # destructuring via iteration

julia> O == F.O && P == F.P
true
```
"""
struct Polar{T,M<:AbstractArray{T},N<:AbstractArray{T}} <: Factorization{T}
    O::M
    P::N
    function Polar{T,M,N}(O, P) where {T,M<:AbstractArray{T},N<:AbstractArray{T}}
        require_one_based_indexing(O, P)
        new{T,M,N}(O, P)
    end
end
Polar{T}(O::AbstractArray{T}, P::AbstractArray{T}) where {T} = Polar{T,typeof(O),typeof(P)}(O,P)

# iteration for destructuring into components
Base.iterate(F::Polar) = (F.O, Val(:P))
Base.iterate(F::Polar, ::Val{:P}) = (F.P, Val(:done))
Base.iterate(F::Polar, ::Val{:done}) = nothing

"""
    polar(S::AbstractMatrix) -> Polar
    polar(S::Symplectic) -> Polar

Compute the polar decomposition of a symplectic matrix `S` and return a `Polar` object.

`O` and `P` can be obtained from the factorization `F` via `F.O` and `F.P`, such that `S = O * P`.
For the symplectic polar decomposition case, `O` is an orthogonal symplectic matrix and `P` is a positive-definite
symmetric symplectic matrix.

Iterating the decomposition produces the components `O` and `P`.

# Examples
```jldoctest
julia> S = [1. 1.; 0. 1.]
2×2 Matrix{Float64}:
 1.0  1.0
 0.0  1.0

julia> issymplectic(BlockForm(1), S)
true

julia> F = polar(S)
Polar{Float64, Matrix{Float64}, Matrix{Float64}}
O factor:
2×2 Matrix{Float64}:
  0.894427  0.447214
 -0.447214  0.894427
P factor:
2×2 Matrix{Float64}:
 0.894427  0.447214
 0.447214  1.34164

julia> isapprox(F.O * F.P, S)
true

julia> O, P = F; # destructuring via iteration

julia> O == F.O && P == F.P
true
```
"""
function polar(x::AbstractMatrix{T}) where {T}
    O, P = _polar(x)
    return Polar{T}(O, P)
end
function polar(x::Symplectic{F,T,D}) where {F<:SymplecticForm,T,D<:AbstractMatrix{T}} 
    O, P = _polar(x.data)
    return Polar{T}(Symplectic(x.form, O), Symplectic(x.form, P))
end
function _polar(x::AbstractMatrix{T}) where {T}
    fact = svd(x)
    dims = size(x)
    O, P = zeros(T, dims), zeros(T, dims)
    mul!(O, fact.U, fact.Vt)
    copyto!(P, fact.V * Diagonal(fact.S) * fact.Vt)
    return O, P
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, F::Polar{<:Any,<:AbstractArray,<:AbstractArray})
    Base.summary(io, F); println(io)
    println(io, "O factor:")
    Base.show(io, mime, F.O)
    println(io, "\nP factor:")
    Base.show(io, mime, F.P)
end