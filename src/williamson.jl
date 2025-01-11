"""
    Williamson <: Factorization

Matrix factorization type of the williamson decomposition of a positive-definite matrix `V`.
This is the return type of [`williamson(_)`](@ref), the corresponding matrix factorization function.

If `F::Williamson` is the factorization object, `S` and `spectrum` can be obtained
via `F.S` and `F.spectrum`.

Iterating the decomposition produces the components `S` and `spectrum`.

# Examples
```
julia> V = [7. 2.; 2. 1.]
2×2 Matrix{Float64}:
 7.0  2.0
 2.0  1.0

julia> isposdef(V)
true

julia> F = williamson(BlockForm(1), V)
Williamson{Float64, Matrix{Float64}, Vector{Float64}}
S factor:
2×2 Matrix{Float64}:
 0.448828  -1.95959
 0.61311   -0.448828
symplectic spectrum:
1-element Vector{Float64}:
 1.7320508075688772

julia> isapprox(F.S * V * F.S', Diagonal(repeat(F.spectrum, 2)))
true

julia> S, spectrum = F; # destructuring via iteration

julia> S == F.S && spectrum == F.spectrum
true
```
"""
struct Williamson{T,M<:AbstractArray{T},N<:AbstractVector{T}} <: Factorization{T}
    S::M
    spectrum::N
    function Williamson{T,M,N}(S, spectrum) where {T,M<:AbstractArray{T},N<:AbstractVector{T}}
        require_one_based_indexing(S, spectrum)
        new{T,M,N}(S, spectrum)
    end
end
Williamson{T}(S::AbstractArray{T}, spectrum::AbstractVector{T},) where {T} = Williamson{T,typeof(S),typeof(spectrum)}(S, spectrum)

# iteration for destructuring into components
Base.iterate(F::Williamson) = (F.S, Val(:spectrum))
Base.iterate(F::Williamson, ::Val{:spectrum}) = (F.spectrum, Val(:done))
Base.iterate(F::Williamson, ::Val{:done}) = nothing

"""
    williamson(form::SymplecticForm, V::AbstractMatrix) -> Williamson
    williamson(::Symplectic, form::SymplecticForm, V::AbstractMatrix) -> Williamson

Compute the williamson decomposition of a positive-definite matrix `V` and return a `Williamson` object.

A symplectic matrix `S` and symplectic spectrum `spectrum` can be obtained
via `F.S` and `F.spectrum`.

Iterating the decomposition produces the components `S` and `spectrum`.

# Examples
```
julia> V = [7. 2.; 2. 1.]
2×2 Matrix{Float64}:
 7.0  2.0
 2.0  1.0

julia> isposdef(V)
true

julia> F = williamson(BlockForm(1), V)
Williamson{Float64, Matrix{Float64}, Vector{Float64}}
S factor:
2×2 Matrix{Float64}:
 0.448828  -1.95959
 0.61311   -0.448828
symplectic spectrum:
1-element Vector{Float64}:
 1.7320508075688772

julia> isapprox(F.S * V * F.S', Diagonal(repeat(F.spectrum, 2)))
true

julia> S, spectrum = F; # destructuring via iteration

julia> S == F.S && spectrum == F.spectrum
true
```
"""
function williamson(form::F, x::AbstractMatrix{T}) where {F<:SymplecticForm, T<:Real}
    S, spectrum = _williamson(form, x)
    return Williamson{T}(S, spectrum)
end
function williamson(::Type{Symplectic}, form::F, x::AbstractMatrix{T}) where {F<:SymplecticForm, T<:Real}
    S, spectrum = _williamson(form, x)
    return Williamson{T}(Symplectic(form, S), spectrum)
end
function _williamson(form::PairForm, x::AbstractMatrix{T}) where {T<:Real}
    J = symplecticform(form)
    spectrum = filter(i -> i > 0, imag.(eigvals(J * x, sortby = λ -> abs(λ))))
    D = Diagonal(repeat(spectrum, inner = 2))
    sqrtV = Symmetric(x)^(-1//2)
    X = sqrtV * J * sqrtV
    U = eigvecs(X, sortby = λ -> -abs(λ))
    G = zeros(ComplexF64, 2*form.n, 2*form.n)
    @inbounds for i in Base.OneTo(form.n)
        G[2i-1, 2i-1] = -im/sqrt(2)
        G[2i-1, 2i] = im/sqrt(2)
        G[2i, 2i-1] = 1/sqrt(2)
        G[2i, 2i] = 1/sqrt(2)
    end
    R = real(G * U')
    S = sqrt(D) * R * sqrtV
    return S, spectrum
end
function _williamson(form::BlockForm, x::AbstractMatrix{T}) where {T<:Real}
    J = symplecticform(form)
    spectrum = filter(i -> i > 0, imag.(eigvals(J * x, sortby = λ -> abs(λ))))
    D = Diagonal(repeat(spectrum, 2))
    sqrtV = Symmetric(x)^(-1//2)
    X = sqrtV * J * sqrtV
    U = eigvecs(X, sortby = λ -> (-sign(imag(λ)), -abs(λ)))
    G = zeros(ComplexF64, 2*form.n, 2*form.n)
    @inbounds for i in Base.OneTo(form.n)
        G[i, i] = -im/sqrt(2)
        G[i, i+form.n] = im/sqrt(2)
        G[i+form.n, i] = 1/sqrt(2)
        G[i+form.n, i+form.n] = 1/sqrt(2)
    end
    R = real(G * U')
    S = sqrt(D) * R * sqrtV
    return S, spectrum
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, F::Williamson{<:Any,<:AbstractArray,<:AbstractVector})
    Base.summary(io, F); println(io)
    println(io, "S factor:")
    Base.show(io, mime, F.S)
    println(io, "\nsymplectic spectrum:")
    Base.show(io, mime, F.spectrum)
end