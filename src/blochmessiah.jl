"""
    BlochMessiah <: Factorization

Matrix factorization type of the Bloch-Messiah/Euler decomposition of a symplectic matrix `S`.
This is the return type of [`blochmessiah(_)`](@ref), the corresponding matrix factorization function.

If `F::BlochMessiah` is the factorization object, `O`, `values` and Q` can be obtained
via `F.O`, `F.values`, and `F.Q`.

Iterating the decomposition produces the components `O`, `values`, and `Q`, in that order.

# Examples
```
julia> S = Symplectic(BlockForm(1), [1. 1.; 0. 1.])
2×2 Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}:
 1.0  1.0
 0.0  1.0

julia> issymplectic(S)
true

julia> F = blochmessiah(S)
BlochMessiah{Float64, Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}, Vector{Float64}}
O factor:
2×2 Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}:
 0.850651  -0.525731
 0.525731   0.850651
values:
1-element Vector{Float64}:
 1.618033988749895
Q factor:
2×2 Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}:
  0.525731  0.850651
 -0.850651  0.525731

julia> isapprox(S, F.O * Diagonal(vcat(F.values, F.values .^ (-1))) * F.Q, atol = 1e-10)
true

julia> issymplectic(F.O, atol = 1e-10) && issymplectic(F.Q, atol = 1e-10)
true

julia> O, values, Q = F; # destructuring via iteration

julia> O == F.O && values == F.values && Q == F.Q
true
```
"""
struct BlochMessiah{T,M<:AbstractArray{T},N<:AbstractVector{T}} <: Factorization{T}
    O::M
    values::N
    Q::M
    function BlochMessiah{T,M,N}(O,values,Q) where {T,M<:AbstractArray{T},N<:AbstractVector{T}}
        require_one_based_indexing(O, values, Q)
        new{T,M,N}(O, values, Q)
    end
end
BlochMessiah{T}(O::AbstractArray{T}, values::AbstractVector{T}, Q::AbstractArray{T}) where {T} = BlochMessiah{T,typeof(O),typeof(values)}(O,values,Q)

# iteration for destructuring into components
Base.iterate(F::BlochMessiah) = (F.O, Val(:values))
Base.iterate(F::BlochMessiah, ::Val{:values}) = (F.values, Val(:Q))
Base.iterate(F::BlochMessiah, ::Val{:Q}) = (F.Q, Val(:done))
Base.iterate(F::BlochMessiah, ::Val{:done}) = nothing

"""
    blochmessiah(form::SymplecticForm, S::AbstractMatrix) -> BlochMessiah
    blochmessiah(S::Symplectic) -> BlochMessiah

Compute the Bloch-Messiah/Euler decomposition of a symplectic matrix `S` and return a `BlockMessiah` object.

The orthogonal symplectic matrices `O` and `Q` as well as the singular values `values` can be obtained
via `F.O`, `F.Q`, and `F.values`, respectively.

Iterating the decomposition produces the components `O`, `values`, and `Q`, in that order.

# Examples
```
julia> S = Symplectic(BlockForm(1), [1. 1.; 0. 1.])
2×2 Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}:
 1.0  1.0
 0.0  1.0

julia> issymplectic(S)
true

julia> F = blochmessiah(S)
BlochMessiah{Float64, Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}, Vector{Float64}}
O factor:
2×2 Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}:
 0.850651  -0.525731
 0.525731   0.850651
values:
1-element Vector{Float64}:
 1.618033988749895
Q factor:
2×2 Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}:
  0.525731  0.850651
 -0.850651  0.525731

julia> isapprox(S, F.O * Diagonal(vcat(F.values, F.values .^ (-1))) * F.Q, atol = 1e-10)
true

julia> issymplectic(F.O, atol = 1e-10) && issymplectic(F.Q, atol = 1e-10)
true

julia> O, values, Q = F; # destructuring via iteration

julia> O == F.O && values == F.values && Q == F.Q
true
```
"""
function blochmessiah(x::Symplectic{F,T,D}) where {F<:SymplecticForm,T,D<:AbstractMatrix{T}} 
    form = x.form
    O, values, Q = _blochmessiah(form, x.data)
    return BlochMessiah{T}(Symplectic(form, O), values, Symplectic(form, Q))
end
function blochmessiah(form::F, x::AbstractMatrix{T}) where {F<:SymplecticForm,T<:Real}
    O, values, Q = _blochmessiah(form, x)
    return BlochMessiah{T}(O, values, Q)
end
function _blochmessiah(form::BlockForm, x::AbstractMatrix{T}) where {T<:Real}
    O, P = polar(x)
    n = form.n
    vals, vecs = eigen(Symmetric(P), sortby = x -> isless(1.0, x) ? -1/x : 1/x)
    @inbounds for i in Base.OneTo(n)
        vecs[1, i] < 0.0 && (@view(vecs[:, i]) .*= -1.0)
        vecs[1, i+n] < 0.0 && (@view(vecs[:, i+n]) .*= -1.0)
    end
    O′ = O * vecs
    Q′ = vecs'
    values′ = vals[1:n]
    return BlochMessiah{T}(O′, values′, Q′)
end
function _blochmessiah(form::PairForm, x::AbstractMatrix{T}) where {T<:Real}
    O, P = polar(x)
    n = form.n
    vals, vecs = eigen(Symmetric(P), sortby = x -> isless(1.0, x) ? -1/x : 1/x)
    Q′ = P
    @inbounds for i in Base.OneTo(n)
        vecs[1, i] < 0.0 && (@view(vecs[:, i]) .*= -1.0)
        vecs[1, i+n] < 0.0 && (@view(vecs[:, i+n]) .*= -1.0)
        Q′[2i-1, :] .= @view(vecs[:, i])
        Q′[2i, :] .= @view(vecs[:, i+n])
    end
    O′ = O * Q′'
    values′ = vals[1:n]
    return BlochMessiah{T}(O′, values′, Q′)
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, F::BlochMessiah{<:Any,<:AbstractArray,<:AbstractVector})
    Base.summary(io, F); println(io)
    println(io, "O factor:")
    Base.show(io, mime, F.O)
    println(io, "\nvalues:")
    Base.show(io, mime, F.values)
    println(io, "\nQ factor:")
    Base.show(io, mime, F.Q)
end