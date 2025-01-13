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
    euler(form::SymplecticForm, x::AbstractMatrix) -> Euler
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
    vals, vecs = eigen(Symmetric(P))
    Q′ = P
    @inbounds for k in Base.OneTo(n)
        Q′[k, :] .= @view(vecs[:, n+k])
        Q′[n+k, :] .= @view(vecs[:, n-k+1])
    end  
    O′ = O * transpose(Q′)
    values′ = vals[n+1:2n]
    return BlochMessiah{T}(O′, values′, Q′)
end
function _blochmessiah(form::PairForm, x::AbstractMatrix{T}) where {T<:Real}
    O, P = polar(x)
    n = form.n
    vals, vecs = eigen(P)
    Q′ = P
    @inbounds for k in Base.OneTo(n)
        Q′[2k-1, :] .= @view(vecs[:, n+k])
        Q′[2k, :] .= @view(vecs[:, n-k+1])
    end
    O′ = O * transpose(Q′)
    values′ = vals[n+1:2n]
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