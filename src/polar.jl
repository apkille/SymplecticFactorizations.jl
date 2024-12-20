struct Polar{T,M<:AbstractArray{T},N<:AbstractArray{T}} <: Factorization{T}
    P::M
    O::N
    function Polar{T,M,N}(P, O) where {T,M<:AbstractArray{T},N<:AbstractArray{T}}
        require_one_based_indexing(P, O)
        new{T,M,N}(P, O)
    end
end
Polar{T}(P::AbstractArray{T}, O::AbstractArray{T}) where {T} = Polar{T,typeof(P),typeof(O)}(P,O)

# iteration for destructuring into components
Base.iterate(F::Polar) = (F.P, Val(:O))
Base.iterate(F::Polar, ::Val{:O}) = (F.O, Val(:done))
Base.iterate(F::Polar, ::Val{:done}) = nothing

function polar(x::AbstractMatrix{T}) where {T}
    fact = svd(x, full = false)
    dims = size(x)
    P, O = zeros(T, dims), zeros(T, dims)
    copyto!(P, fact.V * Diagonal(fact.S) * fact.Vt)
    mul!(O, fact.U, fact.Vt)
    return Polar{T}(P, O)
end